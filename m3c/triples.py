"""M3C N-Triple Generator

Usage:
    m3c generate (-h | --help)
    m3c generate [-x <prev> | --diff=<prev>] <path_to_config>

Options:
    -h --help      Show this message and exit
    -x --diff      See Differential Update.

Differential Update:
    A differential update compares the triples produced by a run with that of
    an older one (<prev>). Two distinct sets of triples are produced:
     - an "Add" set which contains the set of triples that are in the current
       run's triples, but not in the older sets.
     - a "Sub" set which contains the set of triples that are in the older
       set, but not in the current run's.

    The corresponding files are written to add.nt and sub.nt.

Instructions:
    Run the generator where you have access to the Postgres database.
"""

from typing import Dict, IO, Iterable, List, Mapping, Tuple

from datetime import datetime
import getopt
import os
import pathlib
import re
import sys
import traceback
import yaml

import psycopg2

from m3c import config
from m3c import db
from m3c.classes import Dataset
from m3c.classes import Organization
from m3c.classes import Person
from m3c.classes import Photo
from m3c.classes import Project
from m3c.classes import Publication
from m3c.classes import Study
from m3c.classes import Tool
from m3c import prefill
from m3c import tools


def diff(prev_path: str, path: str) -> Tuple[Iterable[str], Iterable[str]]:
    prev = pathlib.Path(prev_path)
    previous = []
    for file in prev.glob("*.nt"):
        if file in ["add.nt", "sub.nt"]:
            continue
        with open(file) as f:
            previous += [line for line in list(f) if line]
    previous.sort()

    curr = pathlib.Path(path)
    current = []
    for file in curr.glob("*.nt"):
        if file in ["add.nt", "sub.nt"]:
            continue
        with open(file) as f:
            current += [line for line in list(f) if line]
    current.sort()

    previous_set = set(previous)
    current_set = set(current)

    add = current_set - previous_set
    sub = previous_set - current_set

    return (add, sub)


def get_organizations(sup_cur):
    print("Gathering Organizations")
    orgs = {}
    sup_cur.execute("""\
                    SELECT id, name, type, parent_id
                    FROM organizations
                    WHERE withheld = FALSE""")
    for row in sup_cur:
        org = Organization(org_id=row[0], name=row[1], type=row[2],
                           parent_id=row[3])
        orgs[org.org_id] = org
    return orgs


def make_organizations(namespace, orgs):
    print("Making Organizations")
    triples = []
    for org in orgs.values():
        triples.extend(org.get_triples(namespace))
    print(f"There are {len(orgs)} organizations.")
    return triples


def get_people(sup_cur: db.Cursor) -> Dict[int, Person]:
    print("Gathering People")
    people: Dict[int, Person] = {}
    records = db.get_people(sup_cur)
    for pid, (first, last, display, email, phone, withheld, overview) in records.items():
        person = Person(person_id=str(pid),
                        first_name=first,
                        last_name=last,
                        display_name=display,
                        email=email,
                        phone=phone,
                        withheld=withheld,
                        overview=overview)
        people[pid] = person
    print(f"There are {len(people)} people.")
    return people


def make_people(namespace, people):
    print("Making People Profiles")
    triples = []
    for person in people.values():
        triples.extend(person.get_triples(namespace))
    print(f"There are {len(people)} people.")
    return triples


def link_people_to_org(namespace: str, sup_cur, people, orgs):
    triples = []
    sup_cur.execute("""\
        SELECT person_id, organization_id
          FROM associations as a,
               organizations as o,
               people as p
         WHERE a.person_id = p.id AND a.organization_id = o.id
           AND p.withheld IS NOT TRUE AND o.withheld IS NOT TRUE
    """)
    for row in sup_cur:
        triples.extend(orgs[row[1]].add_person(namespace, row[0]))
    return triples


def make_photos(namespace: str, photos: list):
    print("Making Photo triples")

    triples = []
    for photo in photos:
        triples.extend(photo.get_triples(namespace))

    print(f"There are {len(photos)} photos.")

    return triples


def get_photos(file_storage_root: str, people):
    photos = []
    for person in people.values():
        photo = Photo(file_storage_root, person.person_id, 'jpg')

        jpg = os.path.join(photo.path(), photo.filename())
        if os.path.isfile(jpg):
            photos.append(photo)
            continue

        photo = Photo(file_storage_root, person.person_id, 'png')
        png = os.path.join(photo.path(), photo.filename())
        if os.path.isfile(png):
            photos.append(photo)
            continue

    return photos


def get_publications(
    sup_cur: db.Cursor,
    permitted_people: Dict[int, Person]
) -> Mapping[str, Publication]:
    print("Gathering publications")

    authorships = db.get_pubmed_authorships(sup_cur)
    pubs = db.get_pubmed_publications(sup_cur)

    publications = {}
    for pmid, xml in pubs.items():
        if pmid not in authorships:
            continue

        try:
            pub = Publication.from_pubmed(xml)
            assert pub and pub.pmid == pmid
            for author in authorships[pmid]:
                if author in permitted_people:
                    pub.add_author(author)
            publications[pmid] = pub
        except Exception:
            traceback.print_exc()
            print(f"Skipping publication {pmid}")

    return publications


def make_publications(namespace: str,
                      pubs: Mapping[str, Publication]
                      ) -> List[str]:
    print("Making Publication triples")
    triples = []
    for pub in pubs.values():
        triples.extend(pub.get_triples(namespace))
    print(f"There are {len(pubs)} publications.")
    return triples


def get_projects(mwb_cur: db.Cursor,
                 sup_cur: db.Cursor,
                 people: Dict[int, Person],
                 orgs: List[Organization]
                 ) -> Mapping[str, Project]:
    print("Gathering Workbench Projects")
    projects = {}
    mwb_cur.execute("""\
        SELECT project_id, project_title, COALESCE(project_type, ''),
               COALESCE(project_summary, ''), COALESCE(doi, ''),
               COALESCE(funding_source, ''),
               last_name, first_name, institute, department, laboratory
          FROM project
    """)
    for row in mwb_cur:
        project = Project(
            project_id=row[0],
            project_title=row[1],
            project_type=row[2],
            summary=row[3],
            doi=row[4],
            funding_source=row[5],
        )
        last_names: str = row[6]
        first_names: str = row[7]
        institutes = row[8]
        departments = row[9]
        labs = row[10]

        institute_list = [inst.strip() for inst in institutes.split(';')]
        try:
            department_list = [dept.strip() for dept in departments.split(';')]
        except AttributeError:
            department_list = []
        try:
            lab_list = [lab.strip() for lab in labs.split(';')]
        except AttributeError:
            lab_list = []
        max_range = len(institute_list)
        if len(department_list) > max_range:
            max_range = len(department_list)
        if len(lab_list) > max_range:
            max_range = len(lab_list)

        for i in range(0, max_range):
            # If there are not enough institutes, default to first
            try:
                sup_cur.execute("""
                    SELECT id, parent_id
                    FROM organizations
                    WHERE name=%s AND type='institute'
                    AND withheld = FALSE
                """, (institute_list[i],))
                try:
                    inst_id = sup_cur.fetchone()[0]
                    project.institutes.append(orgs[inst_id].org_id)
                except TypeError:
                    print("Error: Organization does not exist.")
                    print("Organization for project " + project.project_id)
                    print("Organization name: " + institute_list[i])
                    sys.exit()
            except IndexError:
                sup_cur.execute("""
                    SELECT id, parent_id
                    FROM organizations
                    WHERE name=%s AND type='institute'
                    AND withheld = FALSE
                """, (institute_list[0],))
                inst_id = sup_cur.fetchone()[0]

            # If there are not enough departments, default to first
            if departments:
                try:
                    sup_cur.execute("""
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND type='department'
                        AND withheld = FALSE
                    """, (department_list[i],))
                    try:
                        dept_options = {}
                        for row in sup_cur:
                            dept_options[row[0]] = row[1]
                        for dept_id, parent in dept_options.items():
                            if inst_id == parent:
                                department_id = dept_id
                                org_id = orgs[dept_id].org_id
                                project.departments.append(org_id)
                    except TypeError:
                        print("Error: Organization does not exist.")
                        print("Organization for project " + project.project_id)
                        print("Organization name: " + department_list[i])
                        sys.exit()
                except IndexError:
                    sup_cur.execute("""
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND type='department'
                        AND withheld = FALSE
                    """, (department_list[0],))
                    dept_options = {}
                    for row in sup_cur:
                        dept_options[row[0]] = row[1]
                    for dept_id, parent in dept_options.items():
                        if inst_id == parent:
                            department_id = dept_id
            if labs:
                try:
                    sup_cur.execute("""
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND type='laboratory'
                        AND withheld = FALSE
                    """, (lab_list[i],))
                    try:
                        lab_options = {}
                        for row in sup_cur:
                            lab_options[row[0]] = row[1]
                        for lab_id, parent in lab_options.items():
                            if department_id == parent:
                                project.labs.append(orgs[lab_id].org_id)
                    except TypeError:
                        print("Error: Organization does not exist.")
                        print("Organization for project " + project.project_id)
                        print("Organization name: " + lab_list[i])
                        sys.exit()
                except IndexError:
                    pass

        last_name_list = [ln.strip() for ln in last_names.split(';')]
        first_name_list = [fn.strip() for fn in first_names.split(';')]

        for i in range(0, len(last_name_list)):
            last_name = last_name_list[i]
            first_name = first_name_list[i]
            ids = list(db.get_person(sup_cur, first_name, last_name))
            try:
                person_id = ids[0]
                project.pi.append(people[person_id].person_id)
            except (IndexError, KeyError, TypeError):
                print("Error: Person does not exist.")
                print("PI for project " + project.project_id)
                print("Last name: " + last_name)
                print("First name: " + first_name)
                sys.exit()
        projects[project.project_id] = project
    return projects


def make_projects(namespace, projects: Mapping[str, Project]):
    print("Making Workbench Projects")
    triples = []
    summaries = []
    project_count = 0
    for project in projects.values():
        project_triples, summary_line = project.get_triples(namespace)
        triples.extend(project_triples)
        if summary_line:
            summaries.append(summary_line)
        project_count += 1
    print("There are " + str(project_count) + " projects.")
    return triples, summaries


def is_valid_study(study: Study) -> bool:
    if study.study_id.startswith('ST9'):
        # Studies that start with ST9 are testing studies
        print('Test study found like ST9XXXXX. Skipping...')
        return False
    return True


def get_studies(mwb_cur: db.Cursor,
                sup_cur: db.Cursor,
                people: Dict[int, Person],
                orgs: Dict[int, Organization],
                embargoed: List[str]
                ) -> Dict[str, Study]:
    print("Gathering Workbench Studies")
    studies: Dict[str, Study] = {}
    mwb_cur.execute("""\
        SELECT study.study_id, study.study_title,
               COALESCE(study.study_type, ''),
               COALESCE(study.study_summary, ''), study.submit_date,
               study.project_id, study.last_name, study.first_name,
               study.institute, study.department, study.laboratory
        FROM study, study_status_prod
        WHERE study.study_id = study_status_prod.study_id
          AND study_status_prod.status = 1""")

    for row in mwb_cur:
        submit_date = ""
        if row[4]:
            submit_date = f"{row[4]}T00:00:00"

        study = Study(
            study_id=row[0].replace('\n', ''),
            study_title=row[1].replace('\n', '').replace('"', '\\"'),
            study_type=row[2].replace('\n', ''),
            summary=row[3].replace('\n', '').replace('"', '\\"'),
            submit_date=submit_date,
            project_id=row[5].replace('\n', ''))

        # Exclude embargoed studies.
        if study.study_id in embargoed:
            print(f"Skipping embargoed study: {study.study_id}")
            continue

        # Skip invalid studies
        if not is_valid_study(study):
            continue

        last_names: str = row[6]
        first_names: str = row[7]
        institutes = row[8]
        departments = row[9]
        labs = row[10]

        institute_list = [inst.strip() for inst in institutes.split(';')]
        try:
            department_list = [dept.strip() for dept in departments.split(';')]
        except AttributeError:
            department_list = []
        try:
            lab_list = [lab.strip() for lab in labs.split(';')]
        except AttributeError:
            lab_list = []
        max_range = len(institute_list)
        if len(department_list) > max_range:
            max_range = len(department_list)
        if len(lab_list) > max_range:
            max_range = len(lab_list)

        for i in range(0, max_range):
            # If there are not enough institutes, default to first
            try:
                sup_cur.execute("""
                    SELECT id, parent_id
                    FROM organizations
                    WHERE name=%s AND type='institute'
                    AND withheld = FALSE
                """, (institute_list[i],))
                try:
                    inst_id = sup_cur.fetchone()[0]
                    study.institutes.append(orgs[inst_id].org_id)
                except TypeError:
                    print("Error: Organization does not exist.")
                    print("Organization for study " + study.study_id)
                    print("Organization name: " + institute_list[i])
                    sys.exit()
            except IndexError:
                sup_cur.execute("""
                    SELECT id, parent_id
                    FROM organizations
                    WHERE name=%s AND type='institute'
                    AND withheld = FALSE
                """, (institute_list[0],))
                inst_id = sup_cur.fetchone()[0]

            # If there are not enough departments, default to first
            if departments:
                try:
                    sup_cur.execute("""
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND type='department'
                        AND withheld = FALSE
                    """, (department_list[i],))
                    try:
                        dept_options = {}
                        for row in sup_cur:
                            dept_options[row[0]] = row[1]
                        for dept_id, parent in dept_options.items():
                            if inst_id == parent:
                                department_id = dept_id
                                study.departments.append(orgs[dept_id].org_id)
                    except TypeError:
                        print("Error: Organization does not exist.")
                        print("Organization for study " + study.study_id)
                        print("Organization name: " + department_list[i])
                        sys.exit()
                except IndexError:
                    sup_cur.execute("""
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND type='department'
                        AND withheld = FALSE
                    """, (department_list[0],))
                    dept_options = {}
                    for row in sup_cur:
                        dept_options[row[0]] = row[1]
                    for dept_id, parent in dept_options.items():
                        if inst_id == parent:
                            department_id = dept_id
            if labs:
                try:
                    sup_cur.execute("""
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND type='laboratory'
                        AND withheld = FALSE
                    """, (lab_list[i],))
                    try:
                        lab_options = {}
                        for row in sup_cur:
                            lab_options[row[0]] = row[1]
                        for lab_id, parent in lab_options.items():
                            if department_id == parent:
                                study.labs.append(orgs[lab_id].org_id)
                    except TypeError:
                        print("Error: Organization does not exist.")
                        print("Organization for study " + study.study_id)
                        print("Organization name: " + lab_list[i])
                        sys.exit()
                except IndexError:
                    pass

        last_name_list = [ln.strip() for ln in last_names.split(';')]
        first_name_list = [fn.strip() for fn in first_names.split(';')]

        for i in range(0, len(last_name_list)):
            last_name = last_name_list[i]
            first_name = first_name_list[i]

            ids = list(db.get_person(sup_cur, first_name, last_name))
            try:
                person_id = ids[0]
                study.runner.append(people[person_id].person_id)
            except (IndexError, KeyError, TypeError):
                print("Error: Person does not exist.")
                print("Runner for study " + study.study_id)
                print("Last name: " + last_name + '.')
                print("First name: " + first_name + '.')
                sys.exit()

        studies[study.study_id] = study
    return studies


def make_studies(namespace,
                 studies: Dict[str, Study],
                 projects: Mapping[str, Project]
                 ) -> Tuple[List[str], List[str]]:
    print("Making Workbench Studies")
    triples: List[str] = []
    summaries: List[str] = []
    study_count = 0
    no_proj_study = 0
    for study in studies.values():
        if study.project_id not in projects.keys():
            no_proj_study += 1
        study_triples, summary_line = study.get_triples(namespace)
        triples.extend(study_triples)
        if summary_line:
            summaries.append(summary_line)
        study_count += 1
    print("There are " + str(study_count) + " studies.")
    if no_proj_study > 0:
        print(f"WARNING! There are {no_proj_study} studies without projects")
    return triples, summaries


def get_datasets(mwb_cur):
    print("Gathering Workbench Datasets")
    datasets = {}
    mwb_cur.execute("""\
        SELECT mb_sample_id, study_id, subject_species
        FROM metadata
        INNER JOIN subject
        ON metadata.subject_id = subject.subject_id""")
    for row in mwb_cur:
        dataset = Dataset()
        dataset.mb_sample_id = row[0]
        dataset.study_id = row[1]
        if row[2]:
            dataset.subject_species = row[2].replace('\n', '')
        datasets[dataset.mb_sample_id] = dataset
    return datasets


def make_datasets(namespace, datasets, studies):
    print("Making Workbench Datasets")
    dataset_triples = []
    study_triples = []
    dataset_count = 0
    no_study_datasets = 0
    for dataset in datasets.values():
        dataset.uri = namespace + dataset.mb_sample_id
        if dataset.study_id in studies.keys():
            parent_study = studies[dataset.study_id]
            study_uri = namespace + parent_study.study_id
            if dataset.subject_species not in parent_study.subject_species:
                parent_study.subject_species.append(dataset.subject_species)
        else:
            study_uri = None
            no_study_datasets += 1
        dataset_triples.extend(dataset.get_triples(study_uri))
        dataset_count += 1
    print("There will be " + str(dataset_count) + " new datasets.")
    if no_study_datasets > 0:
        print("There are {} datasets without studies"
              .format(no_study_datasets))
    for study in studies.values():
        study_triples.extend(study.get_species_triples(namespace))
    return dataset_triples, study_triples


def get_authors_pmid(sup_cur: db.Cursor, pmid: str) -> List[Dict[str, str]]:
    try:
        authors = []
        pubs = db.get_pubmed_publications(sup_cur, pmids=[pmid])
        if pmid not in pubs:
            print(f'PMID {pmid}: Could not find PubMed data')
            return []
        author_list = prefill.parse_author_list(pubs[pmid])
        for author in author_list:
            forename = author.findtext('ForeName', '').strip()
            surname = author.findtext('LastName', '').strip()
            auth = {
                'name': f'{forename} {surname}'.strip()
            }
            authors.append(auth)
        return authors
    except Exception:
        traceback.print_exc()
        print(f'Error parsing PubMed Data for tool with PMID {pmid}')
        return []


def get_yaml_tools(cfg: config.Config):
    try:
        tools_path = cfg.get('tools', 'tools.yaml')
        with open(tools_path, 'r') as tools_file:
            t = yaml.safe_load(tools_file)
            tools = []
            for tool_id, data in t.items():
                try:
                    tool = Tool(tool_id, data)
                    tools.append(tool)
                except Exception as e:
                    print(f'{e!r}')
                    print('Error: check configuration for tool "%s"' % tool_id)
                    continue
            return tools
    except Exception:
        print('Error parsing tools config file: %s' % tools_path)
        return []


def fetch_mtw_tools(sup_cur: db.Cursor) -> Iterable[Tool]:
    headers = tools
    for tool in tools.MetabolomicsToolsWiki.tools():
        name = tool[headers.SOFTWARE]
        pmid = tool[headers.PMID]
        authors = []
        if tool[headers.PMID].isnumeric():
            authors = get_authors_pmid(sup_cur, pmid)
        else:
            pmid = ""
        props = {
            "name": name,
            "description": tool[headers.DESCRIPTION],
            "license": {
                "kind": tool[headers.LICENSE],
            },
            "url": tool[headers.WEBSITE],
            "authors": authors,
            "pmid": pmid,
            "approach": tool[headers.APPROACHES],
            "functionality": tool[headers.FUNCTIONALITY],
            "instrumental": tool[headers.INSTRUMENT_DATA_TYPE],
            "language": tool[headers.LANGUAGE],
            "type": tool[headers.SOFTWARE_TYPE],
        }
        yield Tool(name, props)


def make_tools(
    namespace, tools: List[Tool], people, withheld_people, mwb_cur, sup_cur
):
    def find_person(name: str) -> int:
        try:
            first_name, last_name = name.split(' ', 1)
            matches = list(db.get_person(sup_cur, first_name, last_name))
            if len(matches) == 1:
                return matches[0]
        except Exception:
            pass
        return 0

    print("Generating triples for tools")
    triples = []
    tool_count = 0
    for tool in tools:
        # First, find all the authors' URIs
        non_matched_authors = tool.match_authors(find_person, namespace)
        if len(non_matched_authors) != 0:
            print(f"Not all authors matched for Tool: {tool.tool_id}")
            continue
        # Now, generate the triples.
        triples.extend(tool.get_triples(namespace))
        tool_count += 1
    print("There are " + str(tool_count) + " tools.")
    return triples


def print_to_file(triples: List[str], filename: str) -> None:
    with open(filename, 'a+') as file:
        print_to_open_file(triples, file)


def print_to_open_file(triples: List[str], file: IO) -> None:
    for spo in triples:
        # Replace LFs and CRs with escaped equivalent. Since N-Triples uses
        # " .\n" as a record-separator, these absolutely must be escaped.
        # This is mainly for PubMed titles and citations sanitization.
        # See https://www.w3.org/TR/n-triples/
        spo = re.sub(r'\n', r"\\n", spo)
        spo = re.sub(r'\r', r"\\r", spo)
        file.write(f"{spo} .\n")


def generate(config_path: str, old_path: str):
    timestamp = datetime.now()
    path = os.path.join("data_out",
                        timestamp.strftime("%Y"),
                        timestamp.strftime("%m"),
                        timestamp.strftime("%Y_%m_%d"))
    os.makedirs(path, exist_ok=True)

    org_file = os.path.join(path, 'orgs.nt')
    people_file = os.path.join(path, 'people.nt')
    project_file = os.path.join(path, 'projects.nt')
    study_file = os.path.join(path, 'studies.nt')
    dataset_file = os.path.join(path, 'datasets.nt')
    tools_file = os.path.join(path, 'tools.nt')
    photos_file = os.path.join(path, 'photos.nt')
    pubs_file = os.path.join(path, 'pubs.nt')
    add_file = os.path.join(path, 'add.nt')
    sub_file = os.path.join(path, 'sub.nt')

    cfg = config.load(config_path)

    if not cfg.namespace.endswith('/'):
        print(f"WARNING! Namespace doesn't end with '/': {cfg.namespace}")

    mwb_conn: psycopg2.extensions.connection = psycopg2.connect(
        host=cfg.get('mwb_host'), dbname=cfg.get('mwb_database'),
        user=cfg.get('mwb_username'), password=cfg.get('mwb_password'),
        port=cfg.get('mwb_port'))

    sup_conn: psycopg2.extensions.connection = psycopg2.connect(
        host=cfg.get('sup_host'), dbname=cfg.get('sup_database'),
        user=cfg.get('sup_username'), password=cfg.get('sup_password'),
        port=cfg.get('sup_port'))

    with mwb_conn, sup_conn:
        with mwb_conn.cursor() as mwb_cur, sup_conn.cursor() as sup_cur:
            # Organizations
            orgs = get_organizations(sup_cur)
            org_triples = make_organizations(cfg.namespace, orgs)
            print_to_file(org_triples, org_file)

            # People
            all_people = get_people(sup_cur)
            people = {k: v for k, v in all_people.items() if not v.withheld}
            withheld_people = {k: v
                               for k, v in all_people.items() if v.withheld}
            people_triples = make_people(cfg.namespace, people)
            people_triples.extend(
                link_people_to_org(cfg.namespace, sup_cur, people, orgs))
            print_to_file(people_triples, people_file)

            # Photos
            photos = get_photos(cfg.get("picturepath", "."), people)
            photos_triples = make_photos(cfg.namespace, photos)
            print_to_file(photos_triples, photos_file)

            # Tools
            print("Gathering Tools")
            yaml_tools = get_yaml_tools(cfg)
            csv_tools = list(fetch_mtw_tools(sup_cur))
            all_tools = yaml_tools + csv_tools
            tools_triples = make_tools(cfg.namespace, all_tools, people,
                                       withheld_people, mwb_cur, sup_cur)
            print_to_file(tools_triples, tools_file)

            # Publications
            pubs = get_publications(sup_cur, people)
            pubs_triples = make_publications(cfg.namespace, pubs)
            print_to_file(pubs_triples, pubs_file)

            # Projects
            projects = get_projects(mwb_cur, sup_cur, people, orgs)
            project_data = make_projects(cfg.namespace, projects)
            project_triples, project_summaries = project_data
            all_proj_triples = project_triples + project_summaries
            print_to_file(all_proj_triples, project_file)

            # Studies
            # Study file printed after datasets
            embargoed_path = cfg.get('embargoed', '')
            embargoed: List[str] = []
            if embargoed_path:
                with open(embargoed_path) as f:
                    embargoed = [line.strip() for line in f if line]

            studies = get_studies(mwb_cur, sup_cur, people, orgs, embargoed)
            study_triples, study_summaries = \
                make_studies(cfg.namespace, studies, projects)

            # Datasets
            datasets = get_datasets(mwb_cur)
            dataset_triples, study_sup_triples = \
                make_datasets(cfg.namespace, datasets, studies)
            print_to_file(dataset_triples, dataset_file)

            all_study_triples = study_triples + study_summaries \
                + study_sup_triples
            print_to_file(all_study_triples, study_file)

            if old_path:
                add, sub = diff(old_path, path)
                with open(add_file, 'w') as f:
                    f.writelines(add)
                with open(sub_file, 'w') as f:
                    f.writelines(sub)

    sup_conn.close()
    mwb_conn.close()


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      "hx:", ["help", "diff=", "add-devs"])
    except getopt.GetoptError:
        print(__doc__)
        sys.exit(2)

    old_path = ""

    for o, a in optlist:
        if o in ["-h", "--help"]:
            print(__doc__)
            sys.exit()
        elif o in ["-x", "--diff"]:
            old_path = a
            print("Differential update with previous run: " + old_path)
        elif o == "--add-devs":
            print("WARNING! --add-devs has been removed")

    if len(args) != 1:
        print(__doc__)
        sys.exit(2)

    config_path = args[0]
    generate(config_path, old_path)


if __name__ == "__main__":
    main()
