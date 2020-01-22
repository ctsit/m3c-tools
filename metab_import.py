"""
Metab Importer
Usage:
    python3 metab_import.py (-h | --help)
    python3 metab_import.py [-x <prev> | --diff=<prev>] <path_to_config>

Options:
    -h --help      Show this message and exit
    -x --diff      See Differential Update.
    -a --add-devs  Add a Person record for new tools developers.

Differential Update:
    A differential update compares the triples produced by a run with that of
    an older one (<prev>). Two distinct sets of triples are produced:
     - an "Add" set which contains the set of triples that are in the current
       run's triples, but not in the older sets.
     - a "Sub" set which contains the set of triples that are in the older
       set, but not in the current run's.

    The corresponding files are written to add.nt and sub.nt.

Instructions:
    Run the importer where you have access to the postgres metabolomics
    database.
"""

from datetime import datetime
from typing import List
import csv
import getopt
import os
import pathlib
import sys
import time
import traceback
import typing
import yaml

import psycopg2

from aide import Aide
from metab_classes import Dataset
from metab_classes import Organization
from metab_classes import Person
from metab_classes import Photo
from metab_classes import Project
from metab_classes import Study
from metab_classes import Tool
import metab_data


def get_config(config_path):
    try:
        with open(config_path, 'r') as config_file:
            config = yaml.load(config_file.read(), Loader=yaml.FullLoader)
    except Exception as e:
        print("Error: Check config file")
        sys.exit(e)
    return config


def diff(prev_path: str, path: str) -> \
        typing.Tuple[typing.List[str], typing.List[str]]:
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

    previous = set(previous)
    current = set(current)

    add = current - previous
    sub = previous - current

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


def get_people(sup_cur):
    print("Gathering People")
    people = {}
    sup_cur.execute("""\
            SELECT id, first_name, last_name, display_name, email, phone
            FROM people p
            JOIN names n
            ON id=person_id
            WHERE p.withheld = FALSE AND n.withheld = FALSE""")
    for row in sup_cur:
        person = Person(person_id=row[0], first_name=row[1].strip(), last_name=row[2].strip(),
                        display_name=row[3], email=row[4], phone=row[5])
        people[person.person_id] = person
    return people


def add_person(cursor: psycopg2.extensions.cursor,
               first_name: str, last_name: str) -> int:

    statement = '''
        INSERT INTO people (id,      display_name)
             VALUES        (DEFAULT, %s          )
          RETURNING id
    '''

    display_name = '{} {}'.format(first_name.strip(), last_name.strip())
    cursor.execute(statement, (display_name,))

    person_id = cursor.fetchone()[0]

    statement = '''
        INSERT INTO names (person_id, first_name, last_name)
             VALUES       (%s       , %s        , %s       )
    '''

    cursor.execute(statement, (person_id, first_name.strip(), last_name.strip()))

    return person_id


def get_person(sup_cur, person_id):
    person = {}
    sup_cur.execute("""\
            SELECT id, first_name, last_name, display_name, email, phone
            FROM people p
            JOIN names n
            ON id=person_id
            WHERE p.withheld = FALSE AND n.withheld = FALSE AND id = %s""", (person_id,))
    for row in sup_cur:
        person = Person(person_id=row[0], first_name=row[1], last_name=row[2],
                        display_name=row[3], email=row[4], phone=row[5])
    return person


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


def get_projects(mwb_cur, sup_cur,
                 people: List[Person], orgs: List[Organization]):
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
            project_id=row[0].replace('\n', ''),
            project_title=row[1].replace('\n', '').replace('"', '\\"'),
            project_type=row[2].replace('\n', ''),
            summary=row[3].replace('\n', '').replace('"', '\\"'),
            doi=row[4].replace('\n', ''),
            funding_source=row[5].replace('\n', ''))

        last_names: str = row[6]
        first_names: str = row[7]
        institutes = row[8]
        departments = row[9]
        labs = row[10]

        institute_list = [inst for inst in institutes.split(';')]
        try:
            department_list = [dept for dept in departments.split(';')]
        except AttributeError:
            department_list = []
        try:
            lab_list = [lab for lab in labs.split(';')]
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
                sup_cur.execute("""\
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND type='institute'
                        AND withheld = FALSE""",
                            (institute_list[i],))
                try:
                    inst_id = sup_cur.fetchone()[0]
                    project.institutes.append(orgs[inst_id].org_id)
                except TypeError:
                    print("Error: Organization does not exist.")
                    print("Organization for project " + project.project_id)
                    print("Organization name: " + institute_list[i])
                    sys.exit()
            except IndexError:
                sup_cur.execute("""\
                            SELECT id, parent_id
                            FROM organizations
                            WHERE name=%s AND type='institute'
                            AND withheld = FALSE""",
                                (institute_list[0],))
                inst_id = sup_cur.fetchone()[0]

            # If there are not enough departments, default to first
            if departments:
                try:
                    sup_cur.execute("""\
                            SELECT id, parent_id
                            FROM organizations
                            WHERE name=%s AND type='department'
                            AND withheld = FALSE""",
                                (department_list[i],))
                    try:
                        dept_options = {}
                        for row in sup_cur:
                            dept_options[row[0]] = row[1]
                        for dept_id, parent in dept_options.items():
                            if inst_id == parent:
                                department_id = dept_id
                                project.departments.append(orgs[dept_id].org_id)
                    except TypeError:
                        print("Error: Organization does not exist.")
                        print("Organization for project " + project.project_id)
                        print("Organization name: " + department_list[i])
                        sys.exit()
                except IndexError:
                    sup_cur.execute("""\
                                SELECT id, parent_id
                                FROM organizations
                                WHERE name=%s AND type='department'
                                AND withheld = FALSE""",
                                    (department_list[0],))
                    dept_options = {}
                    for row in sup_cur:
                        dept_options[row[0]] = row[1]
                    for dept_id, parent in dept_options.items():
                        if inst_id == parent:
                            department_id = dept_id
            if labs:
                try:
                    sup_cur.execute("""\
                                SELECT id, parent_id
                                FROM organizations
                                WHERE name=%s AND type='laboratory'
                                AND withheld = FALSE""",
                                    (lab_list[i],))
                    try:
                        lab_options = {}
                        for row in sup_cur:
                            lab_options[row[0]] = row[1]
                        for lab_id, parent in lab_options.items():
                            try:
                                if department_id == parent:
                                    project.labs.append(orgs[lab_id].org_id)
                            except Exception:
                                import pdb
                                pdb.set_trace()
                    except TypeError:
                        print("Error: Organization does not exist.")
                        print("Organization for project " + project.project_id)
                        print("Organization name: " + lab_list[i])
                        sys.exit()
                except IndexError:
                    pass

        last_name_list = [ln for ln in last_names.split(';')]
        first_name_list = [fn for fn in first_names.split(';')]

        for i in range(0, len(last_name_list)):
            last_name = last_name_list[i]
            first_name = first_name_list[i]
            ids = metab_data.get_person(sup_cur, first_name, last_name)
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


def make_projects(namespace, projects: typing.Mapping[str, Project]):
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


def get_studies(mwb_cur, sup_cur, people, orgs):
    print("Gathering Workbench Studies")
    studies = {}
    mwb_cur.execute("""\
        SELECT study.study_id, study.study_title,
               COALESCE(study.study_type, ''),
               COALESCE(study.study_summary, ''), study.submit_date,
               study.project_id, study.last_name, study.first_name,
               study.institute, study.department, study.laboratory
        FROM study, study_status
        WHERE study.study_id = study_status.study_id
          AND study_status.status = 1""")

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

        # Skip invalid studies
        if not is_valid_study(study):
            continue

        last_names: str = row[6]
        first_names: str = row[7]
        institutes = row[8]
        departments = row[9]
        labs = row[10]

        institute_list = [inst for inst in institutes.split(';')]
        try:
            department_list = [dept for dept in departments.split(';')]
        except AttributeError:
            department_list = []
        try:
            lab_list = [lab for lab in labs.split(';')]
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
                sup_cur.execute("""\
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND type='institute'
                        AND withheld = FALSE""",
                            (institute_list[i],))
                try:
                    inst_id = sup_cur.fetchone()[0]
                    study.institutes.append(orgs[inst_id].org_id)
                except TypeError:
                    print("Error: Organization does not exist.")
                    print("Organization for study " + study.study_id)
                    print("Organization name: " + institute_list[i])
                    sys.exit()
            except IndexError:
                sup_cur.execute("""\
                            SELECT id, parent_id
                            FROM organizations
                            WHERE name=%s AND type='institute'
                            AND withheld = FALSE""",
                                (institute_list[0],))
                inst_id = sup_cur.fetchone()[0]

            # If there are not enough departments, default to first
            if departments:
                try:
                    sup_cur.execute("""\
                            SELECT id, parent_id
                            FROM organizations
                            WHERE name=%s AND type='department'
                            AND withheld = FALSE""",
                                (department_list[i],))
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
                    sup_cur.execute("""\
                                SELECT id, parent_id
                                FROM organizations
                                WHERE name=%s AND type='department'
                                AND withheld = FALSE""",
                                    (department_list[0],))
                    dept_options = {}
                    for row in sup_cur:
                        dept_options[row[0]] = row[1]
                    for dept_id, parent in dept_options.items():
                        if inst_id == parent:
                            department_id = dept_id
            if labs:
                try:
                    sup_cur.execute("""\
                                SELECT id, parent_id
                                FROM organizations
                                WHERE name=%s AND type='laboratory'
                                AND withheld = FALSE""",
                                    (lab_list[i],))
                    try:
                        lab_options = {}
                        for row in sup_cur:
                            lab_options[row[0]] = row[1]
                        for lab_id, parent in lab_options.items():
                            try:
                                if department_id == parent:
                                    study.labs.append(orgs[lab_id].org_id)
                            except Exception:
                                import pdb
                                pdb.set_trace()
                    except TypeError:
                        print("Error: Organization does not exist.")
                        print("Organization for study " + study.study_id)
                        print("Organization name: " + lab_list[i])
                        sys.exit()
                except IndexError:
                    pass

        last_name_list = [ln for ln in last_names.split(';')]
        first_name_list = [fn for fn in first_names.split(';')]

        for i in range(0, len(last_name_list)):
            last_name = last_name_list[i]
            first_name = first_name_list[i]

            ids = metab_data.get_person(sup_cur, first_name, last_name)
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


def make_studies(namespace, studies: typing.Dict[str, Study], projects):
    print("Making Workbench Studies")
    triples = []
    summaries = []
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


def get_authors_pmid(aide: Aide, pmid: typing.Text) -> typing.List[typing.Dict]:
    count = 0
    authors = []
    redo = True
    while redo and count < 2:
        count += 1
        redo = False
        try:
            data = aide.get_details([pmid])
            for author in data['PubmedArticle'][0]['MedlineCitation']['Article']['AuthorList']:
                authors.append({
                    'name': f"{author['ForeName'].split(' ')[0].strip()} {author['LastName']}"
                })
        except IOError:
            print('Error getting PubMed Data for tool with PMID %s. Trying again in 3 seconds' % pmid)
            redo = True
            time.sleep(3)
        except Exception:
            print('Error parsing PubMed Data for tool with PMID %s' % pmid)
    return authors


def get_yaml_tools(config):
    try:
        tools_path = config.get('tools', 'tools.yaml')
        with open(tools_path, 'r') as tools_file:
            t = yaml.load(tools_file.read(), Loader=yaml.FullLoader)
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


def strip_http(url: typing.Text) -> typing.Text:
    return url.replace('http://', '').replace('https://', '')


def get_csv_tools(config, aide: Aide) -> List[Tool]:
    try:
        csv_tools_path = config.get('tools_csv', 'tools.csv')
        with open(csv_tools_path, 'r') as tools_file:
            t = csv.reader(tools_file)
            # Skip the header row
            next(t)

            tools = []
            for tool_data in t:
                pmid = tool_data[19].strip()
                authors = None
                if pmid.isnumeric():
                    authors = get_authors_pmid(aide, pmid)
                if not tool_data[24].replace('-', '').strip():
                    continue
                tool = Tool(strip_http(tool_data[24]), {
                    'name': tool_data[21],
                    'description': tool_data[1],
                    'url': tool_data[24],
                    'authors': authors,
                    'pmid': pmid,
                    'tags': tool_data[6].split(',')
                })
                tools.append(tool)
            return tools
    except Exception as e:
        print(e)
        print('Error parsing tools config file: %s' % csv_tools_path)
        return []


def make_tools(namespace, tools: List[Tool], people, mwb_cur, sup_cur, add_devs):
    print("Making Tools")
    triples = []
    tool_count = 0
    for tool in tools:
        # First, find all the authors' URIs
        non_matched_authors = tool.match_authors(people, namespace)
        if len(non_matched_authors) != 0:
            if not add_devs:
                print('Not all authors matched for Tool. Use --add-devs to add these people. Skipping...')
                continue
            try:
                print("Not all authors matched for Tool. Inserting new people.")
                for author in non_matched_authors:
                    first_name = author.name.split(' ')[0]
                    last_name = " ".join(author.name.split(' ')[1:])
                    person_id = add_person(sup_cur, first_name, last_name)
                    print(f"Added {first_name} {last_name}: {person_id}")
                    people[person_id] = get_person(sup_cur, person_id)
                print("Trying to match authors again.")
                tool.match_authors(people, namespace)
            except Exception:
                traceback.print_exc()
                print('Error occurred creating new authors for tool. Skipping tool.')
                continue
        # Now, generate the triples.
        triples.extend(tool.get_triples(namespace))
        tool_count += 1
    print("There are " + str(tool_count) + " tools.")
    return triples


def print_to_file(triples, file):
    triples = [t + " ." for t in triples]
    with open(file, 'a+') as rdf:
        rdf.write("\n".join(triples))


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
    add_devs = False

    for o, a in optlist:
        if o in ["-h", "--help"]:
            print(__doc__)
            sys.exit()
        elif o in ["-x", "--diff"]:
            old_path = a
            print("Differential update with previous run: " + old_path)
        elif o == "--add-devs":
            add_devs = True
            print("Creating new developers for tools.")

    if len(args) != 1:
        print(__doc__)
        sys.exit(2)

    config_path = args[0]

    timestamp = datetime.now()
    path = 'data_out/' + timestamp.strftime("%Y") + '/' + \
        timestamp.strftime("%m") + '/' + timestamp.strftime("%Y_%m_%d")
    try:
        os.makedirs(path)
    except FileExistsError:
        pass
    org_file = os.path.join(path, 'orgs.nt')
    people_file = os.path.join(path, 'people.nt')
    project_file = os.path.join(path, 'projects.nt')
    study_file = os.path.join(path, 'studies.nt')
    dataset_file = os.path.join(path, 'datasets.nt')
    tools_file = os.path.join(path, 'tools.nt')
    photos_file = os.path.join(path, 'photos.nt')
    add_file = os.path.join(path, 'add.nt')
    sub_file = os.path.join(path, 'sub.nt')

    config = get_config(config_path)

    aide = Aide(config.get('update_endpoint'),
                config.get('vivo_email'),
                config.get('vivo_password'),
                config.get('namespace'))

    if not aide.namespace.endswith('/'):
        print(f"WARNING! Namespace doesn't end with '/': {aide.namespace}")

    mwb_conn: psycopg2.extensions.connection = psycopg2.connect(
        host=config.get('mwb_host'), dbname=config.get('mwb_database'),
        user=config.get('mwb_username'), password=config.get('mwb_password'),
        port=config.get('mwb_port'))

    sup_conn: psycopg2.extensions.connection = psycopg2.connect(
        host=config.get('sup_host'), dbname=config.get('sup_database'),
        user=config.get('sup_username'), password=config.get('sup_password'),
        port=config.get('sup_port'))

    with mwb_conn, sup_conn:
        with mwb_conn.cursor() as mwb_cur, sup_conn.cursor() as sup_cur:
            # Organizations
            orgs = get_organizations(sup_cur)
            org_triples = make_organizations(aide.namespace, orgs)
            print_to_file(org_triples, org_file)

            # People
            # Don't make the triples yet because tools can create new people.
            people = get_people(sup_cur)

            # Photos
            photos = get_photos(config.get("picturepath", "."), people)
            photos_triples = make_photos(aide.namespace, photos)
            print_to_file(photos_triples, photos_file)

            # Tools
            yaml_tools = get_yaml_tools(config)
            csv_tools = get_csv_tools(config, aide)
            tools_triples = make_tools(aide.namespace, yaml_tools + csv_tools, people, mwb_cur, sup_cur, add_devs)
            print_to_file(tools_triples, tools_file)

            # Projects
            projects = get_projects(mwb_cur, sup_cur, people, orgs)
            project_triples, project_summaries = \
                make_projects(aide.namespace, projects)
            all_proj_triples = project_triples + project_summaries
            print_to_file(all_proj_triples, project_file)

            # Studies
            # Study file printed after datasets
            studies = get_studies(mwb_cur, sup_cur, people, orgs)
            study_triples, study_summaries = \
                make_studies(aide.namespace, studies, projects)

            # Datasets
            datasets = get_datasets(mwb_cur)
            dataset_triples, study_sup_triples = \
                make_datasets(aide.namespace, datasets, studies)
            print_to_file(dataset_triples, dataset_file)

            all_study_triples = study_triples + study_summaries \
                + study_sup_triples
            print_to_file(all_study_triples, study_file)

            # Make People Triples
            people_triples = make_people(aide.namespace, people)
            people_triples.extend(
                link_people_to_org(aide.namespace, sup_cur, people, orgs))
            print_to_file(people_triples, people_file)

            if old_path:
                add, sub = diff(old_path, path)
                with open(add_file, 'w') as f:
                    f.writelines(add)
                with open(sub_file, 'w') as f:
                    f.writelines(sub)

    sup_conn.close()
    mwb_conn.close()


if __name__ == "__main__":
    main()
