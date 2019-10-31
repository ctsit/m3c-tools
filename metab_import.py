"""
Metab Importer
Usage:
    import.py (-h | --help)
    import.py [-d | --dry-run] [-x <prev> | --diff=<prev>] <path_to_config>

Options:
    -h --help     Show this message and exit
    -d --dry-run  Create N-Triples files without deleting and uploading to VIVO
    -x --diff     See Differential Update.

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


def get_config(config_path):
    try:
        with open(config_path, 'r') as config_file:
            config = yaml.load(config_file.read(), Loader=yaml.FullLoader)
    except Exception as e:
        print("Error: Check config file")
        sys.exit(e)
    return config


def connect(host, db, user, pg_password, port):
    conn = psycopg2.connect(host=host, dbname=db, user=user,
                            password=pg_password, port=port)
    cur = conn.cursor()
    return cur


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


def diff_upload(aide: Aide, add: typing.List[str], sub: typing.List[str]):
    lines = []
    for line in sub:
        line = line.strip()
        if line.endswith(" ."):
            line = line[:-2]
        lines.append(line)

    if lines:
        print(f"Differential update: removing {len(lines)} old triples")
        delete(aide, lines)

    lines = []
    for line in add:
        line = line.strip()
        if line.endswith(" ."):
            line = line[:-2]
        lines.append(line)

    if lines:
        print(f"Differential update: adding {len(lines)} new triples")
        insert(aide, lines)


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
        person = Person(person_id=row[0], first_name=row[1], last_name=row[2],
                        display_name=row[3], email=row[4], phone=row[5])
        people[person.person_id] = person
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
    for person in people.values():
        sup_cur.execute("""\
                    SELECT person_id, organization_id
                    FROM associations
                    WHERE person_id=%s""", (person.person_id,))
        for row in sup_cur:
            triples.extend(orgs[row[1]].add_person(namespace, person.uri))
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

        last_name = row[6]
        first_name = row[7]
        institute = row[8]
        department = row[9]
        lab = row[10]

        if institute:
            sup_cur.execute("""\
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND withheld = FALSE""",
                            (institute,))
            try:
                inst_id = sup_cur.fetchone()[0]
                project.institute = orgs[inst_id].org_id
            except TypeError:
                print("Error: Organization does not exist.")
                print("Organization for project " + project.project_id)
                print("Organization name: " + institute)
                sys.exit()
        if department:
            sup_cur.execute("""\
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND withheld = FALSE""",
                            (department,))
            try:
                dept_options = {}
                for row in sup_cur:
                    dept_options[row[0]] = row[1]
                for dept_id, parent in dept_options.items():
                    if inst_id == parent:
                        project.department = orgs[dept_id].org_id
            except TypeError:
                print("Error: Organization does not exist.")
                print("Organization for project " + project.project_id)
                print("Organization name: " + department)
                sys.exit()
        if lab:
            sup_cur.execute("""\
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND withheld = FALSE""",
                            (lab,))
            try:
                lab_options = {}
                for row in sup_cur:
                    lab_options[row[0]] = row[1]
                for lab_id, parent in lab_options.items():
                    if dept_id == parent:
                        project.lab = orgs[lab_id].org_id
            except TypeError:
                print("Error: Organization does not exist.")
                print("Organization for project " + project.project_id)
                print("Organization name: " + lab)
                sys.exit()

        sup_cur.execute("""\
                    SELECT person_id
                    FROM names
                    WHERE last_name=%s AND first_name=%s""",
                        (last_name, first_name))
        try:
            person_id = sup_cur.fetchone()[0]
            project.pi = people[person_id].person_id
        except KeyError:
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


def get_studies(mwb_cur, sup_cur, people, orgs):
    print("Gathering Workbench Studies")
    studies = {}
    mwb_cur.execute("""\
        SELECT study_id, study_title, COALESCE(study_type, ''),
            COALESCE(study_summary, ''), submit_date,
            project_id, last_name, first_name, institute, department,
            laboratory
        FROM study""")
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

        last_name = row[6]
        first_name = row[7]
        institute = row[8]
        department = row[9]
        lab = row[10]

        if institute:
            sup_cur.execute("""\
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND withheld = FALSE""",
                            (institute,))
            try:
                inst_id = sup_cur.fetchone()[0]
                study.institute = orgs[inst_id].org_id
            except TypeError:
                print("Error: Organization does not exist.")
                print("Organization for study " + study.study_id)
                print("Organization name: " + institute)
                sys.exit()
        if department:
            sup_cur.execute("""\
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND withheld = FALSE""",
                            (department,))
            try:
                dept_options = {}
                for row in sup_cur:
                    dept_options[row[0]] = row[1]
                for dept_id, parent in dept_options.items():
                    if inst_id == parent:
                        study.department = orgs[dept_id].org_id
            except TypeError:
                print("Error: Organization does not exist.")
                print("Organization for study " + study.study_id)
                print("Organization name: " + department)
                sys.exit()
        if lab:
            sup_cur.execute("""\
                        SELECT id, parent_id
                        FROM organizations
                        WHERE name=%s AND withheld = FALSE""",
                            (lab,))
            try:
                lab_options = {}
                for row in sup_cur:
                    lab_options[row[0]] = row[1]
                for lab_id, parent in lab_options.items():
                    if dept_id == parent:
                        study.lab = orgs[lab_id].org_id
            except TypeError:
                print("Error: Organization does not exist.")
                print("Organization for study " + study.study_id)
                print("Organization name: " + lab)
                sys.exit()

        sup_cur.execute("""\
                    SELECT person_id
                    FROM names
                    WHERE last_name=%s AND first_name=%s""",
                        (last_name, first_name))
        try:
            person_id = sup_cur.fetchone()[0]
            study.runner = people[person_id].person_id
        except IndexError:
            print("Error: Person does not exist.")
            print("Runner for study " + study.study_id)
            print("Last name: " + last_name)
            print("First name: " + first_name)
            sys.exit()

        studies[study.study_id] = study
    return studies


def make_studies(namespace, studies, projects):
    print("Making Workbench Studies")
    triples = []
    summaries = []
    study_count = 0
    no_proj_study = 0
    for study in studies.values():
        if study.project_id in projects.keys():
            project_uri = namespace + projects[study.project_id].project_id
        else:
            project_uri = None
            no_proj_study += 1
        study_triples, summary_line = study.get_triples(project_uri)
        triples.extend(study_triples)
        if summary_line:
            summaries.append(summary_line)
        study_count += 1
    print("There are " + str(study_count) + " studies.")
    if no_proj_study > 0:
        print("There are " + str(no_proj_study) + " studies without projects")
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


def make_tools(namespace, tools, people):
    print("Making Tools")
    triples = []
    tool_count = 0
    for tool in tools:
        # First, find all the authors' URIs
        if not tool.match_authors(people):
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


def delete(aide, triples, chunk_size=20):
    do_upload(aide, triples, chunk_size, "DELETE")


def insert(aide, triples, chunk_size=20):
    do_upload(aide, triples, chunk_size, "INSERT")


def do_upload(aide, triples, chunk_size=20, upload_type="INSERT"):
    assert upload_type in ["INSERT", "DELETE"]
    chunks = [triples[x:x+chunk_size]
              for x in range(0, len(triples), chunk_size)]
    for chunk in chunks:
        query = upload_type + """
            DATA {{
                GRAPH <http://vitro.mannlib.cornell.edu/default/vitro-kb-2> {{
                    {}
                }}
            }}
        """.format(" . \n".join(chunk))
        aide.do_update(query)


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)

    try:
        optlist, args = getopt.getopt(sys.argv[1:],
                                      "dhx:", ["dry-run", "help", "diff="])
    except getopt.GetoptError:
        print(__doc__)
        sys.exit(2)

    dry_run = False
    old_path = ""

    for o, a in optlist:
        if o in ["-h", "--help"]:
            print(__doc__)
            sys.exit()
        elif o in ["-d", "--dry-run"]:
            dry_run = True
            print("This is a dry run.")
        elif o in ["-x", "--diff"]:
            old_path = a
            print("Differential update with previous run: " + old_path)

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
    mwb_cur = connect(config.get('mwb_host'), config.get('mwb_database'),
                      config.get('mwb_username'), config.get('mwb_password'),
                      config.get('mwb_port'))
    sup_cur = connect(config.get('sup_host'), config.get('sup_database'),
                      config.get('sup_username'), config.get('sup_password'),
                      config.get('sup_port'))

    if not aide.namespace.endswith('/'):
        print(f"WARNING! Namespace doesn't end with '/': {aide.namespace}")

    # Organizations
    orgs = get_organizations(sup_cur)
    org_triples = make_organizations(aide.namespace, orgs)
    print_to_file(org_triples, org_file)

    # People
    people = get_people(sup_cur)
    people_triples = make_people(aide.namespace, people)
    people_triples.extend(link_people_to_org(aide.namespace, sup_cur, people, orgs))
    print_to_file(people_triples, people_file)

    # Photos
    photos = get_photos(config.get("picturepath", "."), people)
    photos_triples = make_photos(aide.namespace, photos)
    print_to_file(photos_triples, photos_file)

    # Tools
    yaml_tools = get_yaml_tools(config)
    csv_tools = get_csv_tools(config, aide)
    tools_triples = make_tools(aide.namespace, yaml_tools + csv_tools, people)
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

    all_study_triples = study_triples + study_summaries + study_sup_triples
    print_to_file(all_study_triples, study_file)

    summary_triples = project_summaries + study_summaries

    if old_path:
        add, sub = diff(old_path, path)
        with open(add_file, 'w') as f:
            f.writelines(add)
        with open(sub_file, 'w') as f:
            f.writelines(sub)

        if not dry_run:
            diff_upload(aide, add, sub)
            return

    if dry_run:
        sys.exit()

    # If you've made it this far, it's time to delete
    aide.do_delete()
    insert(aide, org_triples)
    print("Organizations uploaded")
    insert(aide, people_triples)
    print("People uploaded")
    insert(aide, project_triples)
    print("Projects uploaded")
    insert(aide, study_triples)
    print("Studies uploaded")
    insert(aide, dataset_triples)
    print("Datasets uploaded")
    insert(aide, tools_triples)
    print("Tools uploaded")
    insert(aide, photos_triples)
    print("Photos uploaded")
    insert(aide, summary_triples, 1)
    print("Summaries uploaded")


if __name__ == "__main__":
    main()
