"""
Metab Importer
Usage:
    python metab_prefill.py (-h | --help)
    python metab_prefill.py <path_to_config>

Options:
    -h --help       Show this message and exit

Instructions:
    See README

Example:
    $ python metab_prefill.py config.yaml
"""

import io
import operator
import sys
import typing
import xml.etree.ElementTree as ET

import psycopg2

import db
import metab_import
import tools

get_person = db.get_person  # Facilitate testing using monkeypatching.


INSTITUTE = 'institute'
DEPARTMENT = 'department'
LABORATORY = 'laboratory'


class AmbiguityError(Exception):
    def __init__(self, uid: str, institutes: str, departments: str, labs: str):
        msg = (f"Cannot determine organization hierarchy of {uid}: "
               f'institutes="{institutes}" '
               f'departments="{departments}" '
               f'labs="{labs}"')
        super().__init__(msg)


def add_developers(sup_conn: psycopg2.extensions.connection) -> None:
    with sup_conn.cursor() as cursor:
        pmids = set(tools.MetabolomicsToolsWiki.pmids())
        total = len(pmids)
        if total == 0:
            return

        publications = db.get_pubmed_publications(cursor, pmids)
        pmids = pmids.intersection(publications.keys())
        print(f'Found {len(pmids)} of {total} tools-related publications in '
              'the Supplemental database.')
        if len(pmids) == 0:
            return

        for pmid in pmids:
            author_list = parse_author_list(publications[pmid])
            for author in author_list:
                forename = author.findtext('ForeName', '').strip()
                lastname = author.findtext('LastName', '').strip()

                if not forename:
                    print(f'PMID {pmid}: missing forename of author '
                          f'{lastname}')
                    continue

                if not lastname:
                    print(f'PMID {pmid}: missing surname of author {forename}')
                    continue

                matches = list(db.get_person(cursor, forename, lastname))
                if len(matches) > 1:
                    print(f"PMID {pmid}: WARNING! Found {len(matches)} people "
                          f" named {forename} {lastname}: {' '.join(matches)}")
                    continue

                if matches:
                    pid = matches[0]
                    print(f'PMID {pmid}: found {forename} {lastname}: {pid}')
                else:
                    pid = db.add_person(cursor, forename, lastname, '', '')
                    if not pid:
                        print(f'PMID {pmid}: WARNING failed to add person: '
                              f'{forename} {lastname}')
                        continue
                    print(f'PMID {pmid}: added {forename} {lastname}: {pid}')

                affiliation_list = author.findall('.//Affiliation')
                for affiliation in affiliation_list:
                    affiliation = affiliation.text.strip()
                    if not affiliation:
                        continue
                    print(f"PMID {pmid}: affiliation for {forename} {lastname}"
                          f": {affiliation}")

    return


def add_organizations(mwb_conn: psycopg2.extensions.connection,
                      sup_conn: psycopg2.extensions.connection,
                      embargoed: typing.List[str]) -> None:

    with mwb_conn.cursor() as mwb_cur, sup_conn.cursor() as sup_cur:
        orgs = db.find_organizations(mwb_cur)

        for org in orgs:
            institutes, departments, laboratories, uid = [t.strip()
                                                          for t in org]

            # Exclude embargoed studies.
            if uid in embargoed:
                continue

            if not institutes:
                assert not departments
                assert not laboratories
                continue

            institute_list = [inst.strip() for inst in institutes.split(';')]
            for institute in institute_list:
                institute_id = db.get_organization(sup_cur, INSTITUTE,
                                                   institute)
                if not institute_id:
                    institute_id = db.add_organization(sup_cur, INSTITUTE,
                                                       institute)
                    print("Added institute #{}: {}."
                          .format(institute_id, institute))

            department_id: typing.Optional[int] = None
            if departments:
                department_list = [dept.strip()
                                   for dept in departments.split(';')]
                for department in department_list:
                    department_id = db.get_organization(sup_cur, DEPARTMENT,
                                                        department,
                                                        institute_id)

                    if department and not department_id:
                        department_id = db.add_organization(sup_cur,
                                                            DEPARTMENT,
                                                            department,
                                                            institute_id)
                        print("Added department #{}: {}."
                              .format(department_id, department))

            if not laboratories:
                continue

            laboratory_list = [lab.strip()
                               for lab in laboratories.split(';')]
            for laboratory in laboratory_list:
                laboratory_id = db.get_organization(sup_cur, LABORATORY,
                                                    laboratory, department_id)

                if not laboratory_id:
                    parent_id = department_id
                    if not parent_id:
                        # If there's no Department, the Institute is the parent
                        parent_id = institute_id
                    laboratory_id = db.add_organization(sup_cur, LABORATORY,
                                                        laboratory, parent_id)
                    print("Added laboratory #{}: {}."
                          .format(laboratory_id, laboratory))


def add_people(mwb_conn: psycopg2.extensions.connection,
               sup_conn: psycopg2.extensions.connection,
               embargoed: typing.List[str]) -> None:

    # For each (first_name, last_name) across both project and study:
    #   if (first_name, last_name) in names => skip
    #   add person to people table and name to names table.

    with mwb_conn.cursor() as mwb_cur, sup_conn.cursor() as sup_cur:
        people = db.find_people(mwb_cur)

        for row in people:
            (first_names, last_names,
             institutes, departments, labs,
             uid,
             emails, phones
             ) = tuple(row)

            # Exclude embargoed studies.
            if uid in embargoed:
                continue

            last_name_list = [ln.strip() for ln in last_names.split(';')]
            first_name_list = [fn.strip() for fn in first_names.split(';')]
            emails = [e.strip() for e in emails.split(';')]
            phones = [p.strip() for p in phones.split(';')]

            for i in range(0, len(last_name_list)):
                last_name = last_name_list[i]
                first_name = first_name_list[i]
                try:
                    email = emails[i]
                except IndexError:
                    email = ""  # Allow fewer emails than names.
                try:
                    phone = phones[i]
                except IndexError:
                    phone = ""  # Allow fewer phone numbers than names

                person_ids = get_person(sup_cur, first_name, last_name,
                                        exclude_withheld=False)
                if len(person_ids) > 1:
                    print("ERROR: multiple people with same name",
                          file=sys.stderr)
                    print("first={}. last={}. ids={}."
                          .format(first_name, last_name,
                                  ', '.join(person_ids)),
                          file=sys.stderr)
                    sys.exit(1)

                if len(person_ids) == 1:
                    person_id = person_ids[0]
                    details = db.get_contact_details(sup_cur, person_id)
                    curr_email, curr_phone = details

                    if email == curr_email and phone == curr_phone:
                        print(f"Skipping person #{person_id}: "
                              f"{first_name} {last_name}.")
                    else:
                        updated = db.update_contact_details(sup_cur, person_id,
                                                            email, phone)
                        if updated:
                            print("Updated contact details for person "
                                  f"#{person_id}: {first_name} {last_name}.")
                        else:
                            print("Failed to update contact details for person"
                                  f" #{person_id}: {first_name} {last_name}.")
                else:
                    person_id = db.add_person(sup_cur, first_name, last_name,
                                              email, phone)
                    print(f"Added person #{person_id}: {first_name} "
                          f"{last_name} (email={email}; phone={phone}).")

                # Create associations
                institute_list = [inst.strip()
                                  for inst in institutes.split(';')]
                try:
                    department_list = [dept.strip()
                                       for dept in departments.split(';')]
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
                        institute_id = db.get_organization(
                            sup_cur, INSTITUTE, institute_list[i])
                        assert institute_id
                        db.associate(sup_cur, person_id, institute_id)
                        parent_id = institute_id
                    except IndexError:
                        institute_id = db.get_organization(
                            sup_cur, INSTITUTE, institute_list[0])
                        assert institute_id
                        parent_id = institute_id

                    # If there are not enough departments, default to first
                    if departments:
                        try:
                            department_id = db.get_organization(
                                sup_cur, DEPARTMENT, department_list[i],
                                parent_id)

                            if department_id:
                                db.associate(sup_cur, person_id, department_id)
                                parent_id = department_id
                        except IndexError:
                            department_id = db.get_organization(
                                sup_cur, DEPARTMENT, department_list[0],
                                parent_id)
                            if department_id:
                                parent_id = department_id

                    if labs:
                        try:
                            laboratory_id = db.get_organization(
                                sup_cur, LABORATORY, lab_list[i], parent_id)
                            if laboratory_id:
                                db.associate(sup_cur, person_id, laboratory_id)
                        except IndexError:
                            print('WARNING: There are more institutes/'
                                  'departments than labs\n'
                                  'Institute list: ' + institutes + '\n'
                                  'Department list: ' + departments + '\n'
                                  'Lab list: ' + labs)

                print(('Associated {} {} (person) with '
                       '{} (institute) {} (department) {} (lab)')
                      .format(first_name, last_name, institutes, departments,
                              labs))

    return


def main():
    '''Adds people and organizations to the mwb_supplemental database.'''
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)

    if sys.argv[1] in ["-h", "--help"]:
        print(__doc__)
        sys.exit()

    config_path = sys.argv[1]
    config = metab_import.get_config(config_path)

    mwb_conn: psycopg2.extensions.connection = psycopg2.connect(
        host=config.get('mwb_host'), dbname=config.get('mwb_database'),
        user=config.get('mwb_username'), password=config.get('mwb_password'),
        port=config.get('mwb_port'))

    sup_conn: psycopg2.extensions.connection = psycopg2.connect(
        host=config.get('sup_host'), dbname=config.get('sup_database'),
        user=config.get('sup_username'), password=config.get('sup_password'),
        port=config.get('sup_port'))

    embargoed_path = config.get('embargoed', '')
    embargoed: typing.List[str] = []
    if embargoed_path:
        with open(embargoed_path) as f:
            embargoed = [line.strip() for line in f if line]

    with mwb_conn, sup_conn:
        add_organizations(mwb_conn, sup_conn, embargoed)
        add_people(mwb_conn, sup_conn, embargoed)
        add_developers(sup_conn)

    sup_conn.close()
    mwb_conn.close()


def parse_author_list(xml: str) -> ET.ElementTree:
    file = io.StringIO(xml)
    data = ET.parse(file)
    author_list = data.findall('//Article/AuthorList/Author')
    return author_list


if __name__ == "__main__":
    main()
