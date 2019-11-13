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

import sys
import typing

import psycopg2

import metab_import


INSTITUTE = 'institute'
DEPARTMENT = 'department'
LABORATORY = 'laboratory'


def associate(cursor: psycopg2.extensions.cursor,
              person_id: int, organization_id: int):
    insert_association = '''
        INSERT INTO associations (person_id, organization_id)
             VALUES              (%s       , %s             )
        ON CONFLICT DO NOTHING
    '''

    cursor.execute(insert_association, (person_id, organization_id))


def add_organization(cursor: psycopg2.extensions.cursor, type: str, name: str,
                     parent_id: typing.Optional[int] = None) -> int:
    assert type in [INSTITUTE, DEPARTMENT, LABORATORY]

    insert_org = '''
        INSERT INTO organizations (id     , name, type, parent_id)
             VALUES               (DEFAULT, %s  , %s  , %s       )
        ON CONFLICT (name, type, parent_id)
        DO UPDATE SET name=organizations.name
          RETURNING id
    '''

    cursor.execute(insert_org, (name, type, parent_id))
    assert cursor.rowcount == 1
    row = cursor.fetchone()

    return row[0]


def add_organizations(mwb_conn: psycopg2.extensions.connection,
                      sup_conn: psycopg2.extensions.connection):

    select_orgs = '''
        SELECT institute, department, laboratory
          FROM project
         UNION
        SELECT institute, department, laboratory
          FROM study
    '''

    with mwb_conn.cursor() as mwb_cur, sup_conn.cursor() as sup_cur:
        mwb_cur.execute(select_orgs)

        for row in mwb_cur:
            institutes, departments, laboratories = tuple(row)

            if not institutes:
                assert not departments
                assert not laboratories
                continue

            institute_list = [inst for inst in institutes.split(';')]
            for institute in institute_list:
                institute_id = get_organization(sup_cur, INSTITUTE, institute)
                if not institute_id:
                    institute_id = add_organization(sup_cur, INSTITUTE, institute)
                    print("Added institute #{}: {}."
                        .format(institute_id, institute))

            if departments:
                department_list = [dept for dept in departments.split(';')]
                for department in department_list:
                    department_id = get_organization(sup_cur, DEPARTMENT, department,
                                                    institute_id)

                    if department and not department_id:
                        department_id = add_organization(sup_cur, DEPARTMENT,
                                                        department, institute_id)
                        print("Added department #{}: {}."
                            .format(department_id, department))

            if not laboratories:
                continue

            laboratory_list = [lab for lab in laboratories.split(';')]
            for laboratory in laboratory_list:
                laboratory_id = get_organization(sup_cur, LABORATORY, laboratory,
                                                department_id)

                if not laboratory_id:
                    parent_id = department_id
                    if not parent_id:
                        # If there is no Department, the Institute is the parent.
                        parent_id = institute_id
                    laboratory_id = add_organization(sup_cur, LABORATORY,
                                                    laboratory, parent_id)
                    print("Added laboratory #{}: {}."
                        .format(laboratory_id, laboratory))


def add_people(mwb_conn: psycopg2.extensions.connection,
               sup_conn: psycopg2.extensions.connection):

    # For each (first_name, last_name) across both project and study:
    #   if (first_name, last_name) in names => skip
    #   add person to people table and name to names table.

    select_names = '''
        SELECT first_name, last_name, institute, department, laboratory
          FROM project
         UNION
        SELECT first_name, last_name, institute, department, laboratory
          FROM study
    '''

    with mwb_conn.cursor() as mwb_cur, sup_conn.cursor() as sup_cur:
        mwb_cur.execute(select_names)
        print("Found {} unique names".format(mwb_cur.rowcount), flush=True)

        for row in mwb_cur:
            first_names, last_names, institutes, departments, labs = tuple(row)
            last_name_list = [ln for ln in last_names.split(';')]
            first_name_list = [fn for fn in first_names.split(';')]

            for i in range(0, len(last_name_list)):
                last_name = last_name_list[i]
                first_name = first_name_list[i]

                person_ids = get_person(sup_cur, first_name, last_name)
                if len(person_ids) > 1:
                    print("ERROR: multiple people with same name", file=sys.stderr)
                    print("first={}. last={}. ids={}."
                        .format(first_name, last_name, ', '.join(person_ids)),
                        file=sys.stderr)
                    sys.exit(1)

                if len(person_ids) == 1:
                    person_id = person_ids[0]
                    print("Skipping person #{}: {} {}."
                        .format(person_ids[0], first_name, last_name))
                else:
                    person_id = add_person(sup_cur, first_name, last_name)
                    print("Added person #{}: {} {}."
                        .format(person_id, first_name, last_name))

                # Create associations
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
                        institute_id = get_organization(sup_cur, INSTITUTE, institute_list[i])
                        assert institute_id
                        associate(sup_cur, person_id, institute_id)
                        parent_id = institute_id
                    except IndexError:
                        institute_id = get_organization(sup_cur, INSTITUTE, institute_list[0])
                        assert institute_id
                        parent_id = institute_id
                    
                    if departments:
                        try:
                            department_id = get_organization(sup_cur, DEPARTMENT, department_list[i],
                                                    parent_id)
                            if department_id:
                                associate(sup_cur, person_id, department_id)
                                parent_id = department_id
                        except IndexError:
                            department_id = get_organization(sup_cur, DEPARTMENT, department_list[0],
                                                    parent_id)
                            if department_id:
                                parent_id = department_id

                    if labs:
                        try:
                            laboratory_id = get_organization(sup_cur, LABORATORY, lab_list[i],
                                                    parent_id)
                            if laboratory_id:
                                associate(sup_cur, person_id, laboratory_id)
                        except IndexError:
                            print('WARNING: There are more institutes/departments than labs\n'
                                'Institute list: ' + institutes + '\n'
                                'Department list: ' + departments + '\n'
                                'Lab list: ' + labs)

                print(('Associated {} {} (person) with '
                    '{} (institute) {} (department) {} (lab)')
                    .format(first_name, last_name, institutes, departments, labs))

    return


def add_person(cursor: psycopg2.extensions.cursor,
               first_name: str, last_name: str) -> int:

    statement = '''
        INSERT INTO people (id,      display_name)
             VALUES        (DEFAULT, %s          )
          RETURNING id
    '''

    display_name = '{} {}'.format(first_name, last_name)
    cursor.execute(statement, (display_name,))

    person_id = cursor.fetchone()[0]

    statement = '''
        INSERT INTO names (person_id, first_name, last_name)
             VALUES       (%s       , %s        , %s       )
    '''

    cursor.execute(statement, (person_id, first_name, last_name))

    return person_id


def get_organization(cursor: psycopg2.extensions.cursor, type: str, name: str,
                     parent_id: typing.Optional[int] = None) -> int:
    assert type in [INSTITUTE, DEPARTMENT, LABORATORY]

    if not parent_id:
        select_org = '''
            SELECT id FROM organizations
             WHERE name=%s AND type=%s AND parent_id IS NULL
        '''
        cursor.execute(select_org, (name, type))
    else:
        select_org = '''
            SELECT id FROM organizations
             WHERE name=%s AND type=%s AND parent_id=%s
        '''
        cursor.execute(select_org, (name, type, parent_id))

    assert cursor.rowcount <= 1

    row = cursor.fetchone()
    if not row:
        return 0

    return row[0]


def get_person(cursor: psycopg2.extensions.cursor,
               first_name: str, last_name: str) -> typing.List[int]:

    query = '''
        SELECT person_id
          FROM names
         WHERE first_name=%s
           AND last_name=%s
    '''

    cursor.execute(query, (first_name, last_name))

    ids = []
    for row in cursor:
        ids.append(row[0])

    return ids


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

    with mwb_conn, sup_conn:
        add_organizations(mwb_conn, sup_conn)
        add_people(mwb_conn, sup_conn)

    sup_conn.close()
    mwb_conn.close()


if __name__ == "__main__":
    main()
