import typing

import psycopg2


INSTITUTE = 'institute'
DEPARTMENT = 'department'
LABORATORY = 'laboratory'


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


def find_organizations(cursor: psycopg2.extensions.cursor,
                       embargoed: typing.List[str]) \
        -> typing.Iterable[typing.Tuple[str, str, str, str]]:

    select_orgs = '''
        SELECT institute, department, laboratory, project_id
          FROM project
         UNION
        SELECT study.institute, department, laboratory, study.study_id
          FROM study, study_status_prod
         WHERE study.study_id = study_status_prod.study_id
           AND study_status_prod.status = 1
    '''

    cursor.execute(select_orgs)

    for row in cursor:
        yield tuple(row)


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
