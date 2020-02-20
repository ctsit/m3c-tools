import datetime
import io
import typing

import psycopg2
import psycopg2.extensions


psql_cursor = psycopg2.extensions.cursor


INSTITUTE = 'institute'
DEPARTMENT = 'department'
LABORATORY = 'laboratory'

PUBMED_AUTHORSHIPS = "pubmed_authorships"


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
        SELECT COALESCE(institute, ''), COALESCE(department, ''),
               COALESCE(laboratory, ''), project_id
          FROM project
         UNION
        SELECT COALESCE(study.institute, ''), COALESCE(study.department, ''),
               COALESCE(study.laboratory, ''), study.study_id
          FROM study, study_status_prod
         WHERE study.study_id = study_status_prod.study_id
           AND study_status_prod.status = 1
    '''

    cursor.execute(select_orgs)

    for row in cursor:
        yield tuple(row)


def find_names(cursor: psycopg2.extensions.cursor,
               embargoed: typing.List[str]) \
        -> typing.Iterable[typing.Tuple[str, str, str, str, str, str]]:

    select_names = '''
        SELECT COALESCE(first_name, ''), COALESCE(last_name, ''),
               COALESCE(institute, ''), COALESCE(department, ''),
               COALESCE(laboratory, ''), project_id
          FROM project
         UNION
        SELECT COALESCE(study.first_name, ''), COALESCE(study.last_name, ''),
               COALESCE(study.institute, ''), COALESCE(study.department, ''),
               COALESCE(study.laboratory, ''), study.study_id
          FROM study, study_status_prod
         WHERE study.study_id = study_status_prod.study_id
           AND study_status_prod.status = 1
    '''

    cursor.execute(select_names)

    for row in cursor:
        yield tuple(row)


def get_affiliations(cursor: psql_cursor) -> typing.Mapping[int, str]:
    query = """
        SELECT p.id, o.name
          FROM organizations o,
               people p,
               associations a
         WHERE o.id = a.organization_id
           AND a.person_id = p.id
           AND o.withheld = FALSE
           AND o.type = 'institute'
    """

    cursor.execute(query)

    affiliations = {}
    for row in cursor:
        person_id = int(row[0])
        organization = row[1]
        if person_id not in affiliations:
            affiliations[person_id] = []
        affiliations[person_id].append(organization)

    return affiliations


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


def get_people(cursor: psql_cursor) \
        -> typing.Mapping[int, typing.Tuple[str, str]]:

    cursor.execute("""
        SELECT id, first_name, last_name
          FROM people p, names n
         WHERE p.id=n.person_id
    """)

    people: typing.Mapping[int, typing.Tuple[str, str]] = {}
    for row in cursor:
        person_id = int(row[0])
        first_name = str(row[1])
        last_name = str(row[2])
        people[person_id] = (first_name, last_name)

    return people


def get_pubmed_download_timestamps(cursor: psql_cursor) \
        -> typing.Mapping[str, datetime.datetime]:

    cursor.execute("""
        SELECT pmid, downloaded
          FROM pubmed_publications
    """)

    timestamps: typing.Mapping[int, datetime.datetime] = {}
    for row in cursor:
        pmid = row[0]
        downloaded = row[1]
        timestamps[pmid] = downloaded

    return timestamps


def update_authorships(cursor: psql_cursor,
                       authorships: typing.Dict[int, typing.List[str]]) -> int:
    tsv = io.StringIO()
    for person_id, pmids in authorships.items():
        for pmid in pmids:
            print(person_id, pmid, sep="\t", end="\n", file=tsv)

    tsv.seek(0)

    cursor.execute(f'TRUNCATE "{PUBMED_AUTHORSHIPS}"')
    cursor.copy_from(tsv, PUBMED_AUTHORSHIPS, columns=("person_id", "pmid"))

    return cursor.rowcount


def upsert_publication(cursor: psql_cursor, pmid: str, xml: str):
    ddl = """
        INSERT INTO pubmed_publications (pmid, xml, downloaded)
             VALUES                     (  %s,  %s, DEFAULT)
        ON CONFLICT (pmid)
        DO UPDATE SET xml=EXCLUDED.xml, downloaded=EXCLUDED.downloaded
          RETURNING pmid
    """

    cursor.execute(ddl, (pmid, xml))
    assert cursor.rowcount == 1
