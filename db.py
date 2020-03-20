import datetime
import io
from typing import Iterable, Mapping, Optional, Tuple

import psycopg2
import psycopg2.extensions


Cursor = psycopg2.extensions.cursor
Connection = psycopg2.extensions.connection


INSTITUTE = 'institute'
DEPARTMENT = 'department'
LABORATORY = 'laboratory'


def add_organization(cursor: Cursor, type: str, name: str,
                     parent_id: Optional[int] = None) -> int:
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


def add_person(cursor: Cursor,
               first_name: str, last_name: str, email: str, phone: str) -> int:

    statement = '''
        INSERT INTO people (id,      display_name, email, phone)
             VALUES        (DEFAULT, %s          , %s   , %s)
          RETURNING id
    '''

    first_name = first_name.strip()
    last_name = last_name.strip()

    display_name = f'{first_name} {last_name}'
    cursor.execute(statement, (display_name, email, phone))

    person_id = cursor.fetchone()[0]

    statement = '''
        INSERT INTO names (person_id, first_name, last_name)
             VALUES       (%s       , %s        , %s       )
    '''

    cursor.execute(statement, (person_id, first_name, last_name))

    return person_id


def associate(cursor: Cursor, person_id: int, organization_id: int):
    insert_association = '''
        INSERT INTO associations (person_id, organization_id)
             VALUES              (%s       , %s             )
        ON CONFLICT DO NOTHING
    '''
    cursor.execute(insert_association, (person_id, organization_id))


def find_organizations(cursor: Cursor) \
        -> Iterable[Tuple[str, str, str, str]]:

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


def find_people(cursor: Cursor) \
        -> Iterable[Tuple[str, str, str, str, str, str, str, str]]:

    select_names = '''
        SELECT COALESCE(first_name, ''), COALESCE(last_name, ''),
               COALESCE(institute, ''), COALESCE(department, ''),
               COALESCE(laboratory, ''), project_id,
               COALESCE(email, ''), COALESCE(phone, '')
          FROM project
         UNION
        SELECT COALESCE(study.first_name, ''), COALESCE(study.last_name, ''),
               COALESCE(study.institute, ''), COALESCE(study.department, ''),
               COALESCE(study.laboratory, ''), study.study_id,
               COALESCE(email, ''), COALESCE(phone, '')
          FROM study, study_status_prod
         WHERE study.study_id = study_status_prod.study_id
           AND study_status_prod.status = 1
    '''

    cursor.execute(select_names)

    for row in cursor:
        yield tuple(row)


def get_affiliations(cursor: Cursor) -> Mapping[int, str]:
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


def get_confirmed_publications(cursor: Cursor) \
        -> Mapping[int, Tuple[Iterable[str], Iterable[str]]]:

    publications: Mapping[int, Tuple[Iterable[str], Iterable[str]]] = {}

    select = """
        SELECT person_id, pmid, include
        FROM publications
    """
    cursor.execute(select)

    for row in cursor:
        person_id = int(row[0])
        pmid = str(row[1])
        include = bool(row[2])

        if person_id not in publications:
            publications[person_id] = ([], [])

        if include:
            publications[person_id][0].append(pmid)
        else:
            publications[person_id][1].append(pmid)

    return publications


def get_contact_details(cursor: Cursor, person_id: int) -> Tuple[str, str]:
    query = '''
        SELECT COALESCE(email, ''), COALESCE(phone, '')
          FROM people
         WHERE id=%s
    '''
    cursor.execute(query, (person_id,))
    assert cursor.rowcount == 1
    for row in cursor:
        return row


def get_organization(cursor: Cursor, type: str, name: str,
                     parent_id: Optional[int] = None) -> int:
    assert type in [INSTITUTE, DEPARTMENT, LABORATORY]

    if not parent_id:
        select_org = """
            SELECT id FROM organizations
             WHERE name=%s AND type=%s AND parent_id IS NULL
        """
        cursor.execute(select_org, (name, type))
    else:
        select_org = """
            SELECT id FROM organizations
             WHERE name=%s AND type=%s AND parent_id=%s
        """
        cursor.execute(select_org, (name, type, parent_id))

    assert cursor.rowcount <= 1

    row = cursor.fetchone()
    if not row:
        return 0

    return row[0]


def get_organizations(cursor: Cursor) \
        -> Iterable[Tuple[int, str, str, int, bool]]:
    query = "SELECT id, name, type, parent_id, withheld FROM organizations"
    cursor.execute(query)
    for row in cursor:
        yield row


def get_person(cursor: Cursor, first_name: str, last_name: str,
               exclude_withheld: bool = True) \
        -> Iterable[int]:

    first_name = first_name.strip()
    last_name = last_name.strip()

    assert first_name and last_name

    query = '''
        SELECT person_id
          FROM names
         WHERE CONCAT(first_name, ' ', last_name) = %s
    '''

    if exclude_withheld:
        query = f'{query} AND withheld=FALSE'

    cursor.execute(query, (f'{first_name} {last_name}', ))

    ids = []
    for row in cursor:
        ids.append(row[0])

    return ids


def get_people(cursor: Cursor) \
        -> Mapping[int, Tuple[str, str, str, str, str, bool]]:
    select_names = """
        SELECT id, first_name, last_name, COALESCE(display_name, ''),
               COALESCE(email, ''), COALESCE(phone, ''), p.withheld
          FROM people p, names n
         WHERE p.id=n.person_id
    """

    cursor.execute(select_names)

    people: Mapping[int, Tuple[str, str, str, str, str, bool]] = {}
    for row in cursor:
        person_id = int(row[0])
        people[person_id] = tuple(row[1:7])

    return people


def get_pubmed_authorships(cursor: Cursor) -> Mapping[str, Iterable[int]]:
    select_pubs = """
        SELECT pmid, person_id
          FROM pubmed_authorships
    """

    cursor.execute(select_pubs)

    authorships = {}
    for row in cursor:
        pmid, person_id = row[0:2]
        if pmid not in authorships:
            authorships[pmid] = []
        authorships[pmid].append(person_id)

    return authorships


def get_pubmed_authorships_updates(cursor: Cursor) \
        -> Mapping[str, datetime.datetime]:
    select_pubs = """
        SELECT person_id, updated
          FROM pubmed_authorships_updates
    """

    cursor.execute(select_pubs)

    return {row[0]: row[1] for row in cursor}


def get_pubmed_download_timestamps(cursor: Cursor) \
        -> Mapping[str, datetime.datetime]:
    select_pubs = """
        SELECT pmid, downloaded
          FROM pubmed_publications
    """

    cursor.execute(select_pubs)

    return {row[0]: row[1] for row in cursor}


def get_pubmed_publications(cursor: Cursor,
                            pmids: Optional[Iterable[str]] = []) \
        -> Mapping[str, str]:
    if pmids:
        select_pubs = """
            SELECT pmid, xml
              FROM pubmed_publications
             WHERE pmid = ANY(%s)
        """
        cursor.execute(select_pubs, (list(pmids),))
    else:
        select_pubs = """
            SELECT pmid, xml
            FROM pubmed_publications
        """
        cursor.execute(select_pubs)
    return {row[0]: row[1] for row in cursor}


def update_authorships(cursor: Cursor,
                       authorships: Mapping[int, Iterable[str]]) -> int:
    """Overwrite authorships and timestamp the updates."""

    delete = """
        DELETE FROM pubmed_authorships_updates
              WHERE person_id = ANY(%s);
        DELETE FROM pubmed_authorships
              WHERE person_id = ANY(%s);
    """
    person_ids = list(authorships.keys())
    cursor.execute(delete, (person_ids, person_ids))

    # Update timestamps for those authors who have been updated.
    ids = io.StringIO()
    for person_id in authorships.keys():
        print(person_id, file=ids)
    ids.seek(0)
    cursor.copy_from(ids, "pubmed_authorships_updates", columns=("person_id",))

    # Update publication lists for those authors who have been updated.
    tsv = io.StringIO()
    for person_id, pmids in authorships.items():
        for pmid in pmids:
            print(person_id, pmid, sep="\t", end="\n", file=tsv)
    tsv.seek(0)
    cursor.copy_from(tsv, "pubmed_authorships", columns=("person_id", "pmid"))
    return cursor.rowcount


def update_contact_details(cursor: Cursor,
                           person_id: int, email: str, phone: str) -> bool:
    update = '''
        UPDATE people
           SET email=%s,
               phone=%s
         WHERE id=%s
    '''
    cursor.execute(update, (email, phone, person_id))
    return cursor.rowcount == 1


def upsert_publication(cursor: Cursor, pmid: str, xml: str) -> None:
    insert = """
        INSERT INTO pubmed_publications (pmid, xml, downloaded)
             VALUES                     (  %s,  %s, DEFAULT)
        ON CONFLICT (pmid)
        DO UPDATE SET xml=EXCLUDED.xml, downloaded=EXCLUDED.downloaded
          RETURNING pmid
    """

    cursor.execute(insert, (pmid, xml))
    assert cursor.rowcount == 1
