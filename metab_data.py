"""
Library of common functions used by other metab modules.
"""

import typing

import psycopg2


def get_person(cursor: psycopg2.extensions.cursor,
               first_name: str, last_name: str,
               exclude_withheld: bool = True
               ) -> typing.List[int]:

    first_name = first_name.strip()
    last_name = last_name.strip()

    assert first_name and last_name

    query = '''
        SELECT person_id
          FROM names
         WHERE first_name=%s
           AND last_name=%s
    '''

    if exclude_withheld:
        query = f'{query} AND withheld=FALSE'

    cursor.execute(query, (first_name, last_name))

    ids = []
    for row in cursor:
        ids.append(row[0])

    return ids
