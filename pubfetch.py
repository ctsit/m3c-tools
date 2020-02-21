"""M3C PubMed Fethcher

Usage:
    python3 pubfetch.py -h | --help
    python3 pubfetch.py <config>

Options:
    -h --help    Show this message and exit.

Queries PubMed for a list of PMIDs for authors in the mwb_supplemental
database based on the authors' names and affiliations.

Afterwards, publications are batched and their XML summary is downloaded and
stored in the database along with their authorships.

This is intended to be used by `metab_import.py` to generate Publication and
Authorship triples for VIVO.

Copyright 2020 University of Florida
"""

import datetime
import http
import typing
import sys
import time
import traceback
import urllib.error
import xml.etree.ElementTree as ET

from Bio import Entrez
import psycopg2
import psycopg2.extensions

import catalyst
import db
import metab_import

psql_connection = psycopg2.extensions.connection
psql_cursor = psycopg2.extensions.cursor


def get_pubmed_ids(first_name: str, last_name: str,
                   affiliations: typing.List[str]) \
        -> typing.List[str]:
    """
    Get the PMIDs associated with a person with the passed affiliations.

    Returns an empty array if no affiliations are passed.

    Here is an example of a full query with a person with first name Arthur,
    last name Edison, and two affiliations.

    Edison Arthur[Author - Full] AND
    (University of Florida[Affiliation] OR University of Georgia[Affiliation])

    Which PubMed turns into:

    Edison, Arthur[Full Author Name] AND
    (University of Florida[Affiliation] OR University of Georgia[Affiliation])
    """
    if not affiliations:
        return []

    orgs = [f"{org}[Affiliation]" for org in affiliations]
    affiliation = " OR ".join(orgs)
    query = f"{first_name} {last_name}[Author - Full] AND ({affiliation})"

    try:
        pmids = pubmed_esearch(query)
        return pmids
    except urllib.error.HTTPError as err:
        if err.code != http.HTTPStatus.TOO_MANY_REQUESTS:
            raise err
        log("Too many requests to PubMed. Retrying after 1 second.")
        time.sleep(1)
        return get_pubmed_ids(first_name, last_name, affiliations)


def log(*values):
    """Prints values to stderr."""
    print(*values, file=sys.stderr)


def main():
    """Adds publications and authorships to the mwb_supplemental database."""
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)

    if sys.argv[1] in ["-h", "--help"]:
        print(__doc__)
        sys.exit()

    config_path = sys.argv[1]
    config = metab_import.get_config(config_path)

    Entrez.email = config.get("pubmed_email")
    Entrez.api_key = config.get("pubmed_api_token")

    sup_conn: psql_connection = psycopg2.connect(
        host=config.get("sup_host"),
        dbname=config.get("sup_database"),
        user=config.get("sup_username"),
        password=config.get("sup_password"),
        port=config.get("sup_port"))

    with sup_conn:
        with sup_conn.cursor() as cursor:
            pubfetch(cursor)

    sup_conn.close()


def pubfetch(cursor: psql_cursor):
    pmids_to_download: typing.MutableSet[str] = set()

    people = db.get_people(cursor)
    affiliations = db.get_affiliations(cursor)
    known = db.get_confirmed_publications(cursor)

    authorships: typing.Dict[int, typing.List[str]] = {}

    for person_id, info in people.items():
        person = metab_import.Person(person_id=person_id,
                                     first_name=info[0],
                                     last_name=info[1],
                                     display_name=info[2],
                                     email=info[3],
                                     phone=info[4],
                                     withheld=info[5])
        if person.withheld:
            continue

        log(f"{person_id}: "
            f"fetching PMIDs for {person.first_name} {person.last_name}.")

        if person_id in known:
            pmids = catalyst.fetch_ids(person,
                                       affiliations[person_id],
                                       include_pmids=known[person_id][0],
                                       exclude_pmids=known[person_id][1])
        else:
            pmids = get_pubmed_ids(person.first_name, person.last_name,
                                   affiliations[person_id])

        authorships[person_id] = pmids

        log(f"{person_id}: found {len(pmids)} publications.")

        pmids_to_download.update(pmids)

    count = db.update_authorships(cursor, authorships)
    log(f"Number of authorships updated: {count}")

    downloadts = db.get_pubmed_download_timestamps(cursor)
    skip = [pmid
            for pmid, downloaded
            in downloadts.items()
            if too_recent(downloaded)]
    pmids_to_download -= set(skip)

    log(f"Skipping {len(skip)} publications downloaded recently.")
    if pmids_to_download:
        log(f"Downloading XML for {len(pmids_to_download)} publications.")

    pmids = list(pmids_to_download)
    BATCH_SIZE = 5000
    for i in range(0, len(pmids), BATCH_SIZE):
        try:
            batch = pmids[i:i+BATCH_SIZE]
            log(f"Downloading {i} through "
                f"{min(len(pmids), i+BATCH_SIZE)-1}")
            articles = pubmed_efetch(batch)
            for article in articles.getroot():
                try:
                    pmid = article.find("./MedlineCitation/PMID").text
                    xml = ET.tostring(article).decode("utf-8")
                    db.upsert_publication(cursor, pmid, xml)
                except Exception:
                    traceback.print_exc()
                    continue
            log(f"Batch done.")
        except Exception:
            traceback.print_exc()
            log(f"Error while processing PMIDs: {batch}")
            continue

    return


def pubmed_efetch(id_list: typing.List[str]) -> ET.ElementTree:
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", retmode="xml", id=ids)
    result = ET.parse(handle)
    handle.close()
    return result


def pubmed_esearch(term: str, retstart: int = 0, count_up: int = 0) \
        -> typing.List[str]:
    RETMAX = 100000

    handle = Entrez.esearch(term=term,
                            db="pubmed",
                            retmax=RETMAX,
                            retstart=retstart)
    result = Entrez.read(handle)
    handle.close()

    id_list = result["IdList"]
    total = int(result["Count"])

    count_up += RETMAX
    # if # of results exceeds 100,000, you will need to run the query again
    if count_up < total:
        retstart += RETMAX
        id_list += pubmed_esearch(term, retstart, count_up)

    return id_list


def too_recent(event: datetime.datetime,
               cutoff=datetime.timedelta(days=15)) -> bool:
    """Checks to see if event occurred more recently than cutoff."""
    now = datetime.datetime.now(event.tzinfo)
    return now - cutoff < event


if __name__ == "__main__":
    main()
