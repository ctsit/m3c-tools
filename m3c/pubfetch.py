"""M3C PubMed Fetcher

Usage:
    m3c pubfetch -h | --help
    m3c pubfetch [--authorships] [--delay=DELAY] [--max=MAX] <config>

Options:
    -h --help    Show this help message and exit.

Queries PubMed for a list of PMIDs for authors in the mwb_supplemental
database based on the authors' names and affiliations.

Afterwards, publications are batched and their XML summary is downloaded and
stored in the database along with their authorships.

If `--authorships` is specified, publications will not be downloaded.

`DELAY` is the number of seconds to wait between PubMed requests.

`MAX` is the maximum number of authorship searches to perform.

This is intended to be used by `metab_import.py` to generate Publication and
Authorship triples for VIVO.

Copyright 2020 University of Florida
"""

import datetime
import getopt
import http
import sys
import time
import traceback
import typing
import urllib.error
import xml.etree.ElementTree as ET

from Bio import Entrez
import psycopg2
import psycopg2.extensions

from m3c import catalyst
from m3c import config
from m3c import classes
from m3c import db
from m3c import tools


pubmed_delay: int = 0


def fetch_publications(cursor: db.Cursor):
    authorships = db.get_pubmed_authorships(cursor)
    tools_pmids = tools.MetabolomicsToolsWiki.pmids()
    pmids_to_download = set(tools_pmids).union(authorships.keys())
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
            batch = pmids[i:(i + BATCH_SIZE)]
            log(f"Downloading {i} through "
                f"{min(len(pmids), i+BATCH_SIZE)-1}")
            articles = pubmed_efetch(batch)
            for article in articles.getroot():
                try:
                    if article.tag == "PubmedBookArticle":
                        pmid = article.findtext("./BookDocument/PMID")
                    else:
                        pmid = article.findtext("./MedlineCitation/PMID")
                    assert pmid
                    xml = ET.tostring(article).decode("utf-8")
                    db.upsert_publication(cursor, pmid, xml)
                except Exception:
                    log(ET.tostring(article))
                    traceback.print_exc()
                    continue
            log(f"Batch done.")
        except Exception:
            traceback.print_exc()
            log(f"Error while processing PMIDs: {batch}")
            continue

    return


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
    help, config_path, only_update_authorships, delay, max_authorships = \
        parse_args(sys.argv)

    if help:
        print(__doc__)
        sys.exit()

    if not config_path:
        log(__doc__)
        sys.exit(2)

    pubfetch(config_path, only_update_authorships, delay, max_authorships)


def parse_args(argv) -> typing.Tuple[bool, str, bool, int, int]:
    try:
        opts, args = getopt.getopt(argv[1:],
                                   "h",
                                   ["help", "authorships", "delay=", "max="])
    except getopt.GetoptError:
        log(__doc__)
        sys.exit(2)

    help = False
    authorships = False
    delay = 0
    max_authorships = -1

    for opt, arg in opts:
        if opt in ["-h", "--help"]:
            help = True
            break
        if opt in ["-a", "--authorships"]:
            authorships = True
            continue
        try:
            if opt == "--delay":
                delay = abs(int(arg))
                continue
            if opt == "--max":
                max_authorships = abs(int(arg))
                continue
        except ValueError:
            log(f"error: invalid {opt}: {arg}")
            log(__doc__)
            sys.exit(2)

    if len(args) != 1:
        log("error: missing <config>")
        log(__doc__)
        sys.exit(2)

    config = args[0]

    return (help, config, authorships, delay, max_authorships)


def pubfetch(
    config_path: str,
    only_update_authorships: bool,
    delay: int,
    max_authorships: int
) -> None:
    global pubmed_delay
    pubmed_delay = delay

    cfg = config.load(config_path)

    pubmed_init(email=cfg.get("pubmed_email"),
                api_key=cfg.get("pubmed_api_token"))

    sup_conn: db.Connection
    sup_conn = psycopg2.connect(host=cfg.get("sup_host"),
                                dbname=cfg.get("sup_database"),
                                user=cfg.get("sup_username"),
                                password=cfg.get("sup_password"),
                                port=cfg.get("sup_port"))

    with sup_conn:
        with sup_conn.cursor() as cursor:
            update_authorships(cursor, max_authorships)
            if not only_update_authorships:
                fetch_publications(cursor)

    sup_conn.close()


def pubmed_efetch(id_list: typing.List[str]) -> ET.ElementTree:
    if pubmed_delay:
        log(f"Delaying PubMed efetch for {pubmed_delay} second(s).")
        time.sleep(pubmed_delay)
    ids = ",".join(id_list)
    handle = Entrez.efetch(db="pubmed", retmode="xml", id=ids)
    result = ET.parse(handle)
    handle.close()
    return result


def pubmed_esearch(term: str, retstart: int = 0, count_up: int = 0) \
        -> typing.List[str]:
    RETMAX = 100000

    if pubmed_delay and retstart == 0:
        log(f"Delaying PubMed esearch for {pubmed_delay} second(s).")
        time.sleep(pubmed_delay)

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


def pubmed_init(email: typing.Optional[str], api_key: typing.Optional[str]):
    Entrez.email = email
    Entrez.api_key = api_key


def too_recent(event: datetime.datetime,
               cutoff=datetime.timedelta(days=15)) -> bool:
    """Checks to see if event occurred more recently than cutoff."""
    now = datetime.datetime.now(event.tzinfo)
    return now - cutoff < event


def update_authorships(cursor: db.Cursor, authorships_limit: int = -1):
    people = db.get_people(cursor)
    affiliations = db.get_affiliations(cursor)
    known = db.get_confirmed_publications(cursor)
    authorships_updated = db.get_pubmed_authorships_updates(cursor)

    authorships: typing.Dict[int, typing.List[str]] = {}

    for person_id, info in people.items():
        if authorships_limit == 0:
            log("Reached the limit of authorship searches for this run")
            break

        latest_update = authorships_updated.get(person_id, None)
        if latest_update and too_recent(latest_update):
            log(f"{person_id}: skipping author whose publications were "
                "downloaded recently")
            continue

        person = classes.Person(person_id=str(person_id),
                                first_name=info[0],
                                last_name=info[1],
                                display_name=info[2],
                                email=info[3],
                                phone=info[4],
                                withheld=info[5],
                                overview=info[6])
        if person.withheld:
            log(f"{person_id}: skipping withheld author")
            continue

        if person_id not in affiliations:
            log(f"{person_id}: skipping author without affiliations")
            continue

        log(f"{person_id}: "
            f"fetching PMIDs for {person.first_name} {person.last_name}.")

        if person_id in known:
            pmids = catalyst.fetch_ids(person,
                                       list(affiliations[person_id]),
                                       include_pmids=list(known[person_id][0]),
                                       exclude_pmids=list(known[person_id][1]))
        else:
            pmids = get_pubmed_ids(person.first_name, person.last_name,
                                   list(affiliations[person_id]))

        authorships[person_id] = pmids

        log(f"{person_id}: found {len(pmids)} publications.")

        authorships_limit -= 1

    if not authorships:
        return

    count = db.update_authorships(cursor, authorships)
    log(f"Number of authorships updated: {count}")


if __name__ == "__main__":
    main()
