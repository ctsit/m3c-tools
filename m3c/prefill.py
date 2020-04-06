"""
Metab Importer
Usage:
    python m3c/prefill.py (-h | --help)
    python m3c/prefill.py <path_to_config>

Options:
    -h --help       Show this message and exit

Instructions:
    See README

Example:
    $ python m3c/prefill.py config.yaml
"""

import io
import itertools
import sys
import typing
import xml.etree.ElementTree as ET

import psycopg2

import db
from m3c import mwb
import metab_import
import tools


List = typing.List
Tuple = typing.Tuple


get_person = db.get_person  # Facilitate testing using monkeypatching.


class AmbiguityError(Exception):
    pass


class AmbiguousHierarchyError(AmbiguityError):
    def __init__(self, record: mwb.NameRecord):
        msg = (f"Cannot determine organization hierarchy of {record.psid}: "
               f'institutes="{record.institute}" '
               f'departments="{record.department}" '
               f'labs="{record.laboratory}"')
        super().__init__(msg)


class AmbiguousNamesError(AmbiguityError):
    def __init__(self, record: mwb.NameRecord):
        msg = (f"Number of name parts differs for {record.psid}: "
               f'first_name="{record.first_name}" '
               f'last_name="{record.last_name}"')
        super().__init__(msg)


def add_developers(sup_cur: db.Cursor) -> None:
    pmids = set(tools.MetabolomicsToolsWiki.pmids())
    total = len(pmids)
    if total == 0:
        return

    publications = db.get_pubmed_publications(sup_cur, pmids)
    pmids = pmids.intersection(publications.keys())
    print(f"Found {len(pmids)} of {total} tools-related publications in the "
          "Supplemental database.")
    if len(pmids) == 0:
        return

    for pmid in pmids:
        author_list = parse_author_list(publications[pmid])
        for author in author_list:
            forename = author.findtext("ForeName", "").strip()
            lastname = author.findtext("LastName", "").strip()

            if not forename:
                print(f"PMID {pmid}: missing forename of author {lastname}")
                continue

            if not lastname:
                print(f"PMID {pmid}: missing surname of author {forename}")
                continue

            matches = list(db.get_person(sup_cur, forename, lastname))
            if len(matches) > 1:
                print(f"PMID {pmid}: WARNING! Found {len(matches)} people "
                      f" named {forename} {lastname}: {matches}")
                continue

            if matches:
                pid = matches[0]
                print(f"PMID {pmid}: found {forename} {lastname}: {pid}")
            else:
                pid = db.add_person(sup_cur, forename, lastname, "", "")
                if not pid:
                    print(f"PMID {pmid}: WARNING failed to add person: "
                          f"{forename} {lastname}")
                    continue
                print(f"PMID {pmid}: added {forename} {lastname}: {pid}")

            affiliation_list = author.findall(".//Affiliation")
            for affiliation in affiliation_list:
                affiliation = affiliation.text.strip()
                if not affiliation:
                    continue
                print(f"PMID {pmid}: affiliation for {forename} {lastname}"
                      f": {affiliation}")

    return


def add_organizations(sup_cur: db.Cursor,
                      record: mwb.NameRecord
                      ) -> List[Tuple[int, int, int]]:
    """
    Adds institutes, departments, and labs to the Supplemental database.

    Returns the IDs as a 3-tuple (institute ID, department ID, lab ID).
    """
    if not record.institute:
        assert not record.department
        assert not record.laboratory
        return

    institutes = [inst.strip() for inst in record.institute.split(';')]
    departments = [dept.strip() for dept in record.department.split(';')]
    laboratories = [lab.strip() for lab in record.laboratory.split(';')]

    # Special case: allow department field is be completely empty.
    if record.department.strip() == "":
        departments = [""] * len(institutes)
    # Special case: allow laboratory field is be completely empty.
    if record.laboratory.strip() == "":
        laboratories = [""] * len(institutes)

    icnt = len(institutes)
    dcnt = len(departments)
    lcnt = len(laboratories)

    if not(dcnt in [1, lcnt] and icnt in [1, dcnt]):
        raise AmbiguousHierarchyError(record)

    institute_ids = []
    for i, institute in enumerate(institutes):
        oid = db.get_organization(sup_cur, mwb.INSTITUTE, institute)
        if not oid:
            oid = db.add_organization(sup_cur, mwb.INSTITUTE, institute)
            print(record.psid, f"added institute #{oid}: {institute}.")
        else:
            print(record.psid, f"found institute #{oid}: {institute}.")
        assert oid
        institute_ids.append(oid)

    department_ids = []
    for i, department in enumerate(departments):
        if not department:
            department_ids.append(0)
            continue

        if icnt == 1:
            parent = institute_ids[0]
        elif icnt > 1:
            parent = institute_ids[i]
        assert parent > 0

        oid = db.get_organization(sup_cur, mwb.DEPARTMENT, department, parent)

        if not oid:
            oid = db.add_organization(sup_cur, mwb.DEPARTMENT, department,
                                      parent)
            print(record.psid,
                  f"added department #{oid}: {department} (parent #{parent}).")
        else:
            print(record.psid, f"found department #{oid}: {department}.")

        assert oid
        department_ids.append(oid)

    laboratory_ids = []
    for i, laboratory in enumerate(laboratories):
        if not laboratory:
            laboratory_ids.append(0)
            continue

        if dcnt == 1:
            parent = department_ids[0]
        elif dcnt > 1:
            parent = department_ids[i]
        if parent == 0:
            # If there's no Department, the Institute is the parent.
            if icnt == 1:
                parent = institute_ids[0]
            elif icnt > 0:
                parent = institute_ids[i]
        assert parent > 0

        oid = db.get_organization(sup_cur, mwb.LABORATORY, laboratory, parent)

        if not oid:
            oid = db.add_organization(sup_cur, mwb.LABORATORY, laboratory,
                                      parent)
            print(record.psid,
                  f"added laboratory #{oid}: {laboratory} (parent #{parent}).")
        else:
            print(record.psid, f"found laboratory #{oid}: {laboratory}.")

        assert oid
        laboratory_ids.append(oid)

    ids: List[Tuple[int, int, int]] = []

    for i, laboratory_id in enumerate(laboratory_ids):
        institute_id = institute_ids[0]
        if icnt > 1:
            institute_id = institute_ids[i]

        department_id = department_ids[0]
        if dcnt > 1:
            department_id = department_ids[i]

        ids.append((institute_id, department_id, laboratory_id))

    return ids


def add_people(sup_cur: db.Cursor, record: mwb.NameRecord) -> List[int]:
    """ Adds people to the Supplemental database, returning their IDs. """
    ids: List[int] = []

    psid = record.psid
    last_names = [ln.strip() for ln in record.last_name.split(';')]
    first_names = [fn.strip() for fn in record.first_name.split(';')]
    emails = [e.strip() for e in record.email.split(';')]
    phones = [p.strip() for p in record.phone.split(';')]

    if len(last_names) != len(first_names):
        raise AmbiguousNamesError(record)

    if len(emails) > len(last_names):
        error(psid, "too many email addresses for", record.email)
        emails = []  # Don't try to add ANY email addresses for this record.
    if len(phones) > len(last_names):
        error(psid, "too many phone numbers for", record.phone)
        phones = []  # Don't try to add ANY phone numbers for this record.

    # It's acceptable to have fewer emails and phone numbers than names.
    combined = itertools.zip_longest(last_names, first_names, emails, phones,
                                     fillvalue="")

    for last_name, first_name, email, phone in combined:
        person_ids = list(get_person(sup_cur, first_name, last_name,
                                     exclude_withheld=False))
        if bad_email(email):
            error(psid, f"bad email address: {record.email}")
            email = ""

        if len(person_ids) > 1:
            pid = 0
            error(psid, "multiple people with the same name:",
                  f'first="{first_name}',
                  f'last="{last_name}"',
                  f'ids={person_ids}')
        elif len(person_ids) == 1:
            pid = person_ids[0]
            details = db.get_contact_details(sup_cur, pid)
            curr_email, curr_phone = details
            print(psid, f"found person #{pid}: {first_name} {last_name}.")
            if email != curr_email or phone != curr_phone:
                updated = db.update_contact_details(sup_cur, pid, email, phone)
                if updated:
                    print(psid, "updated contact details for person",
                          f"#{pid}: {first_name} {last_name}.")
                else:
                    error(psid, "failed to update contact details for person",
                          f"#{pid}: {first_name} {last_name}.")
        else:
            pid = db.add_person(sup_cur, first_name, last_name, email, phone)
            print(psid, f"Added person #{pid}: {first_name} {last_name}",
                  f"(email={email}; phone={phone}).")

        ids.append(pid)

    return ids


def associate(sup_cur: db.Cursor, psid: str, person_id: int, institute_id: int,
              department_id: int, laboratory_id: int):
    if person_id == 0:
        return

    templates = ["person #{} already associated with {} #{}.",
                 "associated person #{} with {} #{}."]

    if institute_id > 0:
        t = db.associate(sup_cur, person_id, institute_id)
        tmpl = templates[int(t)]
        print(psid, tmpl.format(person_id, "institute", institute_id))

    if department_id > 0:
        t = db.associate(sup_cur, person_id, department_id)
        tmpl = templates[int(t)]
        print(psid, tmpl.format(person_id, "department", department_id))

    if laboratory_id > 0:
        t = db.associate(sup_cur, person_id, laboratory_id)
        tmpl = templates[int(t)]
        print(psid, tmpl.format(person_id, "laboratory", laboratory_id))


def bad_email(email: str) -> bool:
    return ' ' in email


def error(*values, sep=' ', end='\n', flush=False) -> None:
    """
    Prints values to `stderr` instead of `stdout`.

    Equivalent to `print(*values, file=sys.stderr)`.
    """
    print(*values, sep=sep, end=end, file=sys.stderr, flush=flush)


def main():
    """Adds people and organizations to the mwb_supplemental database."""
    if len(sys.argv) < 2:
        error(__doc__)
        sys.exit(2)

    if sys.argv[1] in ["-h", "--help"]:
        print(__doc__)
        sys.exit()

    config_path = sys.argv[1]
    config = metab_import.get_config(config_path)

    mwb_client = mwb.Client(config.get("mwb_host"), config.get("mwb_port"))
    sup_conn: db.Connection = psycopg2.connect(
        host=config.get("sup_host"), dbname=config.get("sup_database"),
        user=config.get("sup_username"), password=config.get("sup_password"),
        port=config.get("sup_port"))

    embargoed_path = config.get("embargoed", "")
    embargoed: List[str] = []
    if embargoed_path:
        with open(embargoed_path) as f:
            embargoed = [line.strip() for line in f if line]

    with sup_conn:
        with sup_conn.cursor() as sup_cur:
            process_projects_and_studies(mwb_client, sup_cur, embargoed)
            add_developers(sup_cur)

    sup_conn.close()


def parse_author_list(xml: str) -> ET.ElementTree:
    file = io.StringIO(xml)
    data = ET.parse(file)
    author_list = data.findall("//Article/AuthorList/Author")
    return author_list


def process_projects_and_studies(mwb_client: mwb.Client,
                                 sup_cur: db.Cursor,
                                 embargoed: List[str]
                                 ) -> None:
    """
    Process all `project` and `study` records from Metabolomics Workbench.

    Parse names of people and organizations and then associate the right people
    with the right organizations.

    Note: embargoed studies are skipped.
    """
    records = mwb_client.fetch_names()

    for rec in records:
        if rec.pstype == mwb.STUDY and rec.psid in embargoed:
            continue  # Exclude embargoed studies.

        try:
            ppl = add_people(sup_cur, rec)
            orgs = add_organizations(sup_cur, rec)
        except AmbiguityError as e:
            error(rec.psid, type(e).__name__, e)
            continue

        if len(ppl) == 1:
            person_id = ppl[0]
            for (institute_id, dept_id, lab_id) in orgs:
                associate(sup_cur,
                          rec.psid, person_id, institute_id, dept_id, lab_id)
            continue

        if len(orgs) == 1:
            institute_id, dept_id, lab_id = orgs[0]
            for person_id in ppl:
                associate(sup_cur,
                          rec.psid, person_id, institute_id, dept_id, lab_id)
            continue

        if len(ppl) != len(orgs):
            error(rec.psid, f"cannot determine association between {len(ppl)}",
                  f"person(s) and {len(orgs)} organization(s).")
            continue

        for (person_id, org) in zip(ppl, orgs):
            institute_id, dept_id, lab_id = org
            associate(sup_cur,
                      rec.psid, person_id, institute_id, dept_id, lab_id)


if __name__ == "__main__":
    main()
