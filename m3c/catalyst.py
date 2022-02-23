"""
Interface to the Harvard PubMed Disambiguation Tool's Profiles Catalyst service

See http://profiles.catalyst.harvard.edu/docs/ProfilesRNS_DisambiguationEngine.pdf
"""

from typing import List

import sys
import traceback
import xml.etree.ElementTree as ET
import os.path
import requests

from m3c.classes import Person
from m3c.logger import Logger

ENDPOINT = "http://profiles.catalyst.harvard.edu/services/GetPMIDs/default.asp"


def build_catalyst_xml(
    person: Person,
    affiliations: List[str],
    include_pmids: List[str],
    exclude_pmids: List[str]
) -> bytes:
    """
    Build the XML request payload.

    `affiliations` and `include_pmids` must contain at least 1 entry each.
    """
    assert affiliations and len(affiliations) > 0
    assert include_pmids and len(include_pmids) > 0

    root = ET.Element("FindPMIDs")
    name = ET.SubElement(root, "Name")
    first_name = ET.SubElement(name, "First")
    first_name.text = person.first_name

    last_name = ET.SubElement(name, "Last")
    last_name.text = person.last_name

    email_list = ET.SubElement(root, "EmailList")
    if person.email is not None or person.email != "":
        email = ET.SubElement(email_list, "email")
        email.text = person.email

    affiliation_list = ET.SubElement(root, "AffiliationList")
    for affiliation in affiliations:
        affil = ET.SubElement(affiliation_list, "Affiliation")
        affil.text = f"%{affiliation}%"

    local_dup_names = ET.SubElement(root, "LocalDuplicateNames")
    local_dup_names.text = "1"

    require_first_name = ET.SubElement(root, "RequireFirstName")
    require_first_name.text = "false"

    match_threshold = ET.SubElement(root, "MatchThreshold")
    match_threshold.text = "0.98"

    add_list = ET.SubElement(root, "PMIDAddList")
    for pmid in include_pmids:
        pmid_elm = ET.SubElement(add_list, "PMID")
        pmid_elm.text = pmid

    rem_list = ET.SubElement(root, "PMIDExcludeList")
    if exclude_pmids:
        for pmid in exclude_pmids:
            pmid_elm = ET.SubElement(rem_list, "PMID")
            pmid_elm.text = pmid

    xml: bytes = ET.tostring(root)
    return xml


def parse_catalyst_pmids(catalyst_xml: str) -> List[str]:
    """
    Parse out the PMIDs from the Catalyst results XML.

    This XML looks like:
    ```
    <PMIDList>
        <PMID>11707567</PMID>
        <PMID>11788827</PMID>
        <PMID>11815958</PMID>
        <PMID>12209713</PMID>
        <PMID>12463949</PMID>
        <PMID>12659816</PMID>
        <PMID>15219292</PMID>
        <PMID>15226823</PMID>
        <PMID>16359929</PMID>
        <PMID>18541841</PMID>
        <PMID>19567788</PMID>
        <PMID>20190053</PMID>
    </PMIDList>
    ```
    """
    if not catalyst_xml:
        return []
    try:
        root = ET.fromstring(catalyst_xml)
        return [pmid.text if pmid.text else "" for pmid in root]
    except Exception:
        traceback.print_exc()
        return []


def fetch_ids(person: Person, affiliations: List[str],
              include_pmids: List[str], exclude_pmids: List[str]) \
        -> List[str]:
    """
    Gets the disambiguated PubMed publications.

    Must pass affiliations and include_pmids.
    """
    assert len(affiliations) > 0 and len(include_pmids) > 0

    payload_xml = build_catalyst_xml(
        person, affiliations, include_pmids, exclude_pmids)
    headers = {
        "Content-Type": "text/xml"
    }
    resp = requests.post(ENDPOINT, data=payload_xml, headers=headers)
    if resp.status_code != 200:
        log_path = os.path.join('', 'log.txt')
        log = Logger(log_path)

        log("Unexpected response from Catalyst", resp.status_code,
              file=sys.stderr)
        print("Unexpected response from Catalyst", resp.status_code,
              file=sys.stderr)

        return []

    return parse_catalyst_pmids(resp.text)
