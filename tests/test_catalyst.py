import unittest
import xml.etree.ElementTree as ET

import catalyst
from metab_classes import Person


class TestCatalyst(unittest.TestCase):
    def test_parse_pmids_no_pmids(self):
        empty_pmidlist = """
                        <PMIDList>
                        </PMIDList>"""
        pmids = catalyst.parse_catalyst_pmids(empty_pmidlist)
        self.assertListEqual(pmids, [])

    def test_parse_pmids_list_pmids(self):
        pmidlist = """
                        <PMIDList>
                        <PMID>1</PMID>
                        <PMID>2</PMID>
                        <PMID>3</PMID>
                        </PMIDList>"""
        pmids = catalyst.parse_catalyst_pmids(pmidlist)
        self.assertListEqual(pmids, ['1', '2', '3'])

    def test_parse_pmids_none_pmids(self):
        pmidlist = None
        pmids = catalyst.parse_catalyst_pmids(pmidlist)
        self.assertListEqual(pmids, [])

    def test_build_catalyst_xml_none_affiliations(self):
        with self.assertRaises(AssertionError):
            catalyst.build_catalyst_xml(
                Person('1', 'f', 'l', 'd', 'e', 'p'), None, ['1'], ['2'])

    def test_build_catalyst_xml_empty_affiliations(self):
        with self.assertRaises(AssertionError):
            catalyst.build_catalyst_xml(
                Person('1', 'f', 'l', 'd', 'e', 'p'), [], ['1'], ['2'])

    def test_build_catalyst_xml_none_include(self):
        with self.assertRaises(AssertionError):
            catalyst.build_catalyst_xml(
                Person('1', 'f', 'l', 'd', 'e', 'p'), ['a'], None, ['2'])

    def test_build_catalyst_xml_empty_include(self):
        with self.assertRaises(AssertionError):
            catalyst.build_catalyst_xml(
                Person('1', 'f', 'l', 'd', 'e', 'p'), ['a'], [], ['2'])

    def test_build_catalyst_xml_none_exclude(self):
        try:
            catalyst.build_catalyst_xml(
                Person('1', 'f', 'l', 'd', 'e', 'p'), ['a'], ['1'], None)
        except Exception:
            self.fail('None exclude list caused exception')

    def test_build_catalyst_xml_empty_exclude(self):
        try:
            catalyst.build_catalyst_xml(
                Person('1', 'f', 'l', 'd', 'e', 'p'), ['a'], ['1'], [])
        except Exception:
            self.fail('Empty exclude list caused exception')

    def test_build_catalyst_xml_match_example(self):
        xml = """
        <FindPMIDs>
            <Name>
                <First>Griffin</First>
                <Last>Weber</Last>
            </Name>
            <EmailList>
                <email>weber@hms.harvard.edu</email>
            </EmailList>
            <AffiliationList>
                <Affiliation>%Brigham%Women%</Affiliation>
                <Affiliation>%@hms.harvard.edu%</Affiliation>
            </AffiliationList>
            <LocalDuplicateNames>1</LocalDuplicateNames>
            <RequireFirstName>false</RequireFirstName>
            <MatchThreshold>0.98</MatchThreshold>
            <PMIDAddList>
                <PMID>11707567</PMID>
            </PMIDAddList>
            <PMIDExcludeList>
                <PMID>19648504</PMID>
            </PMIDExcludeList>
        </FindPMIDs>
        """.strip().replace('\n', '').replace('\t', '').replace(' ', '')
        xml = ET.fromstring(xml)
        xml = ET.tostring(xml)
        actual_xml = catalyst.build_catalyst_xml(
            Person('1', 'Griffin', 'Weber', 'Griffin Webber',
                   'weber@hms.harvard.edu', '1234567'),
            ['Brigham%Women', '@hms.harvard.edu'],
            ['11707567'],
            ['19648504'])
        self.assertEqual(xml, actual_xml)
