import unittest

from m3c import classes
import m3c.classes as metab_classes


class TestStudy(unittest.TestCase):
    def test_get_triples(self):
        s = metab_classes.Study(
            study_id="STUDY_ID",
            study_title="Study Title",
            study_type="Study Type",
            summary="",
            submit_date="2020-02-07",
            project_id="PROJECT_ID"
        )
        s.institutes.append("INSTITUTE_ID")
        s.institutes.append("OTHER_INSTITUTE_ID")

        actual = s.get_triples(namespace="http://example.com/i/")

        expected = [
            '<http://example.com/i/STUDY_ID> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Study>',
            '<http://example.com/i/STUDY_ID> <http://www.w3.org/2000/01/rdf-schema#label> "Study Title"^^<http://www.w3.org/2001/XMLSchema#string>',
            '<http://example.com/i/STUDY_ID> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#studyId> "STUDY_ID"^^<http://www.w3.org/2001/XMLSchema#string>',
            '<http://example.com/i/STUDY_ID> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#workbenchLink> "https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=STUDY_ID"^^<http://www.w3.org/2001/XMLSchema#string>',
            '<http://example.com/i/STUDY_ID> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#studyType> "Study Type"^^<http://www.w3.org/2001/XMLSchema#string>',
            '<http://example.com/i/STUDY_ID> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#submitted> "2020-02-07"^^<http://www.w3.org/2001/XMLSchema#dateTime>',
            '<http://example.com/i/STUDY_ID> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#inCollection> <http://example.com/i/PROJECT_ID>',
            '<http://example.com/i/PROJECT_ID> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#collectionFor> <http://example.com/i/STUDY_ID>',
            '<http://example.com/i/STUDY_ID> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#directedBy> <http://example.com/i/oINSTITUTE_ID>',
            '<http://example.com/i/oINSTITUTE_ID> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#directs> <http://example.com/i/STUDY_ID>',
            '<http://example.com/i/STUDY_ID> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#directedBy> <http://example.com/i/oOTHER_INSTITUTE_ID>',
            '<http://example.com/i/oOTHER_INSTITUTE_ID> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#directs> <http://example.com/i/STUDY_ID>',
        ]

        self.assertEqual('', actual[1])  # No summary
        self.assertListEqual(expected, actual[0])


class TestProject(unittest.TestCase):
    def test_escapes_doublequote(self):
        summary_with_doublequotes = r'"Evaluating lipid mediator"'
        expected = r' "\"Evaluating lipid mediator\""^^'
        project = metab_classes.Project(
            project_id="PROJECT_ID",
            project_type="",
            project_title="",
            summary=summary_with_doublequotes,
            doi="",
            funding_source="",
        )
        _, actual = project.get_triples("x://test/")
        self.assertIn(expected, actual)

    def test_escapes_newline(self):
        summary_with_newline = "datatrack 1810\nPS project 1"
        expected = ' "datatrack 1810\\nPS project 1"^^'
        project = metab_classes.Project(
            project_id="PROJECT_ID",
            project_type="",
            project_title="",
            summary=summary_with_newline,
            doi="",
            funding_source="",
        )
        _, actual = project.get_triples("x://test/")
        self.assertIn(expected, actual)

    def test_escapes_backslash(self):
        summary_with_backslash = r"In a world where 4\% percent of..."
        expected = r"In a world where 4\\% percent of..."
        project = metab_classes.Project(
            project_id="PROJECT_ID",
            project_type="",
            project_title="",
            summary=summary_with_backslash,
            doi="",
            funding_source="",
        )
        _, actual = project.get_triples("x://test/")
        self.assertIn(expected, actual)


class TestTool(unittest.TestCase):
    def test_match_authors(self):
        pid = 780
        igor = classes.Person(person_id=str(pid),
                              first_name="John",
                              last_name="Doe",
                              display_name="John Doe")
        alias = "John P Doe"

        aliases = [
            (pid, alias)
        ]
        people = {
            pid: igor,
        }

        def find_person(name: str) -> int:
            # This function (and therefore the enitre test) might seem trivial
            # but I'm leaving it to document the expected behavior of
            # find_person and match_authors and prevent a regression.
            # The expected behavior was NOT the actual behavior before SEP-192
            # was fixed.
            for (pid, display_name) in people.items():
                if display_name == name:
                    return pid
            for (pid, alias) in aliases:
                if alias == name:
                    return pid
            return 0

        info = {
            "name": "TestTool",
            "description": "Fake tool for unit-testing",
            "url": "",
            "authors": [{
                "name": alias,
            }],
        }

        tool = classes.Tool(tool_id="TestTool", data=info)

        unmatched = tool.match_authors(find_person, "x://test/")

        self.assertEqual(0, len(unmatched))


if __name__ == "__main__":
    unittest.main()
