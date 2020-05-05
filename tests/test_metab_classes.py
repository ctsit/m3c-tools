import unittest

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

        self.assertEquals('', actual[1])  # No summary
        self.assertListEqual(expected, actual[0])


if __name__ == "__main__":
    unittest.main()
