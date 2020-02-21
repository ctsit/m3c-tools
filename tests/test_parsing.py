import unittest

import metab_classes


class TestParse(unittest.TestCase):
    def test_parse(self):
        """ Confirm citations are being written in the correct format """
        xml = """
            <PubmedArticle>
            <MedlineCitation Status="PubMed-not-MEDLINE" Owner="NLM">
               <PMID Version="1">999999</PMID>
               <Article PubModel="Print">
                   <Journal>
                       <JournalIssue CitedMedium="Internet">
                           <Volume>1</Volume>
                           <Issue>5</Issue>
                           <PubDate>
                               <Year>2007</Year>
                               <Month>Sep</Month>
                           </PubDate>
                       </JournalIssue>
                       <Title>Journal of Example Science</Title>
                   </Journal>
                   <ArticleTitle>This is an example publication.</ArticleTitle>
                   <Pagination>
                       <MedlinePgn>100-105</MedlinePgn>
                   </Pagination>
                   <AuthorList CompleteYN="Y">
                       <Author ValidYN="Y">
                           <LastName>Smith</LastName>
                           <ForeName>John</ForeName>
                           <Initials>J</Initials>
                       </Author>
                       <Author ValidYN="Y">
                           <LastName>Doe</LastName>
                           <ForeName>Jane</ForeName>
                           <Initials>J</Initials>
                       </Author>
                       <Author ValidYN="Y">
                           <LastName>Hamill</LastName>
                           <ForeName>Mark</ForeName>
                           <Initials>M</Initials>
                       </Author>
                   </AuthorList>
               </Article>
            </MedlineCitation>
            <PubmedData>
                <ArticleIdList>
                    <ArticleId IdType="pubmed">999999</ArticleId>
                    <ArticleId IdType="pmc">PMC888888</ArticleId>
                    <ArticleId IdType="doi">12.3456/7890</ArticleId>
                </ArticleIdList>
            </PubmedData>
            </PubmedArticle>
            """

        pub = metab_classes.Publication.from_pubmed(xml.strip())
        expected_citation = (
            "Smith, J., Doe, J., Hamill, M. (2007). This is an example "
            "publication. Journal Of Example Science, 1(5), 100-105. "
            "doi:12.3456/7890")
        self.assertEqual(pub.pmid, '999999')
        self.assertEqual(pub.citation, expected_citation)

    def test_make_pub_handles_medlinedate(self):
        xml = """
            <PubmedArticle>
                <MedlineCitation Status="PubMed-not-MEDLINE" Owner="NLM">
                    <PMID Version="1">999999</PMID>
                    <Article PubModel="Print">
                        <Journal>
                            <JournalIssue CitedMedium="Internet">
                                <Volume>1</Volume>
                                <Issue>5</Issue>
                                <PubDate>
                                   <MedlineDate>1984 Dec-1985 Jan</MedlineDate>
                                </PubDate>
                            </JournalIssue>
                        </Journal>
                    </Article>
                </MedlineCitation>
            </PubmedArticle>
        """
        pub = metab_classes.Publication.from_pubmed(xml.strip())
        self.assertEqual(pub.published.year, 1984)
        self.assertEqual(pub.published.precision, 'year')
