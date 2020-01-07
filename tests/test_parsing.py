import io
import unittest

from Bio import Entrez

import metab_pub_ingest


class TestParse(unittest.TestCase):
    def test_parse(self):
        """ Confirm citations are being written in the correct format """
        results = self.create_citation()
        citation = metab_pub_ingest.Citation(results['PubmedArticle'][0])
        pub = metab_pub_ingest.make_pub(citation)

        expected_citation = "Smith, J., Doe, J., Hamill, M. (2007). This is an example publication. Journal Of Example Science, 1(5), 100-105. doi:12.3456/7890"
        self.assertEqual(pub.pmid, '999999')
        self.assertEqual(pub.citation, expected_citation)

    def test_make_pub_handles_medlinedate(self):
        xml = io.StringIO('''
            <?xml version="1.0"?>
            <!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st January 2019//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_190101.dtd">
            <PubmedArticleSet>
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
            </PubmedArticleSet>
        '''.strip())
        pubmed_articles = Entrez.read(xml)
        article = pubmed_articles['PubmedArticle'][0]
        citation = metab_pub_ingest.Citation(article)
        pub = metab_pub_ingest.make_pub(citation)
        self.assertEqual(pub.published.year, 1984)
        self.assertEqual(pub.published.precision, 'year')

    def create_citation(self):
        handle = io.StringIO()
        handle.write('''\
<?xml version="1.0"?>
<!DOCTYPE PubmedArticleSet PUBLIC "-//NLM//DTD PubMedArticle, 1st January 2019//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_190101.dtd">
<PubmedArticleSet>
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
</PubmedArticleSet>''')

        handle.seek(0)
        results = Entrez.read(handle)
        return results
