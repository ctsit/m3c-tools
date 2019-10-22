import tempfile
import unittest

from Bio import Entrez

from metab_classes import Publication
from metab_pub_ingest import Citation
from metab_pub_ingest import fill_pub


class TestParse(unittest.TestCase):
    def test_parse(self):
        """ Confirm citations are being written in the correct format """
        results = self.create_citation()
        citation = Citation(results['PubmedArticle'][0])
        pub = Publication()
        fill_pub(pub, citation)

        expected_citation = "Smith, J., Doe, J., Hamill, M. (2007). This is an example publication. Journal Of Example Science, 1(5), 100-105. doi:12.3456/7890"
        self.assertEqual(pub.pmid, '999999')
        self.assertEqual(pub.citation, expected_citation)

    def create_citation(self):
        with tempfile.TemporaryFile() as handle:
            handle.write(b'''\
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
