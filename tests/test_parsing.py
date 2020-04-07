import unittest

import metab_classes


class TestParsing(unittest.TestCase):
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

    def test_vernacular_title(self):
        xml = """
            <PubmedArticle>
                <MedlineCitation Status="MEDLINE" Owner="NLM">
                    <PMID Version="1">27759995</PMID>
                    <DateCompleted>
                        <Year>2017</Year>
                        <Month>10</Month>
                        <Day>16</Day>
                    </DateCompleted>
                    <DateRevised>
                        <Year>2018</Year>
                        <Month>10</Month>
                        <Day>23</Day>
                    </DateRevised>
                    <Article PubModel="Electronic">
                        <Journal>
                            <ISSN IssnType="Electronic">1699-5198</ISSN>
                            <JournalIssue CitedMedium="Internet">
                                <Volume>33</Volume>
                                <Issue>5</Issue>
                                <PubDate>
                                    <Year>2016</Year>
                                    <Month>Sep</Month>
                                    <Day>20</Day>
                                </PubDate>
                            </JournalIssue>
                            <Title>Nutricion hospitalaria</Title>
                            <ISOAbbreviation>Nutr Hosp</ISOAbbreviation>
                        </Journal>
                        <ArticleTitle/>
                        <Pagination>
                            <MedlinePgn>591</MedlinePgn>
                        </Pagination>
                        <ELocationID EIdType="doi" ValidYN="Y">10.20960/nh.591</ELocationID>
                        <Abstract>
                            <AbstractText>Introducción: la infertilidad es un problema global en aumento. Se estima que aproximadamente un 15% de las parejas en edad reproductiva tiene dificultades a la hora de concebir. De estas, alrededor de la mitad presentan uno o varios factores masculinos asociados a infertilidad o subfertilidad, aislados o en combinación con problemas de origen femenino. Durante la última década se ha empezado a estudiar la infertilidad desde una perspectiva multifactorial, considerando las interacciones y conexiones entre diferentes situaciones genéticas, epigenéticas, bioquímicas y fisiológicas del paciente.Objetivo: la presente revisión pretende describir mecanismos epigenéticos que pueden ser modulados mediante aspectos nutricionales y que están relacionados con la etiología de la infertilidad masculina y con la herencia transgeneracional de este fenotipo.Material y métodos: se ha realizado una extensa búsqueda de publicaciones científicas en las principales bases de datos electrónicas especializadas: NBCI, Elsevier, Scielo, Scirus y Science Direct.Resultados y conclusión: varios trabajos que muestran la importancia del estado nutricional en la fertilidad del hombre y, más específicamente, la capacidad de los componentes de la dieta para modificar los perfiles epigenéticos que no únicamente pueden afectar a su fertilidad, sino que también pueden ser transmitidos a la descendencia mediante lo que se ha denominado herencia transgeneracional, ocasionándoles problemas de salud diversos entre los que también se hallan problemas en la fertilidad.</AbstractText>
                        </Abstract>
                        <AuthorList CompleteYN="Y">
                            <Author ValidYN="Y">
                                <LastName>Oliver Bonet</LastName>
                                <ForeName>Maria</ForeName>
                                <Initials>M</Initials>
                                <AffiliationInfo>
                                    <Affiliation>Àrea de Ciències de la Salut. Institut Internacional de Postgrau de la Universitat Oberta de Catalunya (UOC). Barcelona. nuria.mach@jouy.inra.fr.</Affiliation>
                                </AffiliationInfo>
                            </Author>
                            <Author ValidYN="Y">
                                <LastName>Mach</LastName>
                                <ForeName>Núria</ForeName>
                                <Initials>N</Initials>
                            </Author>
                        </AuthorList>
                        <Language>spa</Language>
                        <PublicationTypeList>
                            <PublicationType UI="D016428">Journal Article</PublicationType>
                            <PublicationType UI="D016454">Review</PublicationType>
                        </PublicationTypeList>
                        <VernacularTitle>Factores nutricionales y no nutricionales pueden afectar la fertilidad masculina mediante mecanismos epigenéticos.</VernacularTitle>
                        <ArticleDate DateType="Electronic">
                            <Year>2016</Year>
                            <Month>09</Month>
                            <Day>20</Day>
                        </ArticleDate>
                    </Article>
                    <MedlineJournalInfo>
                        <Country>Spain</Country>
                        <MedlineTA>Nutr Hosp</MedlineTA>
                        <NlmUniqueID>9100365</NlmUniqueID>
                        <ISSNLinking>0212-1611</ISSNLinking>
                    </MedlineJournalInfo>
                    <CitationSubset>IM</CitationSubset>
                    <MeshHeadingList>
                        <MeshHeading>
                            <DescriptorName UI="D000818" MajorTopicYN="N">Animals</DescriptorName>
                        </MeshHeading>
                        <MeshHeading>
                            <DescriptorName UI="D006801" MajorTopicYN="N">Humans</DescriptorName>
                        </MeshHeading>
                        <MeshHeading>
                            <DescriptorName UI="D007248" MajorTopicYN="N">Infertility, Male</DescriptorName>
                            <QualifierName UI="Q000235" MajorTopicYN="Y">genetics</QualifierName>
                            <QualifierName UI="Q000503" MajorTopicYN="Y">physiopathology</QualifierName>
                        </MeshHeading>
                        <MeshHeading>
                            <DescriptorName UI="D008297" MajorTopicYN="N">Male</DescriptorName>
                        </MeshHeading>
                        <MeshHeading>
                            <DescriptorName UI="D009752" MajorTopicYN="Y">Nutritional Status</DescriptorName>
                        </MeshHeading>
                    </MeshHeadingList>
                    <KeywordList Owner="NOTNLM">
                        <Keyword MajorTopicYN="N">Infertilidad masculina.  Epigenética.  Nutrición.</Keyword>
                    </KeywordList>
                </MedlineCitation>
                <PubmedData>
                    <History>
                        <PubMedPubDate PubStatus="received">
                            <Year>2016</Year>
                            <Month>09</Month>
                            <Day>20</Day>
                        </PubMedPubDate>
                        <PubMedPubDate PubStatus="pubmed">
                            <Year>2016</Year>
                            <Month>10</Month>
                            <Day>21</Day>
                            <Hour>6</Hour>
                            <Minute>0</Minute>
                        </PubMedPubDate>
                        <PubMedPubDate PubStatus="medline">
                            <Year>2017</Year>
                            <Month>10</Month>
                            <Day>17</Day>
                            <Hour>6</Hour>
                            <Minute>0</Minute>
                        </PubMedPubDate>
                        <PubMedPubDate PubStatus="entrez">
                            <Year>2016</Year>
                            <Month>10</Month>
                            <Day>21</Day>
                            <Hour>6</Hour>
                            <Minute>0</Minute>
                        </PubMedPubDate>
                    </History>
                    <PublicationStatus>epublish</PublicationStatus>
                    <ArticleIdList>
                        <ArticleId IdType="pubmed">27759995</ArticleId>
                        <ArticleId IdType="doi">10.20960/nh.591</ArticleId>
                    </ArticleIdList>
                </PubmedData>
            </PubmedArticle>
        """

        pub = metab_classes.Publication.from_pubmed(xml.strip())
        self.assertEqual(pub.title,
                         "Factores nutricionales y no nutricionales pueden afectar la fertilidad masculina mediante mecanismos epigenéticos.")


if __name__ == "__main__":
    unittest.main()
