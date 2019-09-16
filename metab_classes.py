import json
import re
import typing


class Project(object):
    def __init__(self):
        self.uri = None
        self.project_id = None
        self.project_type = None
        self.summary = None
        self.doi = None
        self.funding_source = None
        self.institute_uri = None
        self.department_uri = None
        self.lab_uri = None
        self.last_name = None
        self.first_name = None
        self.pi_uri = None

    def get_triples(self):
        rdf = []
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Project>".format(self.uri))
        rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.project_title))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#projectId> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.project_id))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#workbenchLink> \"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID={}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.project_id))
        if self.project_type:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#projectType> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.project_type))
        if self.doi:
            rdf.append("<{}> <http://purl.org/ontology/bibo/doi> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.doi))
        if self.institute_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#managedBy> <{}>".format(self.uri, self.institute_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#manages> <{}>".format(self.institute_uri, self.uri))
        if self.department_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#managedBy> <{}>".format(self.uri, self.department_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#manages> <{}>".format(self.department_uri, self.uri))
        if self.lab_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#managedBy> <{}>".format(self.uri, self.lab_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#manages> <{}>".format(self.lab_uri, self.uri))
        if self.pi_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#hasPI> <{}>".format(self.uri, self.pi_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#isPIFor> <{}>".format(self.pi_uri, self.uri))
        if self.summary:
            summary_line = "<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#summary> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.summary)
        else:
            summary_line = None
        return rdf, summary_line


class Study(object):
    def __init__(self):
        self.uri = None
        self.study_id = None
        self.study_title = None
        self.study_type = None
        self.summary = None
        self.submit_date = None
        self.institute_uri = None
        self.department_uri = None
        self.lab_uri = None
        self.last_name = None
        self.first_name = None
        self.runner_uri = None
        self.project_id = None

    def get_triples(self, project_uri=None):
        rdf = []
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Study>".format(self.uri))
        rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.study_title))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#studyId> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.study_id))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#workbenchLink> \"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID={}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.study_id))
        if self.study_type:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#studyType> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self. study_type))
        if self.submit_date:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#submitted> \"{}\"^^<http://www.w3.org/2001/XMLSchema#dateTime>".format(self.uri, self.submit_date))
        if project_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#inCollection> <{}>".format(self.uri, project_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#collectionFor> <{}>".format(project_uri, self.uri))
        if self.institute_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#managedBy> <{}>".format(self.uri, self.institute_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#manages> <{}>".format(self.institute_uri, self.uri))
        if self.department_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#managedBy> <{}>".format(self.uri, self.department_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#manages> <{}>".format(self.department_uri, self.uri))
        if self.lab_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#managedBy> <{}>".format(self.uri, self.lab_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#manages> <{}>".format(self.lab_uri, self.uri))
        if self.runner_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#runBy> <{}>".format(self.uri, self.runner_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#runnerOf> <{}>".format(self.runner_uri, self.uri))
        if self.summary:
            summary_line = "<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#summary> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.summary)
        else:
            summary_line = None
        return rdf, summary_line


class Dataset(object):
    def __init__(self):
        self.uri = None
        self.mb_sample_id = None
        self.subject_species = None
        self.study_id = None

    def get_triples(self, study_uri=None):
        rdf = []
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Dataset>".format(self.uri))
        # rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.mb_sample_id))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#sampleId> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.mb_sample_id))
        if self.subject_species:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#subjectSpecies> \"{}\"".format(self.uri, self.subject_species))
        if study_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#dataFor> <{}>".format(self.uri, study_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#developedFrom> <{}>".format(study_uri, self.uri))
        return rdf


class Person(object):
    def __init__(self):
        self.uri = None
        self.person_id = None
        self.first_name = None
        self.last_name = None
        self.display_name = None
        self.email = None
        self.phone = None

    def make_display_name(self):
        self.display_name = self.first_name + ' ' + self.last_name

    def get_triples(self):
        rdf = []
        vcard_uri = self.uri + "vcard"
        name_uri = vcard_uri + "name"
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://xmlns.com/foaf/0.1/Person>".format(self.uri))
        rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.display_name).replace('\n', ''))
        rdf.append("<{}> <http://purl.obolibrary.org/obo/ARG_2000028> <{}>".format(self.uri, vcard_uri))
        rdf.append("<{}> <http://purl.obolibrary.org/obo/ARG_2000029> <{}>".format(vcard_uri, self.uri))
        rdf.append("<{}> <http://www.w3.org/2006/vcard/ns#hasName> <{}>".format(vcard_uri, name_uri))
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2006/vcard/ns#Name>".format(name_uri))
        rdf.append("<{}> <http://www.w3.org/2006/vcard/ns#familyName> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(vcard_uri, self.last_name).replace('\n', ''))
        rdf.append("<{}> <http://www.w3.org/2006/vcard/ns#givenName> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(vcard_uri, self.first_name.replace('\n', '')))
        if self.email:
            email_uri = vcard_uri + 'email'
            rdf.append("<{}> <http://www.w3.org/2006/vcard/ns#hasEmail>  <{}>".format(vcard_uri, email_uri))
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2006/vcard/ns#Email>".format(email_uri))
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2006/vcard/ns#Work>".format(email_uri))
            rdf.append("<{}> <http://www.w3.org/2006/vcard/ns#email> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(email_uri, self.email))
        if self.phone:
            phone_uri = vcard_uri + 'phone'
            rdf.append("<{}> <http://www.w3.org/2006/vcard/ns#hasTelephone>  <{}>".format(vcard_uri, phone_uri))
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2006/vcard/ns#Telephone>".format(phone_uri))
            rdf.append("<{}> <http://www.w3.org/2006/vcard/ns#telephone> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(phone_uri, self.phone))
        return rdf


class Organization(object):
    def __init__(self):
        self.uri = None
        self.org_id = None
        self.name = None
        self.type = None
        self.parent_id = None
        self.parent_uri = None

    def get_triples(self):
        rdf = []
        if self.type == "institute":
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vivoweb.org/ontology/core#Institute>".format(self.uri))
        if self.type == "department":
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vivoweb.org/ontology/core#Department>".format(self.uri))
        if self.type == "laboratory":
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vivoweb.org/ontology/core#Laboratory>".format(self.uri))
        rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.name.replace('\n', '')))
        if self.parent_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#hasParent> <{}>".format(self.uri, self.parent_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#parentOf> <{}>".format(self.parent_uri, self.uri))
        return rdf

    def add_person(self, person_uri):
        rdf = []
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#associatedWith> <{}>".format(person_uri, self.uri))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#associationFor> <{}>".format(self.uri, person_uri))
        return rdf


class Tool(object):
    class License:
        def __init__(self, kind: typing.Text,
                     url: typing.Optional[typing.Text] = ''):
            self.kind = kind.strip()
            self.url = url.strip()

    class Author:
        def __init__(self, name: str, email: typing.Optional[str] = '',
                     uri: typing.Optional[str] = ''):
            self.name = name.strip()
            self.email = email.strip()
            self.uri = uri.strip()

    def __init__(self, tool_id: typing.Text, data: dict):
        self.tool_id: typing.Text = tool_id.strip()
        self.name: typing.Text = data['name'].replace('\n', ' ').strip()
        self.description: typing.Text = data['description'].strip()
        self.url: typing.Text = data['url'].strip()
        self.authors: typing.List[Tool.Author] = []
        for author in data['authors']:
            self.authors.append(Tool.Author(**author))
        license = data['license']
        self.license: Tool.License = Tool.License(**license)
        self.tags: typing.List[typing.Text] = data.get('tags', [])

    def uri(self, namespace: typing.Text) -> typing.Text:
        encoded = self.tool_id
        encoded = encoded.replace('_', '__')
        encoded = encoded.replace('-', '_d')
        encoded = encoded.replace('/', '_s')

        contains_unhandled_char = re.search('[^A-Za-z0-9._/]', encoded)
        if contains_unhandled_char:
            raise Exception(
                "Unhandled character in tool's ID: %s" % self.tool_id)

        return namespace + 't' + encoded

    def get_triples(self, namespace: typing.Text) -> typing.List[typing.Text]:
        m3c = 'http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#'

        rdf = []
        uri = self.uri(namespace)

        rdf.append('<{uri}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <{m3c}Tool>'
                   .format(uri=uri, m3c=m3c))
        rdf.append('<{uri}> <http://www.w3.org/2000/01/rdf-schema#label> {name}'
                   .format(uri=uri, name=escape(self.name)))
        rdf.append('<{uri}> <{m3c}summary> {desc}'
                   .format(uri=uri, m3c=m3c, desc=escape(self.description)))

        rdf.append('<{uri}> <{m3c}homepage> {desc}'
                   .format(uri=uri, m3c=m3c, desc=escape(self.url)))

        if not self.license or not self.license.kind or not self.license.url:
            raise Exception('Bad license for tool: ' + self.tool_id)

        rdf.append('<{uri}> <{m3c}licenseType> {kind}'
                   .format(uri=uri, m3c=m3c, kind=escape(self.license.kind)))
        rdf.append('<{uri}> <{m3c}licenseUrl> {link}'
                   .format(uri=uri, m3c=m3c, link=escape(self.license.url)))

        for author in self.authors:
            if not author.uri:
                raise Exception('Unknown author "%s" for tool: %s' %
                                (author.name, self.tool_id))
            rdf.append("<{uri}> <{m3c}developedBy> <{author.uri}>"
                       .format(uri=uri, m3c=m3c, author=author))
            rdf.append("<{author.uri}> <{m3c}developerOf> <{uri}>"
                       .format(uri=uri, m3c=m3c, author=author))

        for tag in self.tags:
            tag = tag.strip().lower()
            rdf.append("<{uri}> <{m3c}tag> {tag}"
                       .format(uri=uri, m3c=m3c, tag=escape(tag)))

        return rdf

    def match_authors(self, people: typing.Dict[int, Person]):
        for author in self.authors:
            for person in people.values():
                if person.display_name == author.name:
                    author.uri = person.uri
                    break
            if not author.uri:
                raise Exception('Unknown author "%s" for tool: %s' %
                                (author.name, self.tool_id))


def escape(text):
    return json.dumps(text)
