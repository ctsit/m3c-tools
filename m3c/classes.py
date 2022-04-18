from typing import Callable, Dict, List, Optional, Set, Text

import io
import os
import re
import textwrap

from Bio import Entrez
from m3c.logger import Logger


MONTHS: List[str] = 'Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec'.split()


class Citation(object):
    def __init__(self, data):
        self.data = data

    def check_key(self, paths, data=None):
        if not data:
            data = self.data
        if paths[0] in data:
            trail = data[paths[0]]
            if len(paths) > 1:
                trail = self.check_key(paths[1:], trail)
            return trail
        else:
            return ''


class Project(object):
    def __init__(self, project_id: str, project_type: str, project_title: str,
                 summary: str, doi: str, funding_source: str):
        assert project_id
        self.project_id = project_id.strip()
        self.project_type = project_type
        self.project_title = project_title
        self.summary = summary
        self.doi = doi
        self.funding_source = funding_source
        self.pi: List[str] = []
        self.institutes: List[str] = []
        self.departments: List[str] = []
        self.labs: List[str] = []

    def get_triples(self, namespace: str):
        uri = Project.uri(namespace, self.project_id)
        rdf: List[str] = []
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Project>".format(uri))
        rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, escape(self.project_title)))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#projectId> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, self.project_id))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#workbenchLink> \"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID={}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, self.project_id))
        if self.project_type:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#projectType> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, escape(self.project_type)))
        if self.doi:
            rdf.append("<{}> <http://purl.org/ontology/bibo/doi> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, escape(self.doi)))
        if self.institutes:
            for institute in self.institutes:
                institute_uri = Organization.uri(namespace, institute)
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#managedBy> <{}>".format(uri, institute_uri))
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#manages> <{}>".format(institute_uri, uri))
        if self.departments:
            for department in self.departments:
                dept_id = Organization.uri(namespace, department)
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#managedBy> <{}>".format(uri, dept_id))
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#manages> <{}>".format(dept_id, uri))
        if self.labs:
            for lab in self.labs:
                lab_uri = Organization.uri(namespace, lab)
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#managedBy> <{}>".format(uri, lab_uri))
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#manages> <{}>".format(lab_uri, uri))
        if self.pi:
            for person in self.pi:
                pi_uri = Person.uri(namespace, person)
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#hasPI> <{}>".format(uri, pi_uri))
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#isPIFor> <{}>".format(pi_uri, uri))
        if self.summary:
            summary_line = "<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#summary> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, escape(self.summary))
        else:
            summary_line = ""
        return rdf, summary_line

    @staticmethod
    def uri(namespace: str, project_id: str) -> str:
        return f"{namespace}{project_id}"


class Study(object):
    def __init__(self, study_id: str, study_title: str, study_type: str,
                 summary: str, submit_date: str, project_id: str):
        assert study_id

        self.study_id = study_id.strip()
        self.study_title = study_title
        self.study_type = study_type
        self.summary = summary
        self.submit_date = submit_date
        self.project_id = project_id

        self.runner: List[str] = []
        self.institutes: List[str] = []
        self.departments: List[str] = []
        self.labs: List[str] = []

        self.subject_species: List[str] = []

    def get_triples(self, namespace: str):
        uri = Study.uri(namespace, self.study_id)
        rdf = []
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Study>".format(uri))
        rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, escape(self.study_title)))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#studyId> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, self.study_id))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#workbenchLink> \"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID={}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, self.study_id))
        if self.study_type:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#studyType> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, escape(self.study_type)))
        if self.submit_date:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#submitted> \"{}\"^^<http://www.w3.org/2001/XMLSchema#dateTime>".format(uri, self.submit_date))
        if self.project_id:
            project_uri = Project.uri(namespace, self.project_id)
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#inCollection> <{}>".format(uri, project_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#collectionFor> <{}>".format(project_uri, uri))
        if self.institutes:
            for institute in self.institutes:
                institute_uri = Organization.uri(namespace, institute)
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#directedBy> <{}>".format(uri, institute_uri))
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#directs> <{}>".format(institute_uri, uri))
        if self.departments:
            for department in self.departments:
                dept_uri = Organization.uri(namespace, department)
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#directedBy> <{}>".format(uri, dept_uri))
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#directs> <{}>".format(dept_uri, uri))
        if self.labs:
            for lab in self.labs:
                lab_uri = Organization.uri(namespace, lab)
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#directedBy> <{}>".format(uri, lab_uri))
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#directs> <{}>".format(lab_uri, uri))
        if self.runner:
            for person in self.runner:
                runner_uri = Person.uri(namespace, person)
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#runBy> <{}>".format(uri, runner_uri))
                rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#runnerOf> <{}>".format(runner_uri, uri))
        if self.summary:
            summary_line = "<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#summary> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, escape(self.summary))
        else:
            summary_line = ""
        return rdf, summary_line

    def get_species_triples(self, namespace: str):
        uri = f"{namespace}{self.study_id}"
        rdf = []
        for species in self.subject_species:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#subjectSpecies> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, species))
        return rdf

    @staticmethod
    def uri(namespace: str, study_id: str) -> str:
        return f"{namespace}{study_id}"


class Dataset(object):
    def __init__(self):
        self.uri = None
        self.mb_sample_id = None
        self.subject_species = None
        self.study_id = None

    def get_triples(self, study_uri=None):
        rdf = []
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Dataset>".format(self.uri))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#sampleId> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.mb_sample_id))
        if self.subject_species:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#subjectSpecies> \"{}\"".format(self.uri, self.subject_species))
        if study_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#dataFor> <{}>".format(self.uri, study_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#developedFrom> <{}>".format(study_uri, self.uri))
        return rdf


class Person(object):
    def __init__(self, person_id: str, first_name: str, last_name: str,
                 display_name="", email="", phone="", withheld=False,
                 overview=""):
        assert person_id and first_name and last_name

        self.person_id = person_id
        self.first_name = first_name
        self.last_name = last_name
        self.email = email
        self.phone = phone
        self.display_name = display_name
        self.withheld = withheld
        self.overview = overview

        if not self.display_name:
            self.display_name = f"{self.first_name} {self.last_name}"

    def get_triples(self, namespace: str):
        if self.withheld:
            return []
        uri = Person.uri(namespace, self.person_id)
        rdf = []
        vcard_uri = uri + "vcard"
        name_uri = vcard_uri + "name"
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://xmlns.com/foaf/0.1/Person>".format(uri))
        rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, self.display_name).replace('\n', ''))
        rdf.append("<{}> <http://purl.obolibrary.org/obo/ARG_2000028> <{}>".format(uri, vcard_uri))
        rdf.append("<{}> <http://purl.obolibrary.org/obo/ARG_2000029> <{}>".format(vcard_uri, uri))
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
        if self.overview:
            rdf.append(f'<{uri}> <http://vivoweb.org/ontology/core#overview> "{self.overview}"^^<http://www.w3.org/2001/XMLSchema#string>')
        return rdf

    @staticmethod
    def n_number(person_id: str) -> str:
        return f"p{person_id}"

    @staticmethod
    def uri(namespace: str, person_id: str) -> str:
        return f"{namespace}{Person.n_number(person_id)}"


class Photo(object):
    def __init__(self, file_storage_root: str, person_id: str, extension: str,
                 file_storage_alias: str = "b"):
        self.root = file_storage_root
        self.person_id = person_id
        assert int(self.person_id)
        self.alias = file_storage_alias

        extension = extension.lower()
        assert extension in ("png", "jpeg", "jpg")
        self.extension = "jpg"
        self.mimetype = "image/jpeg"
        if extension == "png":
            self.extension = "png"
            self.mimetype = "image/png"

    def download_url(self) -> str:
        return f"/file/{Person.n_number(self.person_id)}pic/{self.filename()}"

    def filename(self) -> str:
        return f"photo.{self.extension}"

    def get_triples(self, namespace: Text) -> List[Text]:
        person_uri = Person.uri(namespace, self.person_id)
        person = f"<{person_uri}>"
        image = f"<{person_uri}photo>"
        thumb = f"<{person_uri}thumb>"
        image_dl = f"<{person_uri}pic>"
        thumb_dl = f"<{person_uri}tn>"

        rdf = []
        rdf.append(f"{person} <http://vitro.mannlib.cornell.edu/ns/vitro/public#mainImage> {image}")

        rdf.append(f"{image} <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vitro.mannlib.cornell.edu/ns/vitro/public#File>")
        rdf.append(f"{image} <http://vitro.mannlib.cornell.edu/ns/vitro/public#downloadLocation> {image_dl}")
        rdf.append(f'{image} <http://vitro.mannlib.cornell.edu/ns/vitro/public#filename> "{self.filename()}"')
        rdf.append(f'{image} <http://vitro.mannlib.cornell.edu/ns/vitro/public#mimeType> "{self.mimetype}"')
        rdf.append(f"{image} <http://vitro.mannlib.cornell.edu/ns/vitro/public#thumbnailImage> {thumb}")

        rdf.append(f"{image_dl} <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vitro.mannlib.cornell.edu/ns/vitro/public#FileByteStream>")
        rdf.append(f'{image_dl} <http://vitro.mannlib.cornell.edu/ns/vitro/public#directDownloadUrl> "{self.download_url()}"')

        rdf.append(f"{thumb} <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vitro.mannlib.cornell.edu/ns/vitro/public#File>")
        rdf.append(f"{thumb} <http://vitro.mannlib.cornell.edu/ns/vitro/public#downloadLocation> {thumb_dl}")
        rdf.append(f'{thumb} <http://vitro.mannlib.cornell.edu/ns/vitro/public#filename> "{self.filename()}"')
        rdf.append(f'{thumb} <http://vitro.mannlib.cornell.edu/ns/vitro/public#mimeType> "{self.mimetype}"')

        # TODO: actually generate a thumbnail instead of using the full photo
        rdf.append(f"{thumb_dl} <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vitro.mannlib.cornell.edu/ns/vitro/public#FileByteStream>")
        rdf.append(f'{thumb_dl} <http://vitro.mannlib.cornell.edu/ns/vitro/public#directDownloadUrl> "{self.download_url()}"')

        return rdf

    def path(self) -> str:
        """Get the directory path for the specified person with `person_id`."""
        # "b~" is shorthand for https://vivo.metabolomics.info/individual/
        fullpath = f"{self.alias}~{Person.n_number(self.person_id)}pic"
        # VIVO expects each directory to be no longer than 3 characters.
        # See: https://wiki.duraspace.org/display/VIVODOC110x/Image+storage
        split = textwrap.wrap(fullpath, 3)
        path = os.path.join(self.root, *split)
        return path


class Organization(object):
    def __init__(self, org_id: str, name: str, type: str, parent_id: str):
        assert org_id
        self.org_id = org_id
        self.name = name
        self.type = type
        self.parent_id = parent_id

    def get_triples(self, namespace: str):
        uri = Organization.uri(namespace, self.org_id)
        rdf = []
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://xmlns.com/foaf/0.1/Organization>".format(uri))
        if self.type == "institute":
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vivoweb.org/ontology/core#Institute>".format(uri))
        if self.type == "department":
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vivoweb.org/ontology/core#Department>".format(uri))
        if self.type == "laboratory":
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vivoweb.org/ontology/core#Laboratory>".format(uri))
        rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, self.name.replace('\n', '')))
        if self.parent_id:
            parent_uri = Organization.uri(namespace, self.parent_id)
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#hasParent> <{}>".format(uri, parent_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#parentOf> <{}>".format(parent_uri, uri))
        return rdf

    def add_person(self, namespace: str, person_id):
        person_uri = Person.uri(namespace, person_id)
        uri = Organization.uri(namespace, self.org_id)
        rdf = []
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#associatedWith> <{}>".format(person_uri, uri))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#associationFor> <{}>".format(uri, person_uri))
        return rdf

    @staticmethod
    def uri(namespace: str, org_id: str) -> str:
        return f"{namespace}o{org_id}"


class Tool(object):
    class License:
        def __init__(self, kind: str, url: str = ''):
            self.kind = kind.strip()
            self.url = url.strip()

    class Author:
        def __init__(self, name: str, email: str = '', uri: str = ''):
            self.name = name.strip()
            self.email = email.strip()
            self.uri = uri.strip()

    def __init__(self, tool_id: Text, data: dict):
        self.tool_id: Text = tool_id.strip()
        self.name: Text = data['name'].replace('\n', ' ').strip()
        self.description: Text = data['description'].strip()
        self.url: Text = data['url'].strip()
        self.authors: List[Tool.Author] = []
        authors = data.get('authors', None) or []
        for author in authors:
            self.authors.append(Tool.Author(**author))
        license = data.get('license', dict(kind=''))
        self.license: Tool.License = Tool.License(**license)
        self.tags: List[Text] = data.get('tags', [])
        self.pmid: Text = data.get('pmid', None)
        self.approach = data.get('approach', '')
        self.functionality = data.get('functionality', '')
        self.instrumental = data.get('instrumental', '')
        self.language = data.get('language', '')
        self.type = data.get('type', '')

    def uri(self, namespace: Text) -> Text:
        encoded = self.tool_id
        encoded = encoded.replace('_', '__')
        encoded = encoded.replace('&', '_a')
        encoded = encoded.replace(':', '_c')
        encoded = encoded.replace('-', '_d')
        encoded = encoded.replace('=', '_e')
        encoded = encoded.replace('+', '_p')
        encoded = encoded.replace('?', '_q')
        encoded = encoded.replace('/', '_s')
        encoded = encoded.replace(' ', '_w')

        contains_unhandled_char = re.search('[^A-Za-z0-9._/]', encoded)
        if contains_unhandled_char:
            raise Exception(
                "Unhandled character in tool's ID: %s" % self.tool_id)

        return namespace + 't' + encoded

    def get_triples(self, namespace: Text) -> List[Text]:
        m3c = 'http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#'

        rdf = []
        uri = self.uri(namespace)

        rdf.append('<{uri}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <{m3c}Tool>'
                   .format(uri=uri, m3c=m3c))
        rdf.append('<{uri}> <http://www.w3.org/2000/01/rdf-schema#label> \"{name}\"'
                   .format(uri=uri, name=escape(self.name)))
        rdf.append('<{uri}> <{m3c}summary> \"{desc}\"'
                   .format(uri=uri, m3c=m3c, desc=escape(self.description)))

        rdf.append('<{uri}> <{m3c}homepage> \"{desc}\"'
                   .format(uri=uri, m3c=m3c, desc=escape(self.url)))

        if self.license and self.license.kind and self.license.url:
            license = self.license
            rdf.append('<{uri}> <{m3c}licenseType> \"{kind}\"'
                       .format(uri=uri, m3c=m3c, kind=escape(license.kind)))
            rdf.append('<{uri}> <{m3c}licenseUrl> \"{link}\"'
                       .format(uri=uri, m3c=m3c, link=escape(license.url)))

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
            rdf.append("<{uri}> <{m3c}tag> \"{tag}\""
                       .format(uri=uri, m3c=m3c, tag=escape(tag)))

        props: Dict[str, str] = {
            "approach": self.approach,
            "functionality": self.functionality,
            "instrumentalDataType": self.instrumental,
            "programmingLanguage": self.language,
            "softwareType": self.type,
        }

        for prop, values in props.items():
            split = values.replace(',', '\n').replace('/', '\n').split('\n')
            for value in split:
                if value in ["", "-", "?"]:
                    continue
                value = value.strip()
                rdf.append(f'<{uri}> <{m3c}{prop}> "{value}"')

        return rdf

    def match_authors(self, find_person: Callable[[str], int], namespace: str
                      ) -> List["Tool.Author"]:
        non_matched: List["Tool.Author"] = []
        log_path = os.path.join('', 'log.txt')
        log = Logger(log_path)
        for author in self.authors:
            person_id = find_person(author.name)
            if not person_id:
                non_matched.append(author)
                print(f'Unknown author "{author.name}" for tool: {self.tool_id}')
                log(f'Unknown author "{author.name}" for tool: {self.tool_id}')
            else:
                author.uri = Person.uri(namespace, str(person_id))
        return non_matched


class DateTimeValue:
    '''
    Represents a VIVO DateTimeValue.

    If `month` and `day` are specified, the precision is `yearMonthDay`.
    If `month` is specified, the precision is `yearMonth`.
    In all other cases, the precision is `year`.

    See http://vivoweb.org/ontology/core#DateTimeValue
    '''
    def __init__(self, year: int, month: int = 0, day: int = 0):
        self.precision = 'year'
        self.year = year
        self.month = 1
        self.day = 1

        if month:
            self.month = month
            self.precision = 'yearMonth'

            if day:
                self.day = day
                self.precision = 'yearMonthDay'

    def get_triples(self, datetime_value_uri: str) -> List[str]:
        uri = datetime_value_uri

        triples: List[str] = []
        triples.append(f'<{uri}> <http://vivoweb.org/ontology/core#dateTime> "{self.year:04}-{self.month:02}-{self.day:02}T00:00:00"^^<http://www.w3.org/2001/XMLSchema#dateTime>')
        triples.append(f'<{uri}> <http://vivoweb.org/ontology/core#dateTimePrecision> <http://vivoweb.org/ontology/core#{self.precision}Precision>')
        return triples


class Publication(object):
    def __init__(self, pmid: str, title: str,
                 published: Optional[DateTimeValue],
                 doi: str, citation: str):
        self.pmid = pmid
        self.title = title
        self.published = published
        self.doi = doi
        self.citation = citation
        self.authors: Set[str] = set()

    def get_triples(self, namespace):
        uri = Publication.uri(namespace, self.pmid)
        dtv_uri = f"{uri}dtv"
        rdf = []
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://purl.org/ontology/bibo/Article>".format(uri))
        rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, self.title))
        rdf.append("<{}> <http://purl.org/ontology/bibo/pmid> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, self.pmid))
        if self.doi:
            rdf.append("<{}> <http://purl.org/ontology/bibo/doi> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, self.doi))
        if self.published:
            rdf.extend(self.published.get_triples(dtv_uri))
            rdf.append("<{}> <http://vivoweb.org/ontology/core#dateTimeValue> <{}>".format(uri, dtv_uri))
        if self.citation:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#citation> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(uri, self.citation))
        for person_id in self.authors:
            pub_uri = Publication.uri(namespace, self.pmid)
            person_uri = Person.uri(namespace, person_id)
            relation_uri = f"{person_uri}r{self.pmid}"
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vivoweb.org/ontology/core#Authorship>".format(relation_uri))
            rdf.append("<{}> <http://vivoweb.org/ontology/core#relatedBy> <{}>".format(pub_uri, relation_uri))
            rdf.append("<{}> <http://vivoweb.org/ontology/core#relates> <{}>".format(relation_uri, pub_uri))
            rdf.append("<{}> <http://vivoweb.org/ontology/core#relatedBy> <{}>".format(person_uri, relation_uri))
            rdf.append("<{}> <http://vivoweb.org/ontology/core#relates> <{}>".format(relation_uri, person_uri))
        return rdf

    def add_author(self, person_id):
        self.authors.add(person_id)

    @staticmethod
    def from_pubmed(xml: str):
        fullxml = f"""<?xml version="1.0"?>
            <!DOCTYPE PubmedArticle PUBLIC "-//NLM//DTD PubMedArticle, 1st January 2019//EN" "https://dtd.nlm.nih.gov/ncbi/pubmed/out/pubmed_190101.dtd">
            {xml}
        """
        handle = io.BytesIO(fullxml.encode("utf-8"))
        article = Entrez.read(handle)
        citation = Citation(article)
        pub = make_pub(citation)
        if not pub.pmid or not pub.published:
            return None
        return pub

    @staticmethod
    def uri(namespace: str, pmid: str) -> str:
        return f"{namespace}pmid{pmid}"


def escape(text: str):
    text = text.strip()
    text = text.replace('\\', '\\\\')
    text = text.replace("\n", "\\n")
    text = text.replace('"', '\\"')
    return text


def make_pub(citation: Citation) -> Publication:
    title = citation.check_key(['MedlineCitation', 'Article', 'ArticleTitle'])
    if not title:
        title = citation.check_key(
            ['MedlineCitation', 'Article', 'VernacularTitle']
        )
    title = title.replace('"', '\\"')
    assert title

    # For more information on parsing publication dates in PubMed, see:
    #   https://www.nlm.nih.gov/bsd/licensee/elements_descriptions.html#pubdate
    published = None
    pubdate = citation.check_key(
        ['MedlineCitation', 'Article', 'Journal', 'JournalIssue', 'PubDate'])
    if pubdate:
        if 'MedlineDate' in pubdate:
            year = int(pubdate['MedlineDate'][0:4])
            assert 1900 < year and year < 3000
        else:
            year = int(pubdate['Year'])

        try:
            month_text = pubdate['Month']
            month = MONTHS.index(month_text) + 1
        except (KeyError, ValueError):
            month = 0

        try:
            day = int(pubdate['Day'])
        except KeyError:
            day = 0

        published = DateTimeValue(year, month, day)

    pmid = str(citation.check_key(['MedlineCitation', 'PMID']))
    try:
        count = 0
        proto_doi = citation.check_key(['PubmedData', 'ArticleIdList'])[count]
        while proto_doi.attributes['IdType'] != 'doi':
            count += 1
            proto_doi = citation.check_key(['PubmedData',
                                            'ArticleIdList'])[count]
        doi = str(proto_doi)
    except IndexError:
        doi = ''

    # create citation
    author_list = citation.check_key(['MedlineCitation', 'Article',
                                      'AuthorList'])
    names = []
    for author in author_list:
        if 'CollectiveName' in author:
            names.append(author['CollectiveName'])
            continue
        last_name = author['LastName']
        name = last_name
        try:
            initial = author['Initials']
            name = f"{last_name}, {initial}."
        except KeyError:
            name = last_name  # Allow surname-only authors.
        names.append(name)
    volume = citation.check_key(['MedlineCitation', 'Article', 'Journal',
                                 'JournalIssue', 'Volume'])
    issue = citation.check_key(['MedlineCitation', 'Article', 'Journal',
                                'JournalIssue', 'Issue'])
    pages = citation.check_key(['MedlineCitation', 'Article', 'Pagination',
                                'MedlinePgn'])
    journal = citation.check_key(['MedlineCitation', 'Article', 'Journal',
                                  'Title']).title()

    cite = ', '.join(names)
    if published:
        cite += f' ({published.year}). '
    cite += title
    if not cite.endswith('.'):
        cite += '. '
    else:
        cite += ' '
    if journal:
        cite += journal
        if volume or issue:
            cite += ', '
            if volume:
                cite += volume
            if issue:
                cite += '(' + issue + ')'
        if pages:
            cite += ', ' + pages
        cite += '. '
    if doi:
        cite += 'doi:' + doi
    citation_string = cite

    return Publication(pmid, title, published, doi, citation_string)


def parse_api(results: dict) -> Dict[str, Publication]:
    publications: Dict[str, Publication] = {}
    log_path = os.path.join('', 'log.txt')
    log = Logger(log_path)

    for article in results['PubmedArticle']:
        try:
            citation = Citation(article)
            pub = make_pub(citation)
            if not pub.pmid or not pub.published:
                continue
            publications[pub.pmid] = pub
        except Exception:
            citation = Citation(article)
            pmid = str(citation.check_key(['MedlineCitation', 'PMID']))
            print(f'Skipping publication {pmid}')
            log(f'Skipping publication {pmid}')
            continue

    return publications
