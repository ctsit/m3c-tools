class Project(object):
    def __init__(self):
        self.uri = None
        self.project_id = None
        self.project_type = None
        self.summary = None
        self.doi = None
        self.funding_source = None
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
