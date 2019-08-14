#!/usr/bin/python


class Project(object):
    def __init__(self):
        self.uri = None
        self.project_id = None
        self.project_title = None
        self.project_type = None
        self.summary = None
        self.doi = None
        self.funding_source = None
        self.publication_string = None
        self.last_name = None
        self.first_name = None
        self.email = None
        self.institute = None
        self.department = None
        self.laboratory = None

    def check_existence(self, aide):
        query = """\
            SELECT ?uri
            WHERE {{
                ?uri <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Project> .
                ?uri <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#projectId> "{}"^^<http://www.w3.org/2001/XMLSchema#string> .
            }}
        """.format(self.project_id)

        if not self.uri:
            response = aide.do_query(query, True)
            res = response.json()
            uri = None
            try:
                uri = res['results']['bindings'][0]['uri']['value']
                self.uri = uri
                return True
            except IndexError:
                return False
        else:
            return True

    def get_triples(self):
        rdf = []
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Project>".format(self.uri))
        rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.project_title))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#projectId> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.project_id))
        if self.project_type:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#projectType> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self. project_type))
        if self.summary:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#summary> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.summary))
        if self.doi:
            rdf.append("<{}> <http://purl.org/ontology/bibo/doi> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.doi))
        return rdf

    def check_link(self, aide, person_uri):
        query = """\
            ASK {{
                <{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#hasPI> <{}>
            }}
        """.format(self.uri, person_uri)

        response = aide.do_query(query, True)
        res = response.json()
        linked = res['boolean']

        return linked

    def add_person(self, person_uri):
        rdf = []
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#hasPI> <{}>".format(self.uri, person_uri))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#isPIFor> <{}>".format(person_uri, self.uri))
        return rdf
        

class Study(object):
    def __init__(self):
        self.uri = None
        self.study_id = None
        self.study_title = None
        self.study_type = None
        self.summary = None
        self.publication_string = None
        self.last_name = None
        self.first_name = None
        self.email = None
        self.institute = None
        self.department = None
        self.laboratory = None
        self.project_id = None

    def check_existence(self, aide):
        query = """\
            SELECT ?uri
            WHERE {{
                ?uri <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Study> .
                ?uri <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#studyId> "{}"^^<http://www.w3.org/2001/XMLSchema#string> .
            }}
        """.format(self.study_id)

        if not self.uri:
            response = aide.do_query(query, True)
            res = response.json()
            uri = None
            try:
                uri = res['results']['bindings'][0]['uri']['value']
                self.uri = uri
                return True
            except IndexError:
                return False
        else:
            return True

    def get_triples(self, project_uri=None):
        rdf = []
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Study>".format(self.uri))
        rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.study_title))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#studyId> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.study_id))
        if self.study_type:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#studyType> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self. study_type))
        if self.summary:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#summary> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.summary))
        if project_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#inCollection> <{}>".format(self.uri, project_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#collectionFor> <{}>".format(project_uri, self.uri))
        return rdf

    def check_link(self, aide, person_uri):
        query = """\
            ASK {{
                <{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#runBy> <{}>
            }}
        """.format(self.uri, person_uri)

        response = aide.do_query(query, True)
        res = response.json()
        linked = res['boolean']

        return linked

    def add_person(self, person_uri):
        rdf = []
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#runBy> <{}>".format(self.uri, person_uri))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#runnerOf> <{}>".format(person_uri, self.uri))
        return rdf


class Dataset(object):
    def __init__(self):
        self.uri = None
        self.mb_sample_id = None
        self.sample_species = None
        self.study_id = None

    def check_existence(self, aide):
        query = """\
            SELECT ?uri
            WHERE {{
                ?uri <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Dataset> .
                ?uri <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#sampleId> "{}"^^<http://www.w3.org/2001/XMLSchema#string> .
            }}
        """.format(self.mb_sample_id)

        if not self.uri:
            response = aide.do_query(query, True)
            res = response.json()
            uri = None
            try:
                uri = res['results']['bindings'][0]['uri']['value']
                self.uri = uri
                return True
            except IndexError:
                return False
        else:
            return True

    def get_triples(self, study_uri=None):
        rdf = []
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#Dataset>".format(self.uri))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#sampleId> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.mb_sample_id))
        rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#subjectSpecies> \"{}\"".format(self.uri, self.sample_species))
        if study_uri:
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#dataFor> <{}>".format(self.uri, study_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#developedFrom> <{}>".format(study_uri, self.uri))
        return rdf


class Person(object):
    def __init__(self):
        self.uri = None
        self.first_name = None
        self.last_name = None
        self.display_name = None
        self.email = None
        self.institute = None
        self.department = None
        self.laboratory = None

    def make_display_name(self):
        display_name = self.first_name + ' ' + self.last_name
        self.display_name = display_name

    def check_existence(self, aide):
        query = """\
            SELECT ?uri
            WHERE {{
                ?uri <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://xmlns.com/foaf/0.1/Person> .
                ?uri <http://www.w3.org/2000/01/rdf-schema#label> "{}"^^<http://www.w3.org/2001/XMLSchema#string> .
            }}
        """.format(self.display_name)

        if not self.uri:
            response = aide.do_query(query, True)
            res = response.json()
            uri = None
            try:
                uri = res['results']['bindings'][0]['uri']['value']
                self.uri = uri
                return True
            except IndexError:
                return False
        else:
            return True

    def get_triples(self, aide):
        rdf = []
        vcard_uri = aide.make_n_number()
        name_uri = aide.make_n_number()
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://xmlns.com/foaf/0.1/Person>".format(self.uri))
        rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(self.uri, self.display_name))
        rdf.append("<{}> <http://purl.obolibrary.org/obo/ARG_2000028> <{}>".format(self.uri, vcard_uri))
        rdf.append("<{}> <http://purl.obolibrary.org/obo/ARG_2000029> <{}>".format(vcard_uri, self.uri))
        rdf.append("<{}> <http://www.w3.org/2006/vcard/ns#hasName> <{}>".format(vcard_uri, name_uri))
        rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2006/vcard/ns#Name>".format(name_uri))
        rdf.append("<{}> <http://www.w3.org/2006/vcard/ns#familyName> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(vcard_uri, self.last_name))
        rdf.append("<{}> <http://www.w3.org/2006/vcard/ns#givenName> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(vcard_uri, self.first_name))
        if self.email:
            email_uri = aide.make_n_number()
            rdf.append("<{}> <http://www.w3.org/2006/vcard/ns#hasEmail>  <{}>".format(vcard_uri, email_uri))
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2006/vcard/ns#Email>".format(email_uri))
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2006/vcard/ns#Work>".format(email_uri))
            rdf.append("<{}> <http://www.w3.org/2006/vcard/ns#email> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(email_uri, self.email))
        if self.institute:
            inst_uri = aide.make_n_number()
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#associatedWith> <{}>".format(self.uri, inst_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#associationFor> <{}>".format(inst_uri, self.uri))
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vivoweb.org/ontology/core#Institute>".format(inst_uri))
            rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(inst_uri, self.institute))
        if self.department:
            dept_uri = aide.make_n_number()
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#associatedWith> <{}>".format(self.uri, dept_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#associationFor> <{}>".format(dept_uri, self.uri))
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vivoweb.org/ontology/core#Department>".format(dept_uri))
            rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(dept_uri, self.department))
        if self.laboratory:
            lab_uri = aide.make_n_number()
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#associatedWith> <{}>".format(self.uri, lab_uri))
            rdf.append("<{}> <http://www.metabolomics.info/ontologies/2019/metabolomics-consortium#associationFor> <{}>".format(lab_uri, self.uri))
            rdf.append("<{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://vivoweb.org/ontology/core#Laboratory>".format(lab_uri))
            rdf.append("<{}> <http://www.w3.org/2000/01/rdf-schema#label> \"{}\"^^<http://www.w3.org/2001/XMLSchema#string>".format(lab_uri, self.laboratory))
        return rdf 