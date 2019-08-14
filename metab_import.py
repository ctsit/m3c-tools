#!/usr/bin/python

"""
Metab Importer
Usage:
    import.py (-h | --help)
    import.py <path_to_config>

Options:
    -h --help       Show this message and exit

Instructions:
    Run the importer where you have access to the postgres metabolomics database.
"""

from datetime import datetime
import os
import sys
import yaml

import psycopg2

from aide import Aide
from metab_classes import Dataset
from metab_classes import Person
from metab_classes import Project
from metab_classes import Study


def get_config(config_path):
    try:
        with open(config_path, 'r') as config_file:
            config = yaml.load(config_file.read(), Loader=yaml.FullLoader)
    except Exception as e:
        print("Error: Check config file")
        exit(e)
    return config


def connect(host, db, user, pg_password, port):
    conn = psycopg2.connect(host=host,dbname=db,user=user,
                            password=pg_password,port=port)
    cur = conn.cursor()
    return cur


def get_projects(cur):
    print("Gathering Workbench Projects")
    projects = {}
    cur.execute("""\
                SELECT project_id, project_title, project_type, project_summary,
                  doi, funding_source, publications, last_name, first_name,
                  email, institute, department, laboratory
                FROM project""")
    for row in cur:
        project = Project()
        project.project_id = row[0]
        project.project_title = row[1].replace('"', '\\"')
        project.project_type = row[2]
        project.summary = row[3].replace('"', '\\"').replace('\n', ' ')
        project.doi = row[4]
        project.funding_source = row[5]
        project.publication_string = row[6]
        project.last_name = row[7]
        project.first_name = row[8]
        project.email = row[9]
        project.institute = row[10]
        project.department = row[11]
        project.laboratory = row[12]
        projects[project.project_id] = project
    return projects


def make_projects(aide, projects):
    print("Making Workbench Projects")
    triples = []
    project_count = 0

    for project in projects.values():
        if not project.check_existence(aide):
            project_uri = aide.make_n_number()
            project.uri = project_uri
            rdf = project.get_triples()
            triples.extend(rdf)
            project_count += 1

    print("There will be " + str(project_count) + " new projects.")
    return triples


def get_studies(cur):
    print("Gathering Workbench Studies")
    studies = {}
    cur.exectue("""\
        SELECT study_id, study_title, study_type, study_summary, publications,
            project_id, last_name, first_name, email, institute, department,
            laboratory
        FROM study""")
    for row in cur:
        study = Study()
        study.study_id = row[0]
        study.study_title = row[1].replace('"', '\\"')
        study.study_type = row[2]
        study.summary = row[3].replace('"', '\\"').replace('\n', ' ')
        study.publication_string = row[4]
        study.last_name = row[6]
        study.first_name = row[7]
        study.email = row[8]
        study.institute = row[9]
        study.department = row[10]
        study.laboratory = row[11]
        study.project_id = row[5]
        studies[study.study_id] = study
    return studies


def make_studies(aide, studies, projects):
    print("Making Workbench Studies")
    triples = []
    study_count = 0
    for study in studies.values():
        if not study.check_existence(aide):
            study_uri = aide.make_n_number()
            study.uri = study_uri
            try:
                project_uri = projects[study.project_id].uri
            except KeyError:
                project_uri = None
            rdf = study.get_triples(project_uri)
            triples.extend(rdf)
            study_count += 1

    print("There will be " + str(study_count) + " new studies.")
    return triples


def make_people(aide, people, works):
    '''
    Make profiles for people who are related to projects or studies
    '''
    print("Making and linking related people.")
    triples = []
    people_count = 0

    for work in works.values():
        person = Person()
        person.first_name = work.first_name
        person.last_name = work.last_name
        person.email = work.email
        person.institute = work.institute
        person.department = work.department
        person.laboratory = work.laboratory
        person.make_display_name()

        if person.display_name not in people.keys():
            people[person.display_name] = person
        else:
            person = people[person.display_name]

        exists = person.check_existence(aide)
        if not exists:
            person_uri = aide.make_n_number()
            person.uri = person_uri
            rdf = person.get_triples(aide)
            triples.extend(rdf)
            linking_rdf = work.add_person(person.uri)
            triples.extend(linking_rdf)
        else: 
            if not work.check_link(aide, person.uri):
                linking_rdf = work.add_person(person.uri)
                triples.extend(linking_rdf)
    print("There will be " + str(people_count) + " new people in this set.")
    return triples


def get_datasets(cur):
    print("Gathering Workbench Datasets")
    datasets = {}
    cur.execute("""\
        SELECT mb_sample_id, study_id
        FROM metadata""")
    for row in cur:
        dataset = Dataset()
        dataset.mb_sample_id = row[0]
        dataset.study_id = row[1]
        datasets[dataset.mb_sample_id] = dataset
    return datasets


def make_datasets(aide, datasets, studies):
    print("Making Workbench Datasets")
    triples = []
    dataset_count = 0
    for dataset in datasets.values():
        if not dataset.check_existence(aide):
            dataset_uri = aide.make_n_number()
            dataset.uri = dataset_uri
            try:
                study_uri = studies[dataset.study_id].uri
            except KeyError:
                study_uri = None
            rdf = dataset.get_triples(study_uri)
            triples.extend(rdf)
            dataset_count += 1

    print("There will be " + str(dataset_count) + " new datasets.")
    return triples


def print_to_file(triples, file):
    with open(file, 'a+') as rdf:
        rdf.write(" . \n".join(triples))
        rdf.write(" . \n")


def main(config_path):
    timestamp = datetime.now()
    path = 'data_out/' + timestamp.strftime("%Y") + '/' +\
            timestamp.strftime("%m") + '/' + timestamp.strftime("%Y_%m_%d")
    try:
        os.makedirs(path)
    except FileExistsError:
        pass
    project_file = os.path.join(path, 'projects.rdf')
    study_file = os.path.join(path, 'studies.rdf')
    dataset_file = os.path.join(path, 'datasets.rdf')
    people_file = os.path.join(path, 'people.rdf')

    config = get_config(config_path)
    aide = Aide(config.get('query_endpoint'),
                config.get('vivo_email'),
                config.get('vivo_password'),
                config.get('namespace'))
    cur = connect(config.get('host'), config.get('database'),
                  config.get('pg_username'), config.get('pg_password'),
                  config.get('port'))
    people = {}

    # get all projects
    projects = get_projects(cur)
    project_triples = make_projects(aide, projects)
    if project_triples:
        print_to_file(project_triples, project_file)
    project_people_triples = make_people(aide, people, projects)
    if project_people_triples:
        print_to_file(project_people_triples, people_file)
    # get all studies
    studies = get_studies(cur)
    study_triples = make_studies(aide, studies, projects)
    if study_triples:
        print_to_file(study_triples, study_file)
    study_people_triples = make_people(aide, people, studies)
    if study_people_triples:
        print_to_file(study_people_triples, people_file)
    # get all datasets
    datasets = get_datasets(cur)
    dataset_triples = make_datasets(aide, datasets, studies)
    if dataset_triples:
        print_to_file(dataset_triples, dataset_file)

if __name__ == '__main__':
    main(sys.argv[1])