"""
Metab Pub Ingest
Usage:
    metab_pub_ingest.py (-h | --help)
    metab_pub_ingest.py [-id <id_number>] <path_to_config>

Options:
    -h --help       Show this message and exit
    -id             Run ingest for single person_id

Example:
    $ python metab_pub_ingest.py config.yaml
"""

from datetime import datetime
import os
import sys
import yaml

import psycopg2

from aide import Aide
from metab_classes import Person
from metab_classes import Publication


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


def get_config(config_path):
    try:
        with open(config_path, 'r') as config_file:
            config = yaml.load(config_file.read(), Loader=yaml.FullLoader)
    except Exception as e:
        print("Error: Check config file")
        sys.exit(e)
    return config


def connect(host, db, user, pg_password, port):
    conn = psycopg2.connect(host=host, dbname=db, user=user,
                            password=pg_password, port=port)
    cur = conn.cursor()
    return cur


def get_people(cur, person_id=None):
    people = {}
    if person_id:
        cur.execute("""\
                SELECT id, first_name, last_name
                FROM people
                JOIN names
                ON id=person_id
                WHERE id=%s""", (int(person_id),))
    else:
        cur.execute("""\
                SELECT id, first_name, last_name
                FROM people
                JOIN names
                ON id=person_id""")
    for row in cur:
        person = Person()
        person.person_id = row[0]
        person.first_name = row[1]
        person.last_name = row[2]
        people[row[0]] = person
    return people


def get_supplementals(cur, person_id=None):
    extras = {}
    exceptions = {}
    if person_id:
        cur.execute("""\
            SELECT pmid, person_id, include
            FROM publications
            WHERE person_id=%s""", (person_id,))
    else:
        cur.execute("""\
            SELECT pmid, person_id, include
            FROM publications""")
    for row in cur:
        pmid = row[0]
        person_id = int(row[1])
        include = row[2]
        if include:
            if person_id in extras.keys():
                extras[person_id].append(pmid)
            else:
                extras[person_id] = [pmid]
        else:
            if person_id in exceptions.keys():
                exceptions[person_id].append(pmid)
            else:
                exceptions[person_id] = [pmid]

    return extras, exceptions


def get_ids(aide, person):
    query = person.last_name + ', ' + person.first_name + ' [Full Author Name]'
    id_list = aide.get_id_list(query)
    return id_list


def parse_api(results, namespace):
    publications = {}
    for citing in results['PubmedArticle']:
        pub = Publication()
        citation = Citation(citing)

        fill_pub(pub, citation)

        if pub.pmid:
            pub.uri = namespace + pub.pmid
            publications[pub.pmid] = pub
    return publications


def fill_pub(pub, citation):
    pub.title = (citation.check_key(['MedlineCitation', 'Article', 'ArticleTitle'])).replace('"', '\\"')
    pub.publication_year = (citation.check_key(['MedlineCitation', 'Article', 'Journal',
                            'JournalIssue', 'PubDate', 'Year']))
    pub.pmid = str(citation.check_key(['MedlineCitation', 'PMID']))
    try:
        count = 0
        proto_doi = citation.check_key(['PubmedData', 'ArticleIdList'])[count]
        while proto_doi.attributes['IdType'] != 'doi':
            count += 1
            proto_doi = citation.check_key(['PubmedData',
                                            'ArticleIdList'])[count]
        pub.doi = str(proto_doi)
    except IndexError:
        pub.doi = ''
    # create citation
    author_list = citation.check_key(['MedlineCitation', 'Article', 'AuthorList'])
    names = []
    for author in author_list:
        last_name = author['LastName']
        initial = author['Initials']
        name = last_name + ", " + initial + "."
        names.append(name)
    volume = citation.check_key(['MedlineCitation', 'Article', 'Journal', 'JournalIssue',
                                'Volume'])
    issue = citation.check_key(['MedlineCitation', 'Article', 'Journal', 'JournalIssue',
                                'Issue'])
    pages = citation.check_key(['MedlineCitation', 'Article', 'Pagination', 'MedlinePgn'])
    journal = citation.check_key(['MedlineCitation', 'Article', 'Journal', 'Title']).title()

    cite = ', '.join(names)
    if pub.publication_year:
        cite += ' (' + pub.publication_year + '). '
    cite += pub.title
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
    if pub.doi:
        cite += 'doi:' + pub.doi
    pub.citation = cite
    return


def write_triples(aide, person, pubs):
    rdf = []
    for pub in pubs.values():
        rdf.extend(pub.add_person(aide.namespace, person.person_id))
    return rdf


def print_to_file(triples, file):
    triples = [t + " ." for t in triples]
    with open(file, 'a+') as rdf:
        rdf.write("\n".join(triples))


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)

    if sys.argv[1] in ["-h", "--help"]:
        print(__doc__)
        sys.exit()

    if sys.argv[1] == "-id":
        person_id = sys.argv[2]
        config_path = sys.argv[3]
    else:
        person_id = None
        config_path = sys.argv[1]

    timestamp = datetime.now()
    path = 'data_out/' + timestamp.strftime("%Y") + '/' + \
        timestamp.strftime("%m") + '/' + timestamp.strftime("%Y_%m_%d")
    try:
        os.makedirs(path)
    except FileExistsError:
        pass
    if sys.argv[1] == "-id":
        pub_file = os.path.join(path, person_id + "_pubs.rdf")
    else:
        pub_file = os.path.join(path, 'pubs.rdf')

    config = get_config(config_path)

    aide = Aide(config.get('update_endpoint'),
                config.get('vivo_email'),
                config.get('vivo_password'),
                config.get('namespace'))

    cur = connect(config.get('sup_host'), config.get('sup_database'),
                  config.get('sup_username'), config.get('sup_password'),
                  config.get('sup_port'))

    people = get_people(cur, person_id)
    triples = []
    pub_collective = {}

    extras, exceptions = get_supplementals(cur, person_id)
    for person in people.values():
        pmids = get_ids(aide, person)
        if person.person_id in extras.keys():
            for pub in extras[person.person_id]:
                if pub not in pmids:
                    pmids.append(pub)
        if person.person_id in exceptions.keys():
            for pub in exceptions[person.person_id]:
                if pub in pmids:
                    pmids.remove(pub)
        if pmids:
            results = aide.get_details(pmids)
            pubs = parse_api(results, aide.namespace)
            pub_collective.update(pubs)
            triples.extend(write_triples(aide, person, pubs))

    pub_count = 0
    for pub in pub_collective.values():
        triples.extend(pub.get_triples())
        pub_count += 1
    print("There are " + str(pub_count) + " publications.")

    print_to_file(triples, pub_file)


if __name__ == "__main__":
    main()