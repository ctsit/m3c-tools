"""
Metab Importer
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
    people = []
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
        people.append(person)
    return people


def query_pubmed(aide, person):
    query = person.last_name + ', ' + person.first_name + ' [Full Author Name]'
    id_list = aide.get_id_list(query)
    results = aide.get_details(id_list)
    return results


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
    for person in people:
        results = query_pubmed(aide, person)
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
