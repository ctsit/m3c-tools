#!/usr/bin/python

import random
import requests

class Aide(object):
    def __init__(self, query_endpoint, email, password, namespace):
        self.endpoint = query_endpoint
        self.email = email
        self.password = password
        self.namespace = namespace
        self.used_n_numbers = []

    def do_query(self, query, silent=False):
        if not silent:
            print("Query:\n" + query)
        payload = {
            'email': self.email,
            'password': self.password,
            'query': query
        }
        headers = {'Accept': 'application/sparql-results+json'}
        try:
            response = requests.get(self.endpoint, params=payload, headers=headers, verify=False)
            if response.status_code == 400:
                print(query)
                exit("400 Error: check query")
            elif response.status_code == 403:
                print(query)
                exit("403 Error: check credentials")
            elif response.status_code == 406:
                print(query)
                exit("406 Error: check accept header")
            if response.status_code == 500:
                print(query)
                exit("500 Error: check server")
        except requests.exceptions.ConnectionError as e:
            print("Error with this query.")
            if silent:
                print("Query:\n" + query)
            print(e)
            response = None
        return response

    def make_n_number(self):
        bad_n = True
        while bad_n:
            # create an n number
            uri = self.namespace + "n" + str(random.randint(1,999999999))
            # check if number is taken
            if uri in self.used_n_numbers:
                bad_n = True
            else:
                bad_n = self.check_n_number(uri)
        self.used_n_numbers.append(uri)
        return uri

    def check_n_number(self, uri):
        query = """\
            ASK {{
                <{}> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#Thing>
            }} 
        """.format(uri)

        response = self.do_query(query, True)
        res = response.json()
        exists = res['boolean']

        return exists