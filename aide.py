#!/usr/bin/python

import random  # maybe don't need-- increment number + related id?
import requests

class Aide(object):
    def __init__(self, endpoint, email, password, namespace):
        # Endpoint must be update endpoint, not query
        self.endpoint = endpoint
        self.email = email
        self.password = password
        self.namespace = namespace
        self.used_n_numbers = []

    def do_update(self, query):
        payload = {
            'email': self.email,
            'password': self.password,
            'query': query
        }
        response = requests.post(endpoint, params=payload, verify=False)
        if response.status_code != 200:
            print('Error with following query:')
            print(query)
            print(response.status_code)
            return False
        else:
            return True

    def do_delete(self):
        query='''\
            DELETE {?s ?p ?o}
            WHERE {
                ?s ?p ?o
            }
            '''
        self.do_update(query)