from Bio import Entrez


class Aide(object):
    def __init__(self, endpoint, email, password, namespace):
        # VIVO related
        # Endpoint must be update endpoint, not query
        self.endpoint = endpoint
        self.email = email
        self.password = password
        self.namespace = namespace
        self.used_n_numbers = []

        # Pubmed related
        self.retstart = 0
        self.count_up = 0

    def get_id_list(self, term):
        Entrez.email = self.email
        res = Entrez.esearch(term=term,
                             db='pubmed',
                             retmax=100000,
                             retstart=self.retstart)
        result = Entrez.read(res)
        id_list = result['IdList']
        total = result['Count']

        self.count_up += 100000
        # if # of results exceeds 100,000, you will need to run the query again
        if self.count_up < int(total):
            self.retstart += 100000
            id_list += self.get_id_list(term)

        return id_list

    def get_details(self, id_list):
        ids = ','.join(id_list)
        Entrez.email = self.email
        handle = Entrez.efetch(db='pubmed',
                               retmode='xml',
                               id=ids)
        results = Entrez.read(handle)
        return results
