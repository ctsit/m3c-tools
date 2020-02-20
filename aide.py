from Bio import Entrez
import xml.etree.ElementTree as ET


class Aide(object):
    def __init__(self, endpoint: str, email: str, password: str,
                 namespace: str, pubmed_email: str, pubmed_api_key):
        # VIVO related
        # Endpoint must be update endpoint, not query
        self.endpoint = endpoint
        self.email = email
        self.password = password
        self.namespace = namespace
        self.used_n_numbers = []

        # Pubmed related
        Entrez.email = pubmed_email
        Entrez.api_key = pubmed_api_key

    def get_id_list(self, term: str, retstart: int = 0, count_up: int = 0):
        RETMAX = 100000
        handle = Entrez.esearch(term=term,
                                db="pubmed",
                                retmax=RETMAX,
                                retstart=retstart)
        result = Entrez.read(handle)
        handle.close()
        id_list = result["IdList"]
        total = int(result["Count"])

        count_up += RETMAX
        # if # of results exceeds 100,000, you will need to run the query again
        if count_up < total:
            retstart += RETMAX
            id_list += self.get_id_list(term, retstart, count_up)

        return id_list

    def get_details(self, id_list):
        ids = ",".join(id_list)
        handle = Entrez.efetch(db="pubmed", retmode="xml", id=ids)
        result = ET.parse(handle)
        handle.close()
        return result
