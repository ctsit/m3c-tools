class Aide(object):
    def __init__(self, endpoint: str, email: str, password: str, prefix: str):
        # Endpoint must be update endpoint, not query
        self.endpoint = endpoint
        self.email = email
        self.password = password
        self.namespace = prefix
        self.used_n_numbers = []
