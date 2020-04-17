import typing
import yaml

Dict = typing.Dict
Optional = typing.Optional


class Config:
    def __init__(self, endpoint: str, email: str, password: str, prefix: str,
                 data: Optional[Dict[str, str]] = {}):
        # Endpoint must be the Update API's endpoint, not the Query API's.
        self.endpoint = endpoint
        self.email = email
        self.password = password
        self.namespace = prefix
        self.used_n_numbers = []
        self._data = data or {}

    def get(self, prop, default=None):
        return self._data.get(prop, default)


def load(yamlfile: str):
    try:
        with open(yamlfile, "r") as fp:
            data = yaml.load(fp.read(), Loader=yaml.FullLoader)
        config = Config(data.get("update_endpoint"),
                        data.get("vivo_email"),
                        data.get("vivo_password"),
                        data.get("namespace"),
                        data)
        return config
    except Exception:
        print("Error: Check config file")
        raise
