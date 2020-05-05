from typing import Dict, Optional

import traceback
import yaml


class Config:
    def __init__(self, endpoint: str, email: str, password: str, prefix: str,
                 data: Optional[Dict[str, str]] = None):
        # Endpoint must be the Update API's endpoint, not the Query API's.
        self.endpoint = endpoint
        self.email = email
        self.password = password
        self.namespace = prefix
        self._data = data or {}

    def get(self, prop, default=None):
        return self._data.get(prop, default)


def load(yamlfile: str) -> Config:
    try:
        with open(yamlfile, "r") as fp:
            data = yaml.safe_load(fp)
        config = Config(data.get("update_endpoint"),
                        data.get("vivo_email"),
                        data.get("vivo_password"),
                        data.get("namespace"),
                        data)
        return config
    except Exception:
        traceback.print_exc()
        print("Error: Check config file")
        raise
