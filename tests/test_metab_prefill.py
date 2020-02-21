import unittest

import metab_prefill


class TestMetabPrefill(unittest.TestCase):
    def setUp(self):
        self.olddb = [
            metab_prefill.db.add_organization,
            metab_prefill.db.find_organizations,
            metab_prefill.db.get_organization,
        ]

        metab_prefill.db.add_organization = add_organization
        metab_prefill.db.find_organizations = find_organizations
        metab_prefill.db.get_organization = get_organization

    def tearDown(self):
        metab_prefill.db.add_organization,
        metab_prefill.db.find_organizations,
        metab_prefill.db.get_organization = self.olddb
        del self.olddb

        organizations.clear()
        projects.clear()

    def test_monkeypatch(self):
        expected = "The Corporation"
        projects.append(expected)

        conn = MockDatabaseConnection()
        metab_prefill.add_organizations(conn, conn, [])

        self.assertEqual(len(organizations), 1)
        self.assertEqual(organizations[0], expected)

    def test_strip_names(self):
        expected = ["The Corporation", "College University"]
        projects.append("The Corporation; College University")

        conn = MockDatabaseConnection()
        metab_prefill.add_organizations(conn, conn, [])

        self.assertListEqual(organizations, expected)


projects = []
organizations = []


def add_organization(cursor, type, name, parent_id=None):
    organizations.append(name)
    return len(organizations)


def find_organizations(cursor):
    for institute in projects:
        yield [institute, "", "", "PROJECT_ID"]


def get_organization(cursor, type, name, parent_id=None):
    try:
        return organizations.index(name)+1
    except ValueError:
        return 0


class MockDatabaseConnection:
    def cursor(self):
        return self

    def __enter__(self):
        pass

    def __exit__(self, a, b, c):
        pass


if __name__ == "__main__":
    unittest.main()
