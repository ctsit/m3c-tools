from typing import List, Tuple
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
        projects.append((expected, "", ""))

        conn = MockDatabaseConnection()
        metab_prefill.add_organizations(conn, conn, [])

        self.assertEqual(len(organizations), 1)
        self.assertEqual(organizations[0], expected)

    def test_strip_names(self):
        expected = ["The Corporation", "College University"]
        projects.append(("The Corporation; College University", "", ""))

        conn = MockDatabaseConnection()
        metab_prefill.add_organizations(conn, conn, [])

        self.assertListEqual(organizations, expected)

    def test_multiname_single_institute_department_lab(self):
        expected = ["The Institute", "The Department", "The Lab"]
        projects.append(tuple(expected))
        conn = MockDatabaseConnection()
        metab_prefill.add_organizations(conn, conn, [])
        self.assertListEqual(organizations, expected)

    def test_multiname_same_number_of_institutes_departments_and_labs(self):
        projects.append(("UF  ; FSU", "Chem; Chem", "Smith;Jones"))
        conn = MockDatabaseConnection()
        metab_prefill.add_organizations(conn, conn, [])
        expected = ["UF", "FSU", "Chem", "Chem", "Smith", "Jones"]
        self.assertListEqual(organizations, expected)

    def test_multiname_single_institute_and_dept_many_labs(self):
        projects.append(("UF", "Computers", "Teabeau; Clueknee"))
        conn = MockDatabaseConnection()
        metab_prefill.add_organizations(conn, conn, [])
        expected = ["UF", "Computers", "Teabeau", "Clueknee"]
        self.assertListEqual(organizations, expected)

    def test_multiname_single_institute_same_number_of_depts_and_labs(self):
        projects.append(("UF", "Chemistry; Taste and Smell", "Smith;Akkbar"))
        conn = MockDatabaseConnection()
        metab_prefill.add_organizations(conn, conn, [])
        expected = ["UF", "Chemistry", "Taste and Smell", "Smith", "Akkbar"]
        self.assertListEqual(organizations, expected)

    def test_multiname_fewer_depts_than_institutes_errors(self):
        # Ambiguity: which institute does the department belong to?
        projects.append(("UF;FSU", "Biology", "Bobby"))
        conn = MockDatabaseConnection()
        with self.assertRaises(metab_prefill.AmbiguityError):
            metab_prefill.add_organizations(conn, conn, [])

    def test_multiname_fewer_labs_than_departments_errors(self):
        # Ambiguity: which department does the lab belong to?
        projects.append(("UF", "Biology;Chem", "Bobby"))
        conn = MockDatabaseConnection()
        with self.assertRaises(metab_prefill.AmbiguityError):
            metab_prefill.add_organizations(conn, conn, [])

    def test_multiname_too_many_labs(self):
        # Ambiguity: which department does the last lab belong to?
        projects.append(("UF;FSU", "Biology;Chem", "Bobby;Jones;Davis"))
        conn = MockDatabaseConnection()
        with self.assertRaises(metab_prefill.AmbiguityError):
            metab_prefill.add_organizations(conn, conn, [])

    def test_multiname_too_many_departments(self):
        # Ambiguity: which institute does the last department belong to?
        projects.append(("UF;FSU", "Biology;Chem;Yo", "Bobby;Jones;Davis"))
        conn = MockDatabaseConnection()
        with self.assertRaises(metab_prefill.AmbiguityError):
            metab_prefill.add_organizations(conn, conn, [])


projects: List[Tuple[str, str, str]] = []
organizations: List[Tuple[str, int]] = []


def add_organization(cursor, type, name, parent_id=None):
    organizations.append(name)
    return len(organizations)


def find_organizations(cursor):
    for (institute, department, lab) in projects:
        yield [institute, department, lab, "PROJECT_ID"]


def get_organization(cursor, type, name, parent_id=None):
    if not parent_id:
        parent_id = 0
    try:
        return organizations.index((name, parent_id))+1
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
