import typing
import unittest

from m3c import mwb
from m3c import prefill


List = typing.List


class TestPrefill(unittest.TestCase):
    def setUp(self):
        self.olddb = [
            prefill.db.add_organization,
            prefill.db.find_organizations,
            prefill.db.get_organization,
            prefill.db.add_person,
            prefill.get_person,
        ]

        prefill.db.add_organization = add_organization
        prefill.db.find_organizations = find_organizations
        prefill.db.get_organization = get_organization
        prefill.db.add_person = add_person
        prefill.get_person = get_person

    def tearDown(self):
        prefill.db.add_organization,
        prefill.db.find_organizations,
        prefill.db.get_organization,
        prefill.db.add_person,
        prefill.get_person = self.olddb
        del self.olddb

        organizations.clear()

    def test_monkeypatch(self):
        expected = "The Corporation"
        cursor = MockDatabaseConnection().cursor()
        rec = make_record(institute=expected)
        prefill.add_organizations(cursor, rec)
        self.assertEqual(len(organizations), 1)
        self.assertEqual(organizations[0], expected)

    def test_strip_names(self):
        expected = ["The Corporation", "College University"]
        cursor = MockDatabaseConnection().cursor()
        rec = make_record(institute="The Corporation; College University")
        prefill.add_organizations(cursor, rec)
        self.assertListEqual(organizations, expected)

    def test_multiname_single_institute_department_lab(self):
        expected = ["The Institute", "The Department", "The Lab"]
        cursor = MockDatabaseConnection().cursor()
        rec = make_record(institute=expected[0], department=expected[1],
                          laboratory=expected[2])
        prefill.add_organizations(cursor, rec)
        self.assertListEqual(organizations, expected)

    def test_multiname_same_number_of_institutes_departments_and_labs(self):
        cursor = MockDatabaseConnection().cursor()
        rec = make_record(institute="UF  ; FSU",
                          department="Chem; Chem",
                          laboratory="Smith;Jones")
        prefill.add_organizations(cursor, rec)
        expected = ["UF", "FSU", "Chem", "Chem", "Smith", "Jones"]
        self.assertListEqual(organizations, expected)

    def test_multiname_single_institute_and_dept_many_labs(self):
        cursor = MockDatabaseConnection().cursor()
        rec = make_record(institute="UF",
                          department="Computers",
                          laboratory="Teabeau; Clueknee")
        prefill.add_organizations(cursor, rec)
        expected = ["UF", "Computers", "Teabeau", "Clueknee"]
        self.assertListEqual(organizations, expected)

    def test_multiname_single_institute_same_number_of_depts_and_labs(self):
        cursor = MockDatabaseConnection().cursor()
        rec = make_record(institute="UF",
                          department="Chemistry; Taste and Smell",
                          laboratory="Smith;Akkbar")
        prefill.add_organizations(cursor, rec)
        expected = ["UF", "Chemistry", "Taste and Smell", "Smith", "Akkbar"]
        self.assertListEqual(organizations, expected)

    def test_multiname_fewer_depts_than_institutes_errors(self):
        # Ambiguity: which institute does the department belong to?
        rec = make_record(institute="UF;FSU",
                          department="Biology",
                          laboratory="Bobby")
        cursor = MockDatabaseConnection().cursor()
        with self.assertRaises(prefill.AmbiguityError):
            prefill.add_organizations(cursor, rec)

    def test_multiname_fewer_labs_than_departments_errors(self):
        # Ambiguity: which department does the lab belong to?
        rec = make_record(institute="UF",
                          department="Biology;Chem",
                          laboratory="Bobby")
        cursor = MockDatabaseConnection().cursor()
        with self.assertRaises(prefill.AmbiguityError):
            prefill.add_organizations(cursor, rec)

    def test_multiname_too_many_labs(self):
        # Ambiguity: which department does the last lab belong to?
        rec = make_record(institute="UF;FSU",
                          department="Biology;Chem",
                          laboratory="Bobby;Jones;Davis")
        cursor = MockDatabaseConnection().cursor()
        with self.assertRaises(prefill.AmbiguityError):
            prefill.add_organizations(cursor, rec)

    def test_multiname_too_many_departments(self):
        # Ambiguity: which institute does the last department belong to?
        rec = make_record(institute="UF;FSU",
                          department="Biology;Chem;Yo",
                          laboratory="Bobby;Jones;Davis")
        cursor = MockDatabaseConnection().cursor()
        with self.assertRaises(prefill.AmbiguityError):
            prefill.add_organizations(cursor, rec)

    def test_multiname_single_person(self):
        rec = make_record(last_name="Bond", first_name="James")
        cursor = MockDatabaseConnection().cursor()
        actual = prefill.add_people(cursor, rec)
        self.assertEqual(len(actual), 1)
        self.assertEqual(people[0], "James Bond")

    def test_multiname_too_few_surnames(self):
        rec = make_record(last_name="Bond", first_name="James;Michael")
        cursor = MockDatabaseConnection().cursor()
        with self.assertRaises(prefill.AmbiguousNamesError):
            prefill.add_people(cursor, rec)

    def test_multiname_too_many_emails(self):
        rec = make_record(last_name="Bond", first_name="James",
                          email="007@secret.gov.uk and foo@example.com")
        cursor = MockDatabaseConnection().cursor()
        prefill.add_people(cursor, rec)
        self.assertEqual(emails[0], "")


organizations: List[str] = []
people: List[str] = []
emails: List[str] = []


def add_organization(cursor, type, name, parent_id=None):
    organizations.append(name)
    return len(organizations)


def add_person(cursor, first_name, last_name, email, phone):
    people.append(f"{first_name} {last_name}")
    emails.append(email)
    return len(people)


def find_organizations(cursor):
    return []


def get_organization(cursor, type, name, parent_id=None):
    if not parent_id:
        parent_id = 0
    try:
        return organizations.index((name, parent_id))+1
    except ValueError:
        return 0


def get_person(cursor, first_name, last_name, exclude_withheld=True):
    return []


class MockDatabaseConnection:
    def cursor(self):
        return self

    def __enter__(self):
        pass

    def __exit__(self, a, b, c):
        pass


def make_record(psid="PR123", pstype=mwb.PROJECT, first_name="", last_name="",
                institute="", department="", laboratory="", email="", phone=""
                ) -> mwb.NameRecord:
    return mwb.NameRecord(psid, pstype, first_name, last_name, institute,
                          department, laboratory, email, phone)


if __name__ == "__main__":
    unittest.main()
