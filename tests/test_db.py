import sqlite3
import unittest

from m3c import db


class TestDb(unittest.TestCase):
    def setUp(self):
        self.conn = sqlite3.connect(":memory:")
        self.conn.executescript("""
            CREATE TABLE names (
                person_id  INT,
                first_name TEXT,
                last_name  TEXT,
                withheld   BOOLEAN
            );
            INSERT INTO names VALUES (7, "James", "Bond", 0);
        """)

    def tearDown(self):
        self.conn.close()

    def test_case_insensitive_name_matching(self):
        cursor = self.conn.cursor()
        actual = list(db.get_person(cursor, "james", "bond", False))
        expected = [7]
        self.assertListEqual(expected, actual)


if __name__ == "__main__":
    unittest.main()
