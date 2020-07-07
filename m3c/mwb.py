"""
Metabolomics Workbench library
"""

from typing import Iterable, Optional

import psycopg2
import psycopg2.extensions


Connection = psycopg2.Connection
Cursor = psycopg2.Cursor


DEPARTMENT = "department"
INSTITUTE = "institute"
LABORATORY = "laboratory"
PROJECT = "project"
STUDY = "study"


class NameRecord:
    def __init__(self, psid: str, pstype: str, first_name: str, last_name: str,
                 institute: str, department: str, laboratory: str, email: str,
                 phone: str):
        assert psid[0:len("PR")] in ["PR", "ST"]
        assert pstype in [PROJECT, STUDY]
        self.psid = psid
        self.pstype = pstype
        self.first_name = first_name
        self.last_name = last_name
        self.institute = institute
        self.department = department
        self.laboratory = laboratory
        self.email = email
        self.phone = phone


class Client:
    def __init__(self,
                 host: str = "localhost",
                 port: str = "5432"):
        self.host = host
        self.port = int(port)
        self.conn: Optional[Connection] = None

    def __del__(self):
        self.disconnect()

    def connect(self):
        if not self.conn:
            self.conn = psycopg2.connect(database="mb",
                                         host=self.host, port=self.port,
                                         user="massbank", password="password")
        return self.conn

    def disconnect(self):
        if not self.conn:
            return
        self.conn.close()
        self.conn = None

    def fetch_names(self) -> Iterable[NameRecord]:
        select_names = f"""
            SELECT p.project_id AS psid, '{PROJECT}' as pstype,
                COALESCE(p.first_name, ''), COALESCE(p.last_name, ''),
                COALESCE(p.institute, ''), COALESCE(p.department, ''),
                COALESCE(p.laboratory, ''),
                COALESCE(p.email, ''), COALESCE(p.phone, '')
            FROM project AS p
            UNION
            SELECT s.study_id AS psid, '{STUDY}' as pstype,
                COALESCE(s.first_name, ''), COALESCE(s.last_name, ''),
                COALESCE(s.institute, ''), COALESCE(s.department, ''),
                COALESCE(s.laboratory, ''),
                COALESCE(s.email, ''), COALESCE(s.phone, '')
            FROM study AS s, study_status_prod
            WHERE s.study_id = study_status_prod.study_id
            AND study_status_prod.status = 1
        """
        with self.connect().cursor() as cursor:
            cursor.execute(select_names)
            for row in cursor:
                yield NameRecord(*row)
