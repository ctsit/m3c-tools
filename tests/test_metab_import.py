import io
import unittest

import metab_import


class TestMetabImport(unittest.TestCase):
    def test_print_to_open_file_escapes(self):
        with io.StringIO() as file:
            metab_import.print_to_open_file([sentence], file)
            actual = file.getvalue()
        expected = 'THIS IS "multiple\\n lines" .\n'
        self.assertEqual(expected, actual)


sentence = '''THIS IS "multiple
 lines"'''

if __name__ == "__main__":
    unittest.main()
