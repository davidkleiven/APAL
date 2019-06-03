import unittest
import os
from apal_cxx import pyread_from_file, pyinit_keys_from_entry

class TestTrackValueLogger(unittest.TestCase):
    def test_read_from_file(self):
        self.assertTrue(pyread_from_file())
        os.remove("track_values_logger.txt")

    def test_init_keys_from_entry(self):
        pyinit_keys_from_entry()

        expected = [
            "# Iter,Value1,Value2,Value3\n",
            "1,0.1,-1.2,0.4\n"
        ]

        with open("init_key_from_entry.txt", 'r') as infile:
            lines = infile.readlines()
        self.assertEqual(lines, expected)
        os.remove("init_key_from_entry.txt")
if __name__ == "__main__":
    unittest.main()