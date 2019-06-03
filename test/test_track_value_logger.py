import unittest
import os
from apal_cxx import pyread_from_file

class TestTrackValueLogger(unittest.TestCase):
    def test_read_from_file(self):
        self.assertTrue(pyread_from_file())
        os.remove("track_values_logger.txt")

if __name__ == "__main__":
    unittest.main()