import unittest
import pandas as pd
from basic_analysis import parsefile
from utils import result_respir

class TestDataExpl(unittest.TestCase):

    @classmethod
    def setUpClass(tb):
        print("Loading dataset")
        tb._df = parse_file('../data/tb_functions.pl')

    def test_result_respir(self):
        print("Starting test_result_respir")
        self.assertEqual(result_respir(self._df, "respiration"), 0)