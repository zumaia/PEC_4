import unittest
from utils import numero_protein, numero_lista
from utils import parse_file, result_respir, numero_hydro
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


class TestDataExpl(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print("Loading dataset")
        cls._df = parse_file("../data/tb_functions.pl")

    def test_result_respir(self):
        print("Starting test_result_respir")
        self.assertEqual(result_respir(self._df), 0)

    def test_numero_hydro(self):
        print("Starting test_numero_hydro")
        self.assertEqual(numero_hydro(self._df), 21)

    def test_numero_protein(self):
        print("Starting test_numero_protein")
        self.assertEqual(numero_protein(self._df), 43)

    def test_numero_lista(self):
        print("Starting test_numero_lista")
        self.assertEqual(numero_lista(self._df), [123])