import unittest, omnical._omnical as _O

class TestImport(unittest.TestCase):
    def test_test(self):
        self.assertEqual(_O.test(1,2), 3)

if __name__ == '__main__':
    unittest.main()
