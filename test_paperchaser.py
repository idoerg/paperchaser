
import sys, os
sys.path.append(os.path.dirname(__file__))
import io
import paperchaser as pc
import unittest
class CompFiles(unittest.TestCase):
    def compare_files(self,tst_path,ref_path):
        self.assertListEqual(list(io.open(tst_path)),
            list(io.open(ref_path)))

def test_paperchaser():
    compfiles = CompFiles()
    fout = open("foo.txt","w")
    pc.main(["Iddo Friedberg"], affiliation=None,
            email="idoerg@iastate.edu", years=[2015, 2017], outfile=fout,
            exclude_file=None,conflicts=None,verbose=1)
    fout.close()
    compfiles.compare_files(fout.name, "bar.txt")

