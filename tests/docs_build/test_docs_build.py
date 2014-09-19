import os, sys
import subprocess
import pytest

def docs_build():
   sys.path.insert(0, os.path.abspath('.'))
   print sys.path

   with open("docs_build.log", "w") as out:
      process = subprocess.Popen("cd docs; make clean; make html; cd ../", shell=True, stdout=out, stderr=sys.stdout.fileno())
      exit_code = process.wait()
   
   with open("docs_build.log", "r") as f:
      log = f.read()
   
   return log, exit_code

def test_docs_build():
   log, exit_code = docs_build()
   
   print log
   print exit_code
   
   assert (exit_code == 0)
   assert ("build succeeded" in log) # Ignoring warnings for now.
   #assert ("build succeeded." in log) # Note that the full stop is needed here, otherwise warning messages will not be noticed.
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
