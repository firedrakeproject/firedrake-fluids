import os, sys
import subprocess
import pytest

cwd = os.path.dirname(os.path.abspath(__file__))
print cwd

def manual():
   process = subprocess.Popen("make manual", shell=True, stdout=subprocess.PIPE, stderr=sys.stdout.fileno())
   log = process.stdout.read()
   exit_code = process.wait()
   return log, exit_code

def test_manual():
   log, exit_code = manual()
   
   print log
   print exit_code
   
   assert (exit_code == 0)
   assert (not "There were undefined citations." in log)
   assert (not "There were undefined references." in log)
   assert (not "There were multiply-defined labels." in log)
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
