#!/usr/bin/env python

import sys
import subprocess

def setup(project):

   procs = [1, 2, 4, 8, 16, 24, 48, 96, 192, 384, 768]
   
   for p in procs:
      if(p <= 24):
         nodes = 1
         procs_per_node = p
      else:
         nodes = int(p/24)
         procs_per_node = 24
         
      subprocess.call('mkdir %s' % "p"+str(p), shell=True)
      subprocess.call('cp src/*.msh %s' % "./p"+str(p), shell=True)
      
      src = open("src/submit_archer.pbs", "r")
      out = open("p"+str(p)+"/submit_archer.pbs", "w")
      
      for line in src:
         if("__JOBNAME__" in line):
            s = line.replace("__JOBNAME__", "p"+str(p))
         elif("__PROJECT__" in line):
            s = line.replace("__PROJECT__", project)
         elif("__NODES__" in line):
            s = line.replace("__NODES__", str(nodes))
         elif("__PROCS__" in line):
            s = line.replace("__PROCS__", str(p))
            s = s.replace("__MAXPROCS__", str(procs_per_node))
         else:
            s = line
         out.write(s)
      
      src.close()
      out.close()
       
if __name__ == "__main__":
   setup(sys.argv[1])
