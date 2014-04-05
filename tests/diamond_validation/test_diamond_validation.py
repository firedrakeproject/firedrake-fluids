import os
import pytest
import numpy
from firedrake import *

def diamond_validation():
   # This test checks the validity of the simulation configuration files in the 'tests' directory
   # (in XML format) against the corresponding schema.
   # Note: most of the code for this test is from the diamond_validation test for the Fluidity CFD code,
   # also developed at Imperial College London.

   import glob
   import os
   import sys
   import xml.dom.minidom
   import xml.parsers.expat

   import diamond.debug as debug
   import diamond.schema as schema

   debug.SetDebugLevel(0)

   class DiamondValidator:
     def __init__(self, rootDir):
       self._rootDir = rootDir
       self.Reset()
       
       return
       
     def Reset(self):
       self._passes = 0
       self._optionErrors = {}
       
       return
       
     def TestXmlFiles(self, testDir, depth):
       debug.dprint("Checking xml files:", 0)
       for filename in self._TestFiles("xml", testDir, depth):
         try:
           xmlParse = xml.dom.minidom.parse(filename)
           debug.dprint(filename + " : Pass", 0)
           self._passes += 1
         except xml.parsers.expat.ExpatError:
           debug.dprint(filename + " : Fail", 0)
           self._optionErrors[filename] = xml.parsers.expat.ExpatError
       
       return
       
     def ValidateOptionsFiles(self, schemafile, testDir, depth, extension = None, xmlRootNode = None):
       debug.dprint("Validating options file against schema: " + schemafile, 0)
     
       schemafile = os.path.join(self._rootDir, schemafile)
       sch = schema.Schema(schemafile)

       if not extension is None:
         debug.dprint("Testing files with extension: " + extension, 0)
         for filename in self._TestFiles(extension, testDir, depth): 
           optionsTree = sch.read(filename)
           lost_eles, added_eles, lost_attrs, added_attrs = sch.read_errors()
           if len(lost_eles) + len(added_eles) + len(lost_attrs) + len(added_attrs) == 0 and optionsTree.valid:
             debug.dprint(filename + " : Pass", 0)
             self._passes += 1
           else:
             debug.dprint(filename + " : Fail", 0)
             self._optionErrors[filename] = (lost_eles, added_eles, lost_attrs, added_attrs)
             
       if not xmlRootNode is None:
         debug.dprint("Testing xml files with root node: " + xmlRootNode, 0)
         for filename in self._TestFiles("xml", testDir, depth):
           try:
             xmlParse = xml.dom.minidom.parse(filename)
           except xml.parsers.expat.ExpatError:
             continue
           rootEles = xmlParse.getElementsByTagName(xmlRootNode)
           if len(rootEles) == 0:
             continue
           optionsTree = sch.read(filename)
           lost_eles, added_eles, lost_attrs, added_attrs = sch.read_errors()
           if len(lost_eles) + len(added_eles) + len(lost_attrs) + len(added_attrs) == 0 and optionsTree.valid:
             debug.dprint(filename + " : Pass", 0)
             self._passes += 1
           else:
             debug.dprint(filename + " : Fail", 0)
             self._optionErrors[filename] = (lost_eles, added_eles, lost_attrs, added_attrs)
       
       return
       
     def _TestFiles(self, extension, testDir, depth):
       filenames = []
       baseDir = os.path.join(self._rootDir, testDir)
       for i in range(depth + 1):
         filenames += glob.glob(os.path.join(baseDir, "*." + extension))
         baseDir = os.path.join(baseDir, "*")
       
       return filenames
       
     def Passes(self):
       return self._passes
       
     def OptionErrors(self):
       return self._optionErrors

   validator = DiamondValidator(rootDir = os.path.join(os.path.pardir, os.path.pardir))

   # Shallow water model-related tests in Firedrake-Fluids.
   validator.TestXmlFiles(testDir = "tests", depth = 2)
   validator.ValidateOptionsFiles(schemafile = os.path.join("schema", "shallow_water.rng"), testDir = "tests", depth = 2, extension = "swml")
   
   passes = validator.Passes()
   optionErrors = validator.OptionErrors()
       
   print "Summary of options files with failures:"
   for filename in optionErrors.keys():
     print filename
   print "Passes: " + str(passes)
   print "Failures: " + str(len(optionErrors))
   
   return optionErrors

def test_diamond_validation():
   optionErrors = diamond_validation()
   print optionErrors.keys()
   
   failures = []
   for filename in optionErrors.keys():
     failures.append(filename)
   print "Summary of options files with failures:"
   for filename in failures:
     print filename
   print "Failures: " + str(len(failures))
   assert(len(failures) == 0)
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))


