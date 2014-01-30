# Copyright 2014 Imperial College London. All rights reserved.

import vtktools
import numpy

class Detectors:
   """ A class for writing out field values at particular coordinates. """
   def __init__(self, locations_file_name, values_file_name, fields):
   
      # A list of all the fields that the detectors consider.
      self.fields = fields
         
      # The locations file contains the coordinates of all the detectors
      try:
         self.locations_file = open(locations_file_name, "r")
      except IOError:
         print "Error: Could not open detector locations file for reading. Check file path and read permissions?"
        
      # The values file contains the values of each field (in the 'fields' list) at each detector location.
      try:
         self.values_file = open(values_file_name, "w")
      except IOError:
         print "Error: Could not open detector values file for writing. Check write permissions?"
         
      # Read in all the coordinates to a list of tuples.
      self.coordinates = []
      for line in self.locations_file:
         xyz = line.split()
         if(len(xyz) == 1):
            self.coordinates.append((float(xyz[0]), 0.0, 0.0))
         elif(len(xyz) == 2):
            self.coordinates.append((float(xyz[0]), float(xyz[1]), 0.0))
         elif(len(xyz) == 3):
            self.coordinates.append((float(xyz[0]), float(xyz[1]), float(xyz[2])))
         else:
            print "Warning: Coordinates file is empty!"
            
      return

   def write(self, simulation_name, t, dt):
      
      # FIXME: Don't probe the VTU files after the fields have been written. Instead, pass in the fields and just interpolate them directly.
      vtu_files = []
      for field in self.fields:
         vtu_files.append(vtktools.vtu("%s_%s_%d.vtu" % (simulation_name, field, int(t/dt))))
         
      self.values_file.write(str(t) + " ") # The current time is the first entry in each row.
      for xyz in self.coordinates:
         for i in range(0, len(self.fields)):
            self.values_file.write(str(vtktools.vtu.ProbeData(vtu_files[i], numpy.array([xyz]), self.fields[i])[0][0]) + " ")
      self.values_file.write("\n")
      self.values_file.flush()
      return
      
   def finalise(self):
      self.locations_file.close()
      self.values_file.close()
      return
      
