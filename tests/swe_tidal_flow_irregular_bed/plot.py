import os
import numpy
import pylab
import vtktools


def bed_height(x):
   if(x <= 50):
      return 0
   elif(x <= 100):
      return (x-50)*(2.5 - 0.0)/50.0
   elif(x <= 150):
      return 2.5 + (x-100)*(5.0-2.5)/50.0
   elif(x <= 250):
      return 5.0
   elif(x <= 300):
      return 5.0 + (x-250)*(3.0-5.0)/50.0
   elif(x <= 350):
      return 3.0 + (x-300)*(5.0-3.0)/50.0 
   elif(x <= 400):
      return 5
   elif(x <= 425):
      return 5.0 + (x-400)*(7.5-5.0)/25.0
   elif(x <= 435):
      return 7.5 + (x-425)*(8.0-7.5)/10.0
   elif(x <= 450):
      return 8.0 + (x-435)*(9.0-8.0)/15.0
   elif(x <= 475):
      return 9
   elif(x <= 500):
      return 9.0 + (x-475)*(9.1-9.0)/25.0
   elif(x <= 505):
      return 9.1 + (x-500)*(9.0-9.1)/5.0
   elif(x <= 530):
      return 9
   elif(x <= 550):
      return 9.0 + (x-530)*(6.0-9.0)/20.0
   elif(x <= 565):
      return 6.0 + (x-550)*(5.5-6.0)/15.0
   elif(x <= 575):
      return 5.5
   elif(x <= 600):
      return 5.5 + (x-575)*(5.0-5.5)/25.0
   elif(x <= 650):
      return 5.0 + (x-600)*(4.0-5.0)/50.0
   elif(x <= 700):
      return 4.0 + (x-650)*(3.0-4.0)/50.0
   elif(x <= 750):
      return 3
   elif(x <= 800):
      return 3.0 + (x-750)*(2.3-3.0)/50.0
   elif(x <= 820):
      return 2.3 + (x-800)*(2.0-2.3)/20.0
   elif(x <= 900):
      return 2.0 + (x-820)*(1.2-2.0)/80.0
   elif(x <= 950):
      return 1.2 + (x-900)*(0.4-1.2)/50.0
   elif(x <= 1000):
      return 0.4 + (x-950)*(0.0-0.4)/50.0
   else:
      return 0
      
L = 1500.0
t = 10800.0

x_array = numpy.linspace(0, 1500, 201)

f = vtktools.vtu('swe_tidal_flow_irregular_bed_Velocity_36000.vtu')
f.GetFieldNames()
temp = f.GetVectorField('Velocity')

u_numerical = []
for i in range(len(temp)):
   u_numerical.append(temp[i][0])
print len(u_numerical)
print len(x_array)
h_analytical = []
u_analytical = []
for x in x_array:
   b = bed_height(x)

   H = 20.0 - b
   
   h = H - 4*numpy.sin(numpy.pi*(4*t/86400.0 + 0.5))
   u_analytical.append( ((x - L)*numpy.pi/(5400.0*h))*numpy.cos(numpy.pi*(4*t/86400.0 + 0.5)) )
   h_analytical.append(h)


pylab.figure(0)
pylab.plot(x_array, u_analytical, 'g--', label="Analytical")
pylab.plot(x_array, u_numerical, 'g-', label="Numerical (Firedrake-Fluids)")
pylab.legend()
pylab.xlabel("x (m)")
pylab.ylabel("Velocity (m/s)")
pylab.savefig('swe_tidal_flow_irregular_bed_velocity.png')


