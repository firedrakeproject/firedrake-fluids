import os
import pytest
import numpy
import pylab
f = open('/home/christian/Documents/svn/scripts/rcparams.py', 'r')
exec(f.read())
import vtktools

import riemann_swe

# FreeSurfacePerturbation
pylab.figure(1)
f = vtktools.vtu("swe_dam_break_1d_FreeSurfacePerturbation_240.vtu")
h_numerical = f.GetScalarField('FreeSurfacePerturbation')
coordinates = f.GetLocations()
x = coordinates[:,0]

print len(h_numerical), len(x)

t = 60.0
h_analytical = []
for point in x:
   h_analytical.append(riemann_swe.solution(point-1000, t)[0] - 5.0)

pylab.plot(x, h_numerical, 'g-', label="Numerical (Firedrake)")
pylab.plot(x, h_analytical, 'g--', label="Analytical")
pylab.legend()
pylab.xlabel(r"$x$ (m)")
pylab.ylabel(r"$h$ (m)")
pylab.axis([0, 2000, 0, 5])
pylab.savefig('h.pdf')

# Velocity
pylab.figure(2)
f = vtktools.vtu("swe_dam_break_1d_Velocity_240.vtu")
temp = f.GetVectorField('Velocity')

u_numerical = []
for u in temp:
   u_numerical.append(u[0])

coordinates = f.GetLocations()
x = coordinates[:,0]

print len(u_numerical), len(x)
print u_numerical
t = 60.0
u_analytical = []
for point in x:
   u_analytical.append(riemann_swe.solution(point-1000, t)[1])

pylab.plot(x, u_numerical, 'g-', label="Numerical (Firedrake)")
pylab.plot(x, u_analytical, 'g--', label="Analytical")
pylab.legend()
pylab.xlabel(r"$x$ (m)")
pylab.ylabel(r"$u$ (ms$^{-1}$)")
pylab.axis([0, 2000, 0, 3])
pylab.savefig('u.pdf')

