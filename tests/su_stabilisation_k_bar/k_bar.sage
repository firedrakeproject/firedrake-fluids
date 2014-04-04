y = var('y')
cellsize = var('cellsize')
nu = var('nu')

from math import pi

u = x
v = y

magnitude = sqrt(u*u + v*v)
print str(magnitude).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')

grid_pe = (magnitude*cellsize)/(2.0*nu)

k_bar = ( (1.0/tanh(grid_pe)) - (1.0/grid_pe) ) * cellsize * magnitude
print str(k_bar).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')
