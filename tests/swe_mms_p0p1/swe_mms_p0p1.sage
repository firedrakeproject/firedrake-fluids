y = var('y')
from math import pi

h = sin(x)*sin(y)
u = cos(x)*sin(y)
v = sin(x*x) + cos(y)

g = 9.8

print (diff(u,x) + diff(v,y))

nu = 0.7
tau_xx = nu*diff(u,x)
tau_xy = nu*diff(u,y)
tau_yy = nu*diff(v,y)
tau_yx = nu*diff(v,x)

Su = u*diff(u,x) + v*diff(u,y) + g*diff(h,x) - diff(tau_xx, x) - diff(tau_xy, y)
Sv = u*diff(v,x) + v*diff(v,y) + g*diff(h,y) - diff(tau_yy, y) - diff(tau_yx, x) 

h_mean = 50.0
Sh = diff((h_mean + h)*u, x) + diff((h_mean + h)*v, y)

print 'from math import sin, cos, tanh, pi, sqrt'
print ''
print 'def u(X):'
print '    return', str(u).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')
print ''
print 'def v(X):'
print '    return', str(v).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')
print ''  
print 'def h(X):'
print '    return', str(h).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')
print '' 
print 'def forcing_u(X):'
print '    return', str(Su).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')
print ''
print 'def forcing_v(X):'
print '    return', str(Sv).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')
print ''
print 'def forcing_h(X):'
print '    return', str(Sh).replace('e^', 'exp').replace('^', '**').replace('000000000000', '').replace('x', 'x[0]').replace('y', 'x[1]')
print ''
print 'def velocity(X):'
print '   return [u(X), v(X)]'
print ''
print 'def forcing_velocity(X):'
print '   return [forcing_u(X), forcing_v(X)]'
print ''
