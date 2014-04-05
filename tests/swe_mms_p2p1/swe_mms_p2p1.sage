y = var('y')
from math import pi

def function(phi_0, phi_x, phi_y, phi_xy, 
             f_sin_x, f_cos_x, f_sin_y, f_cos_y, f_sin_xy, f_cos_xy, 
             alpha_x, alpha_y, alpha_xy):
    
    f_0 = phi_0 
    f_x = phi_x*(f_sin_x*sin(alpha_x*x) + f_cos_x*cos(alpha_x*x)) 
    f_y = phi_y*(f_sin_y*sin(alpha_y*y) + f_cos_y*cos(alpha_y*y)) 
    f_xy = phi_xy*(f_sin_xy*sin(alpha_xy*x*y/pi) + f_cos_xy*cos(alpha_xy*x*y/pi)) 
    f = f_0 + f_x + f_y + f_xy
    return f

h = sin(x)*sin(y)
u = cos(x)*sin(y)
v = sin(x*x) + cos(y)

g = 9.8

mu = 0.7
tau_xx = mu*diff(u,x)            
tau_yy = mu*diff(v,y)
tau_xy = mu*diff(u,y)
tau_yx = mu*diff(v,x)

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
