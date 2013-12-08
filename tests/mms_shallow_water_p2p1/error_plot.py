import pylab

dx = [0.1, 0.05, 0.025]
h_error = [1.176096, 0.405122, 0.191630] # Error in L2 norm
first_order = [(1.0/2.0)**i for i in range(0, len(h_error))]
second_order = [(1.0/4.0)**i for i in range(0, len(h_error))]

pylab.loglog(dx, h_error, label="Free surface height")
pylab.loglog(dx, first_order, label="First-order convergence")
pylab.loglog(dx, second_order, label="Second-order convergence")
pylab.legend()
pylab.xlabel("Characteristic element length")
pylab.ylabel("Error in L2 norm")
pylab.axis([0.01, 1, 0.01, 10])
pylab.savefig("mms.png")
