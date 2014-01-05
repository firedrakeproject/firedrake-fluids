import pylab

dx = [0.2, 0.1, 0.05]
h_error = [0.00348261433835, 0.00113851249434, 0.000398910404953] # Error in L2 norm
first_order = [1.0e-2*(1.0/2.0)**i for i in range(0, len(h_error))]
second_order = [1.0e-2*(1.0/4.0)**i for i in range(0, len(h_error))]

pylab.loglog(dx, h_error, label="Free surface height")
pylab.loglog(dx, first_order, label="First-order convergence")
pylab.loglog(dx, second_order, label="Second-order convergence")
pylab.legend()
pylab.xlabel("Characteristic element length")
pylab.ylabel("Error in L2 norm")
pylab.axis([0.01, 1.0, 0.0001, 0.1])
pylab.savefig("swe_mms_p2p1_error.png")
