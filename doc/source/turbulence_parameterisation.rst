Turbulence parameterisation
===========================

This chapter describes the turbulence models that are available in Firedrake-Fluids.

Large Eddy Simulation (LES)
---------------------------

The UFL implementation of all LES models can be found in the file
``les.py``.

Smagorinsky model
~~~~~~~~~~~~~~~~~

The model calculates an eddy viscosity :math:`\nu^\prime`

.. math:: \nu^\prime = \left(C_s\Delta_e\right)^2|\mathbb{S}|,

where :math:`C_s` is the Smagorinsky coefficient which is typically
between 0.1 and 0.2 , and :math:`\Delta_e` is some measure of the
element size. Here it is given by the square root of the element’s area
in 2D, or cube root of the element’s volume in 3D. The term
:math:`|\mathbb{S}|` is the modulus of the strain rate tensor
:math:`\mathbb{S}`\ :

.. math:: \mathbb{S} = \frac{1}{2}\left(\nabla\mathbf{u} + \nabla\mathbf{u}^{\mathrm{T}}\right).

The eddy viscosity :math:`\nu^\prime` is added to the
background/physical viscosity of the fluid, thereby contributing to the
stress term in the momentum equation.

