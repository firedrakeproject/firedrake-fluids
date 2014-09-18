Stabilisation methods
=====================

[chap:stabilisation] When using a continuous Galerkin discretisation in
advection-dominated problems, it may be necessary to stabilise the
advection term in the momentum equation.

The implementation of the stabilisation methods can be found in the file
``stabilisation.py``.

Streamline upwind
-----------------

This method adds some upwind diffusion in the direction of the
streamlines. The term is given by

.. math:: \int_{\Omega} \frac{\bar{k}}{||\mathbf{u}||^2}(\mathbf{u}\cdot\nabla\mathbf{w})(\mathbf{u}\cdot\nabla\mathbf{u})

which is added to the LHS of the momentum equation. The term
:math:`\bar{k}` takes the form

.. math:: \bar{k} = \frac{1}{2}\left(\frac{1}{\tanh(\mathrm{Pe})} - \frac{1}{\mathrm{Pe}}\right)||\mathbf{u}||\Delta x

where

.. math:: \mathrm{Pe} = \frac{||\mathbf{u}||\Delta x}{2\nu}

is the Peclet number, and :math:`\Delta x` is the size of each element.

