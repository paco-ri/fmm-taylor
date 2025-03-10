.. fmm-taylor documentation master file, created by
   sphinx-quickstart on Tue Jan 21 12:34:17 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

############################
``fmm-taylor`` documentation
############################

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   classes
   methods
   examples

Introduction
============
``fmm-taylor`` is a code that solves for magnetohydrodynamic (MHD) equilibrium in stellarators using boundary integral methods and the fast multipole method (FMM). The vacuum chamber of a stellarator has a boundary that is diffeomorphic to a torus, so in this documentation, we abuse the term 'torus' to refer to the boundary of a stellarator-like geometry.

Let

- :math:`\Omega` be the interior of one torus (:math:`N_k=1`) or two nested tori (:math:`N_k=2`),
- :math:`\Gamma := \partial \Omega`,
- :math:`\mathbf{n}` the outward unit normal to :math:`\Omega`,
- :math:`S_k` a cross-section of :math:`\Omega`,
- :math:`\lambda \geq 0`, and
- :math:`\Phi_k \in \mathbb{R}` for each :math:`k`.

``fmm-taylor`` is a MATLAB code that solves the boundary-value problem (BVP)

.. math::
   \begin{align}
   \nabla \times \mathbf{B} &= \lambda \mathbf{B} \ \text{in} \ \Omega \\
   \mathbf{B} \cdot \mathbf{n} &= 0 \ \text{on} \ \Gamma \\
   \oint_{S_k} \mathbf{B} \cdot \, \text{d}\mathbf{a} &= \Phi_k, \quad k = 1, \dots, N_k.
   \end{align}
   :label: eq:bvp

Here :math:`\mathbf{B}` is the magnetic field and :math:`\Phi_k` is the magnetic flux through :math:`S_k`. This BVP arises when computing MHD equilibrium in a stellarator according to the multi-region relaxed MHD model. When :math:`\mathbf{B}` solves :math:numref:`eq:bvp`, we say that the plasma is in a Taylor state.

An integral representation for the Taylor state BVP
---------------------------------------------------

The method that this code uses to solve :math:numref:`eq:bvp` arises by writing the equations in the form of integrals over the boundary :math:`\Gamma`. To start, :math:`\mathbf{B}(\mathbf{x}), \mathbf{x} \in \Omega` can be written as 

.. math::
   \mathbf{B} = i\lambda \mathcal{S}_\lambda[\mathbf{m}] - \nabla \mathcal{S}_\lambda[\sigma] + i \nabla \times \mathcal{S}_\lambda[\mathbf{m}], 
   :label: eq:B

where

.. math::
   \mathbf{m} = i\lambda (\nabla_\Gamma \Delta_\Gamma^{-1} \sigma - i \mathbf{n} \times \nabla_\Gamma \Delta_\Gamma^{-1} \sigma) + \sum_{k=1}^{N_k} \alpha_k \mathbf{m}_H^k, 
   :label: eq:m

:math:`\mathcal{S}_\lambda[f]` is the Helmholtz single-layer potential of a scalar or vector function :math:`f`, given by

.. math::
   \mathcal{S}_\lambda[f](\mathbf{x}) = \int_\Gamma f(\mathbf{y}) \frac{e^{i\lambda|\mathbf{x}-\mathbf{y}|}}{4\pi|\mathbf{x}-\mathbf{y}|} \, \text{d}s(\mathbf{y})
   :label: eq:singlelayer

and :math:`\sigma, \alpha_k, \mathbf{m}, \mathbf{n}, \mathbf{m}_H^k` are functions of :math:`\mathbf{x} \in \Omega`. In particular, :math:`\left\{m_H^k\right\}_{k=1}^{N_k}` comprise a basis of harmonic vector fields on the surface :math:`\Gamma`. :math:`\mathbf{B}` given by :math:numref:`eq:B` automatically solves the PDE :math:`\nabla \times \mathbf{B} = \lambda \mathbf{B}`, so we solve for :math:`\sigma` and :math:`\alpha := (\alpha_1, \dots, \alpha_k)` by substituting :math:numref:`eq:B` into the boundary condition of :math:numref:`eq:bvp`. This requires taking a limit of :math:`\mathbf{B}` as :math:`\mathbf{x} \to \Gamma, \mathbf{x} \in \Omega`. This limit is

.. math::
   \mathbf{B} = \frac{-\sigma}{2} \mathbf{n} + i \frac{\mathbf{n} \times \mathbf{m}}{2} + i\lambda\mathcal{S}_\lambda[\mathbf{m}] - \nabla \mathcal{S}_\lambda[\sigma] + i\nabla \times \mathcal{S}_\lambda[\mathbf{m}], 
   :label: eq:limB

evaluated at :math:`\mathbf{x} \in \Gamma`, and taking a dot product with :math:`\mathbf{n}`, the boundary condition becomes

.. math::
   0 = \frac{-\sigma}{2} + i\lambda\mathbf{n} \cdot \mathcal{S}_\lambda[\mathbf{m}] - \mathbf{n} \cdot \nabla \mathcal{S}_\lambda[\sigma] + i \mathbf{n} \cdot \nabla \times \mathcal{S}_\lambda[\mathbf{m}].
   :label: eq:ndotB

Finally, to convert the flux conditions in :math:numref:`eq:bvp` into boundary integrals, we carry out different procedures depending on the value of :math:`\lambda`. If :math:`\lambda \neq 0`, then by the fact that :math:`\mathbf{B} = \nabla \times \mathbf{B}/\lambda` and Stokes' theorem,

.. math::
   \oint_{S_k} \mathbf{B} \cdot \, \text{d}\mathbf{a} = \frac{1}{\lambda} \oint_{S_k} \nabla \times \mathbf{B} \cdot \, \text{d}\mathbf{a} = \frac{1}{\lambda} \oint_{\partial S_k} \mathbf{B} \cdot \text{d}\boldsymbol{\ell},
   :label: eq:stokes

where the last integral is a line integral over the boundary of :math:`S_k`, and we can substitute :math:`\mathbf{B}` using :math:numref:`eq:B` to obtain the flux condition in the form of a boundary integral involving :math:`\sigma` and :math:`\alpha`. If :math:`\lambda = 0`, we observe that :math:`\mathbf{B} = \nabla \times \mathcal{S}_0[\mathbf{n} \times \mathbf{B}]` and use Stokes' theorem to write

.. math::
   \oint_{S_k} \mathbf{B} \cdot \, \text{d}\mathbf{a} = \oint_{S_k} \nabla \times \mathcal{S}_0[\mathbf{n} \times \mathbf{B}] \cdot \, \text{d}\mathbf{a} = \oint_{\partial S_k} \mathcal{S}_0[\mathbf{n} \times \mathbf{B}] \cdot \text{d}\boldsymbol{\ell}.
   :label: eq:stokes0

The integral equation we set out to solve, when :math:`\lambda \neq 0`, is

.. math::
   \begin{align}
   \frac{-\sigma}{2} + i\lambda\mathbf{n} \cdot \mathcal{S}_\lambda[\mathbf{m}] - \mathbf{n} \cdot \nabla \mathcal{S}_\lambda[\sigma] + i \mathbf{n} \cdot \nabla \times \mathcal{S}_\lambda[\mathbf{m}] &= 0 \\
   \oint_{\partial S_k} \left(i\mathcal{S}_\lambda[\mathbf{m}] - \lambda^{-1} \nabla \mathcal{S}_\lambda[\sigma] + i \lambda^{-1} \nabla \times \mathcal{S}_\lambda[\mathbf{m}]\right) \cdot \text{d}\boldsymbol{\ell} &= \Phi_k, \quad k = 1, \dots, N_k.
   \end{align}
   :label: eq:bie

We can write this informally as the "linear system"

.. math::
   \begin{bmatrix}
   A_{11} & A_{12} \\
   A_{21} & A_{22}
   \end{bmatrix} \begin{bmatrix}
   \sigma \\
   \alpha
   \end{bmatrix} = \begin{bmatrix}
   0 \\
   \Phi_k
   \end{bmatrix},
   :label: eq:linsys

where the entries of the matrix are the operators

.. math::
   \begin{align}
   A_{11} &:= \frac{-I}{2} - \mathbf{n} \cdot \nabla \mathcal{S}_\lambda - \lambda^2 \mathbf{n} \cdot \mathcal{S}_\lambda \circ (\nabla_\Gamma \Delta_\Gamma^{-1} - i \mathbf{n} \times \nabla_\Gamma \Delta_\Gamma^{-1}) \\
   &\qquad - \lambda \mathbf{n} \cdot (\nabla \times \mathcal{S}_\lambda) \circ (\nabla_\Gamma \Delta_\Gamma^{-1} - i \mathbf{n} \times \nabla_\Gamma \Delta_\Gamma^{-1}), \\
   A_{12} &:= i\lambda \mathbf{n} \cdot \mathcal{S}_\lambda [\mathbf{m}_H^k] + i \mathbf{n} \cdot \nabla \times \mathcal{S}_k[\mathbf{m}_H^k], \\
   A_{21} &:= \oint_{\partial S_k} \left(-\lambda \mathcal{S}_\lambda \circ (\nabla_\Gamma \Delta_\Gamma^{-1} \sigma - i \mathbf{n} \times \nabla_\Gamma \Delta_\Gamma^{-1} \sigma) - \lambda^{-1} \nabla \mathcal{S}_\lambda \right. \\
   &\qquad \left. - (\nabla \times \mathcal{S}_\lambda) \circ (\nabla_\Gamma \Delta_\Gamma^{-1} \sigma - i \mathbf{n} \times \nabla_\Gamma \Delta_\Gamma^{-1} \sigma) \right) \cdot \text{d}\boldsymbol{\ell}, \\
   A_{22} &:= \oint_{\partial S_k} \left(i\mathcal{S}_\lambda[\mathbf{m}_H^k] + i \lambda^{-1} \nabla \times \mathcal{S}_\lambda[\mathbf{m}_H^k]\right) \cdot \text{d}\boldsymbol{\ell}.
   \end{align}

Overview of the algorithm (put this at the end)
-----------------------------------------------

Equipped with the boundary integral representation of :math:numref:`eq:bvp`, we summarize the algorithm for computing its solution. Given a toroidal domain :math:`\Omega` with boundary :math:`\Gamma`,

1. do a surface PDE solve to find :math:`\mathbf{m}_H`;
2. discretize the BIE and flux condition, writing these as

.. math::
   \begin{bmatrix}
   A_{11} & A_{12} \\
   A_{21} & A_{22}
   \end{bmatrix} \begin{bmatrix}
   \boldsymbol{\sigma} \\
   \alpha
   \end{bmatrix} = \begin{bmatrix}
   \mathbf{0} \\
   \Phi
   \end{bmatrix};

3. solve this system using GMRES, which requires an FMM and a surface PDE solve at each iteration; and
4. use :math:`\boldsymbol{\sigma}` and :math:`\alpha` to evaluate :math:`\mathbf{B}` at any point in :math:`\overline{\Omega}`.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
