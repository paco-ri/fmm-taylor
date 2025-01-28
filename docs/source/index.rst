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

Introduction
============
``fmm-taylor`` is a code that solves for magnetohydrodynamic (MHD) equilibrium in stellarators. The vacuum chamber of a stellarator has a boundary that is diffeomorphic to a torus, so in this documentation, we abuse the term 'torus' to refer to the boundary of a stellarator-like geometry.

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
   \oint_{S_k} \mathbf{B} \cdot \, \text{d}\mathbf{a} &= \Phi_k, k = 1, \dots, N_k.
   \end{align}
   :label: eq:bvp

Here :math:`\mathbf{B}` is the magnetic field and :math:`\Phi_k` is the magnetic flux through :math:`S_k`. This BVP arises when computing MHD equilibrium in a stellarator according to the multi-region relaxed MHD model. When :math:`\mathbf{B}` solves :math:numref:`eq:bvp`, we say that the plasma is in a Taylor state.

This code solves :math:numref:`eq:bvp` by recasting by writing :math:`\mathbf{B}(\mathbf{x})` as 

.. math::
   \mathbf{B} = i\lambda S_\lambda[\mathbf{m}] - \nabla S_\lambda[\sigma] + i \nabla \times S_\lambda[\mathbf{m}], 
   :label: eq:B

where

.. math::
   \mathbf{m} = i\lambda (\nabla_\Gamma \Delta_\Gamma^{-1} \sigma - i \mathbf{n} \times \nabla_\Gamma \Delta_\Gamma^{-1} \sigma) + \sum_{k=1}^{N_k} \alpha_k \mathbf{m}_H^k
   :label: eq:m

and :math:`\sigma, \alpha_k, \mathbf{m}, \mathbf{n}, \mathbf{m}_H^k` are functions of :math:`\mathbf{x} \in \Omega`.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
