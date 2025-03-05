************
Installation
************

This code is currently only supported for Linux.

Dependencies
============

Before installing ``fmm-taylor``, install the following packages:

- `fmm3d <https://fmm3d.readthedocs.io/en/latest/>`_, the FMM library,
- `fmm3dbie <https://fmm3dbie.readthedocs.io/en/latest/index.html>`_, the FMM boundary integral solver library, and
- `this fork <https://github.com/paco-ri/surfacefun-plus>`_ of `surfacefun <https://surfacefun.readthedocs.io/en/latest/index.html>`_, the surface PDE solver.

Installation instructions
=========================

These are instructions for installing from source.

1. Clone/download the `source code <https://github.com/paco-ri/fmm-taylor-matlab/tree/main>`_; for example, if you have ``git``,::

     git clone https://github.com/paco-ri/fmm-taylor-matlab.git

2. Change into ``fmm-taylor``
3. Run::

     make install

4. Run::

     make matlab

