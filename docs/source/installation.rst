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

The FMM libraries are Fortran codes, so subroutines in those libraries, as well as additional Fortran subroutines written just for ``fmm-taylor``, are called from MATLAB using MEX files. If you are an expert and plan to edit the MATLAB MEX functions that wrap and call Fortran subroutines, also install `mwrap <https://github.com/zgimbutas/mwrap>`_, which is a code that significantly automates the process of creating MEX files.

Installation instructions
=========================

These are instructions for installing from source.

1. Clone/download the `source code <https://github.com/paco-ri/fmm-taylor-matlab/tree/main>`_; for example, if you have ``git``,::

     git clone https://github.com/paco-ri/fmm-taylor-matlab.git

2. Change into the source code directory.
3. Run::

     make install

4. Run::

     make matlab

5. Open MATLAB and run the setup script by running::

     setup

   in the Command Window.

Notes:

- ``make install`` and ``make matlab`` are separate commands because if a MATLAB function that calls a Fortran subroutine is edited, but the subroutine itself is unchanged, one only needs to rerun ``make matlab``.
- To avoid having to run ``setup`` every time ``fmm-taylor`` is used, consider adding the call to ``setup`` to your MATLAB startup file. To do this, in MATLAB, open ``startup.m`` by running::

    edit(fullfile(userpath,'startup.m'))

  in the command window, and add ``setup`` to this file.
