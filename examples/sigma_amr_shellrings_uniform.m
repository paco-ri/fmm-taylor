% Uniform refinement on the geometry and reference field of
% sigma_amr_shellrings.m: every patch is split into four at every level
% (frac = 1), with everything else as there. Writes the same outputs as
% sigma_amr_shellrings.m with BASE = sigma_amr_shellrings_uniform.
%
% Variables read if set before running: setup_only, nv, rmax and figext, as
% in sigma_amr_shellrings.m.

frac = 1;
sigma_amr_shellrings
