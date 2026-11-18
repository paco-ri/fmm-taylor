% Uniform refinement on the geometry and reference field of
% sigma_amr_torusring.m: every patch is split into four at every level
% (frac = 1), with everything else as there. Writes the same outputs as
% sigma_amr_torusring.m with BASE = sigma_amr_torusring_uniform.
%
% Variables read if set before running: setup_only, nv, rmax and figext, as
% in sigma_amr_torusring.m.

frac = 1;
sigma_amr_torusring
