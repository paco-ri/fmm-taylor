function [doms, fluxnodes, fluxwts] = prepare_square_torus(n, nu, nv, scale, nt)
% prepare_square_torus Create the square-torus surface mesh and boundary flux data
%
% The geometry matches @surfacemesh/square_torus.m:
%   1. outer cylindrical ring
%   2. top annulus
%   3. inner cylindrical ring
%   4. bottom annulus

router = scale;
rinner = 0.5 * scale;
zb = -0.25 * scale;
zt = 0.25 * scale;

% Use the validated surfacemesh constructor for the final domain object.
doms = {surfacemesh.square_torus(n, nu, nv)};

% A single quadrature set is used around the square cross-section.
% The angle is fixed and the poloidal direction is parameterized by the
% boundary of the (r,z) rectangle in the x-z plane at phi = 0.
[t, tw] = chebpts(nt, [0 1], 1);

fluxnodes = cell(1,1);
fluxwts = cell(1,1);
fluxnodes{1} = zeros(3, 4*nt);
fluxwts{1} = zeros(1, 4*nt);

for i = 1:nt
    s = t(i);
    ws = tw(i);

    % outer vertical segment: (r,z) = (router, zb -> zt)
    ij = i;
    fluxnodes{1}(:,ij) = [router; 0; zb + (zt - zb)*s];
    fluxwts{1}(ij) = ws * (zt - zb);

    % top horizontal segment: (r,z) = (router -> rinner, zt)
    ij = nt + i;
    fluxnodes{1}(:,ij) = [router - (router - rinner)*s; 0; zt];
    fluxwts{1}(ij) = -ws * (router - rinner);

    % inner vertical segment: (r,z) = (rinner, zt -> zb)
    ij = 2*nt + i;
    fluxnodes{1}(:,ij) = [rinner; 0; zt - (zt - zb)*s];
    fluxwts{1}(ij) = -ws * (zt - zb);

    % bottom horizontal segment: (r,z) = (rinner -> router, zb)
    ij = 3*nt + i;
    fluxnodes{1}(:,ij) = [rinner + (router - rinner)*s; 0; zb];
    fluxwts{1}(ij) = ws * (router - rinner);
end

end