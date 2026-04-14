% test_array_surfacefun_roundtrip.m
%
% Generate a surfacemesh.stellarator, define a surfacefun on it,
% convert to array via surfacefun_to_array, convert back via
% array_to_surfacefun, and verify the result matches the original.

% --- geometry parameters ---
n  = 6;   % polynomial order per patch
nu = 4;   % patches in poloidal direction
nv = 3*nu; % patches in toroidal direction

dom = surfacemesh.stellarator(n, nu, nv);
S   = surfer.surfacemesh_to_surfer(dom);

% --- define a smooth scalar surfacefun ---
f = surfacefun(@(x,y,z) sin(x).*cos(y) + z.^2, dom);

% --- round-trip: surfacefun -> array -> surfacefun ---
fvals  = surfacefun_to_array(f, dom, S);   % (npts x 1) double array
f_back = array_to_surfacefun(fvals, dom, S); % surfacefun

% --- check error ---
err = norm(f - f_back);
fprintf('Scalar round-trip error: %.4e\n', err);
assert(err < 1e-12, ...
    sprintf('Scalar round-trip failed: error = %.4e', err));

% --- define a smooth vector surfacefunv ---
gx = surfacefun(@(x,y,z) cos(z).*x, dom);
gy = surfacefun(@(x,y,z) sin(x).*y, dom);
gz = surfacefun(@(x,y,z) cos(y).*z, dom);
g  = surfacefunv(gx, gy, gz);

% --- round-trip: surfacefunv -> array -> surfacefunv ---
gvals  = surfacefun_to_array(g, dom, S);    % (npts x 3) double array
g_back = array_to_surfacefun(gvals, dom, S); % surfacefunv
gvals_back = surfacefun_to_array(g_back, dom, S);

% --- check error ---
errv = norm(gvals - gvals_back);
fprintf('Vector round-trip error: %.4e\n', errv);
assert(errv < 1e-12, ...
    sprintf('Vector round-trip failed: error = %.4e', errv));

fprintf('All round-trip tests passed.\n');
