function [u, v, w, curlfree, divfree] = hodge_inward(f)
%HODGE_INWARD   Hodge decomposition of a SURFACEFUNV.
%   Version of the surfacefun function where the domain has an inward-
%   pointing normal

pdo = [];
pdo.lap = 1;
dom = domain(f);
n = -normal(dom);

% Curl-free component: Solve lap(u) = div(f)
L = surfaceop(dom, pdo, div(f));
L.rankdef = true;
u = solve(L);

% Divergence-free component: Solve lap(v) = -div(n x f)
L.rhs = -div(cross(n, f));
v = solve(L);

% Harmonic component: w = f - grad(u) - n x grad(v)
curlfree = grad(u);
divfree = cross(n, grad(v));
w = f - curlfree - divfree;

end