% Two-surface RefTaylorState solve on UNEQUAL, mixed-level meshes.
%
% Refines only the OUTER surface, so the two surfaces carry different numbers
% of points.
% 
% Reference errors vs B0 (tol = 1e-8):
%   uniform (control)                 [675 675]   1.7939e-01
%   outer refined ALL, inner none     [2700 675]  5.2250e-02
%   outer refined PART, inner none    [1725 675]  1.3168e-01
fprintf('=== test_asymmetric start ===\n');

% geometry parameters
n = 5; nv = 3; nu = 3*nv; ao = 1.0; ai = 0.6;
domparams = [n, nu, nv];
zk = 0.5; ntheta = 200; rmin = 2.0; rmaj = 2.0; jmag = 1.0;
nr = 16; nt = 40; np = 40;
tol = 1e-8;

% two toroidal domains, quadrature nodes and weights for fluxes
[dom, qnodes, qweights] = prepare_torus(n,nu,nv,n,nu,nv,ao,ai,nr,nt,np);

% compute fluxes from reference magnetic field B0e
flux = zeros(1,2);
for i = 1:nr*nt
    B0e = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes{1}(:,i));
    flux(1) = flux(1) + B0e(2)*qweights{1}(i);
end
for i = 1:nr*np
    B0e = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes{2}(:,i));
    flux(2) = flux(2) - B0e(3)*qweights{2}(i);
end
fprintf('flux = [%.10e %.10e]\n', flux(1), flux(2));

P = struct('n',n,'nu',nu,'nv',nv,'ntheta',ntheta,'rmin',rmin,'rmaj',rmaj, ...
    'jmag',jmag,'zk',zk,'domparams',domparams,'flux',flux, ...
    'qnodes',{qnodes},'qweights',{qweights},'tol',tol);

rmax = 1; % max quadtree depth
np1 = length(dom{1}.x); np2 = length(dom{2}.x); % number of patches on each surface
p1 = [(1:np1).' zeros(np1,2)]; 
p2 = [(1:np2).' zeros(np2,2)];

% different sets of marked patches
cases = { {[],            'uniform (control)'}, ...
          {(1:np1).',     'outer refined ALL, inner none'}, ...
          {(1:2:np1).',   'outer refined PART, inner none'} };

dm = cell(1,2); qf = cell(1,2); pq = cell(1,2);
% store domain and quadforest information for the inner surface (no refinement)
[dm{2}, qf{2}, pq{2}] = surfacemesh.refine_leaves(dom{2}, p2, [], rmax);

for c = 1:numel(cases)
    mk1 = cases{c}{1}; tag = cases{c}{2};
    % store domain and quadforest information for the outer surface (refined)
    [dm{1}, qf{1}, pq{1}] = surfacemesh.refine_leaves(dom{1}, p1, mk1, rmax);
    try
        [e, D, alpha] = referr(dm, qf, pq, P); % relative B error, domain, alpha
        lv1 = mat2str(unique(pq{1}(:,2)).');
        fprintf(['%-32s npat=%3d+%-3d npts=%s lvls(out)=%-7s ' ...
                 'err=%.4e alpha=[%.6e %.6e]\n'], tag, ...
            length(dm{1}.x), length(dm{2}.x), mat2str(D.nptspersurf), ...
            lv1, e, real(alpha(1)), real(alpha(2)));
    catch ME
        fprintf('%-32s FAILED: %s\n', tag, ME.message);
        for q=1:numel(ME.stack), fprintf('    at %s:%d\n', ME.stack(q).name, ME.stack(q).line); end
    end
end
fprintf('DONE_SENTINEL\n');

function [e, D, alpha] = referr(dd, qf, p2q, P)
B0 = cell(1,2);
for s = 1:2
    B0{s} = reftaylorsurffun(dd{s}, P.n, P.ntheta, P.rmin, ...
        P.rmaj, P.jmag, P.zk);
end
D = Domain(dd, P.domparams, qf, p2q);
ts = RefTaylorState(D, P.domparams, P.zk, P.flux, B0, P.qnodes, P.qweights, P.tol);
ts = ts.solve(false);
B = ts.surface_B();
num = 0; den = 0;
for s = 1:2
    num = max(num, vecinfnorm(B0{s} - B{s}));
    den = max(den, vecinfnorm(B0{s}));
end
e = num/den; alpha = ts.alpha;
end
