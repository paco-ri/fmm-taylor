%% Near-quadrature correction self-convergence near the sharp edges
%
% Evaluates curl S0[mH] on-surface (the operator used throughout
% mtxfluxalpha/mtxBalpha/gmresA) on my_square_torus at two values of n,
% with near-quadrature corrections from
% taylor.static.get_quadrature_correction. Reports the per-patch L2 norm
% of curl S0[mH] at the lower resolution, split into patches that touch
% one of the 4 sharp edges (v = 0,1,2,3 in the cross-section
% parametrization) vs interior patches, and how those compare to the
% higher-resolution run.
%
% This isolates whether the near-singular quadrature correction near the
% sharp edges is a/the accuracy bottleneck, independently of GMRES and the
% flux/alpha bookkeeping.

clear

zk = 0;
epstaylor = 1e-3; % lower precision -> much smaller near-quadrature correction
m = 3;            % fewer patches per face -> smaller problem (nu=9, nv=12)
nu = 3*m;
nv = 4*m;

% Patches touching a sharp edge: cross-section index kv in {1, m, m+1,
% 2m, 2m+1, 3m, 3m+1, nv} (the two patches on either side of each of the
% 4 face-to-face transitions, including the v=0/v=nv wrap).
edge_kv = unique([1, m, m+1, 2*m, 2*m+1, 3*m, 3*m+1, nv]);

results = struct('n',{},'edge_norm',{},'interior_norm',{});

for n = [3 5]
    clear dom S domain mHvals targinfo opts Q curlS0mH curlS0mH_func
    dom = my_square_torus(n, nu, nv);
    S = surfer.surfacemesh_to_surfer(dom);
    domain = Domain(dom, [n, nu, nv]);

    mHvals = surfacefun_to_array(domain.mH{1}, dom, S);
    mHvals = mHvals.';

    targinfo = [];
    targinfo.r = S.r;
    opts = [];
    opts.format = 'rsc';
    Q = taylor.static.get_quadrature_correction(S, epstaylor, targinfo, opts);
    opts.precomp_quadrature = Q;

    curlS0mH = taylor.static.eval_curlS0(S, mHvals, epstaylor, targinfo, opts);
    curlS0mH_func = array_to_surfacefun(curlS0mH.', dom, S);

    edge_norm = 0;
    edge_npts = 0;
    interior_norm = 0;
    interior_npts = 0;
    for ku = 1:nu
        for kv = 1:nv
            k = (ku-1)*nv + kv;
            vals = [curlS0mH_func.components{1}.vals{k}(:), ...
                    curlS0mH_func.components{2}.vals{k}(:), ...
                    curlS0mH_func.components{3}.vals{k}(:)];
            patchnorm2 = sum(abs(vals(:)).^2);
            if ismember(kv, edge_kv)
                edge_norm = edge_norm + patchnorm2;
                edge_npts = edge_npts + numel(vals);
            else
                interior_norm = interior_norm + patchnorm2;
                interior_npts = interior_npts + numel(vals);
            end
        end
    end
    edge_norm = sqrt(edge_norm/edge_npts);
    interior_norm = sqrt(interior_norm/interior_npts);

    fprintf('n = %2d: rms |curl S0[mH]|, edge patches = %.6e, interior patches = %.6e (ratio %.3f)\n', ...
        n, edge_norm, interior_norm, edge_norm/interior_norm);

    results(end+1) = struct('n', n, 'edge_norm', edge_norm, 'interior_norm', interior_norm); %#ok<SAGROW>
end

fprintf(['\ntestquadcorr_square_torus: curl S0[mH] is the curl-free density ' ...
    'mH passed through the curl-S0 layer potential -- with no boundary ' ...
    'data it is not analytically zero, but its *edge/interior ratio* and ' ...
    'how the absolute values change with n indicate whether error is ' ...
    'concentrated at the sharp edges and whether it shrinks under ' ...
    'refinement.\n']);
