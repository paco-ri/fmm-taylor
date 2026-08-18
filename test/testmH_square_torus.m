%% vn / mH sanity checks on the (reparametrized) my_square_torus
%
% Checks, for the new single-ring my_square_torus(n,nu,nv):
%   1. vn is finite and unit-norm everywhere (Domain construction already
%      enforces outward orientation via computeNormals' volume check, so
%      this is mainly a regression check that the new parametrization
%      doesn't produce degenerate normals at the 4 sharp edges).
%   2. mH (Domain.mH, the surface-harmonic vector field) is finite
%      everywhere, and is continuous across patch boundaries to within a
%      tolerance that should shrink under refinement (away from the
%      sharp edges) but may saturate at the 4 edges where the geometry is
%      only C0.
%   3. mH's A-cycle circulation (via TaylorState.intacyc, same routine
%      used in mtxfluxalpha/mtxfluxsigma) is O(1) and has the same sign
%      and order of magnitude as on a smooth circular torus.

clear
addpath('../examples')

ns = [5 7];
m  = 6;

circ_nu = 3*8;
circ_nv = 8;
circ_rmin = 0.5;
circ_rmaj = 1.0;
circdom = circulartorus(5, circ_nu, circ_nv, circ_rmin, circ_rmaj);
circdomain = Domain(circdom, [5, circ_nu, circ_nv]);
circ_circ = TaylorState.intacyc(circdomain.mH{1}, 5, circ_nu, circ_nv);
fprintf('circulartorus:    mH A-cycle circulation = % .6e\n', circ_circ);

for n = ns
    nu = 3*m;
    nv = 4*m;
    dom = my_square_torus(n, nu, nv);
    domparams = [n, nu, nv];
    domain = Domain(dom, domparams);

    %% 1. vn finite, unit norm
    maxnormerr = 0;
    for k = 1:numel(dom.x)
        vnx = domain.vn{1}.components{1}.vals{k};
        vny = domain.vn{1}.components{2}.vals{k};
        vnz = domain.vn{1}.components{3}.vals{k};
        assert(~any(isnan(vnx(:)) | isinf(vnx(:))), 'vn has NaN/Inf');
        nrm = sqrt(vnx.^2 + vny.^2 + vnz.^2);
        maxnormerr = max(maxnormerr, max(abs(nrm(:) - 1)));
    end
    fprintf('n = %d, m = %d: max |vn| - 1 over all patches = %.3e\n', ...
        n, m, maxnormerr);
    assert(maxnormerr < 1e-8, 'vn is not unit norm');

    %% 1b. vn direction on the outer wall (v in [0,1], r = router = 1)
    % There, the surface normal must be purely radial: +rhat = (x,y,0)
    % (router=1) if outward, -rhat if (incorrectly) inward.
    router = 1;
    k = 1; % ku = 1, kv = 1 (outer wall)
    xk = dom.x{k};
    yk = dom.y{k};
    vnx = domain.vn{1}.components{1}.vals{k};
    vny = domain.vn{1}.components{2}.vals{k};
    vnz = domain.vn{1}.components{3}.vals{k};
    rdotvn = (xk.*vnx + yk.*vny)/router; % +1 if vn = +rhat, -1 if vn = -rhat
    fprintf('n = %d, m = %d: outer-wall vn . rhat = % .6f (max |vn_z| = %.3e)\n', ...
        n, m, mean(rdotvn(:)), max(abs(vnz(:))));
    assert(abs(mean(rdotvn(:))) > 0.99, ...
        'vn on the outer wall is not aligned with +/- rhat as expected');
    assert(mean(rdotvn(:)) > 0, ...
        'vn on the outer wall points INWARD (-rhat) instead of outward (+rhat)');

    %% 2. mH finite, continuous across patch edges
    maxjump = 0;
    maxjump_interior = 0; % excludes the 4 sharp-edge boundaries (v = 0,1,2,3)
    for ku = 1:nu
        for kv = 1:nv
            k = (ku-1)*nv + kv;
            kvnext = mod(kv, nv) + 1;
            kvn = (ku-1)*nv + kvnext;
            kunext = mod(ku, nu) + 1;
            kun = (kunext-1)*nv + kv;
            for c = 1:3
                valsk = domain.mH{1}.components{c}.vals{k};
                assert(~any(isnan(valsk(:)) | isinf(valsk(:))), 'mH has NaN/Inf');

                % v-direction neighbor (rows): last row of k vs first row of kvn
                jump_v = max(abs(valsk(end,:) - domain.mH{1}.components{c}.vals{kvn}(1,:)));
                % u-direction neighbor (cols): last col of k vs first col of kun
                jump_u = max(abs(valsk(:,end) - domain.mH{1}.components{c}.vals{kun}(:,1)));

                maxjump = max([maxjump, jump_v, jump_u]);
                if mod(kv, m) ~= 0 % not a sharp-edge boundary in v
                    maxjump_interior = max([maxjump_interior, jump_v, jump_u]);
                end
            end
        end
    end
    fprintf('n = %d, m = %d: max mH jump across patch edges = %.3e (away from sharp edges: %.3e)\n', ...
        n, m, maxjump, maxjump_interior);

    %% 3. mH A-cycle circulation
    sq_circ = TaylorState.intacyc(domain.mH{1}, n, nu, nv);
    fprintf('n = %d, m = %d: square torus mH A-cycle circulation = % .6e\n', ...
        n, m, sq_circ);
end

fprintf(['testmH_square_torus: compare the circulation values above ' ...
    '(same sign/order of magnitude expected); jumps at sharp edges ' ...
    '(v = 0,1,2,3) are expected to be larger than interior jumps but ' ...
    'should not blow up or be NaN.\n']);
