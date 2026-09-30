function convstell_sweep(outcsv)
%CONVSTELL_SWEEP Uniform-refinement convergence on the one-surface stellarator.
%
%   Reference Taylor state with B0 from the reftaylor current ring
%   (rmin = rmaj = 2), which links the stellarator tube with 0.44 clearance
%   to the surface. Each mesh is solved with two cross-section rules for the
%   flux constraint (two independent solves):
%     'orig'  prepare_stellarator's fan from (2,0,0); the u = 0 section is not
%             star-shaped about that point, so 28 of 600 nodes lie outside
%             the domain
%     'fix'   fan from (2.4,0,0), about which the section is star-shaped
%
%   Setup: zk = 0, nu = 3*nv, nr = 12, nt = 50, ntheta = 200, tol = 1e-8.
%   Error is relative L^inf of (B - B0) over the surface, max over components.

if nargin < 1
    outcsv = fullfile(fileparts(mfilename('fullpath')), 'convdata_stellarator.csv');
end

zk     = 0;
nr     = 12;
nt     = 50;
ntheta = 200;
rmin   = 2.0;
rmaj   = 2.0;
jmag   = 1.0;
tol    = 1e-8;
xscfix = 2.4;
runs   = {5, [4 6 8 12 16]
          7, [4 6 8 12 16]
          9, [4 6 8 12]};

fid = fopen(outcsv, 'w');
fprintf(fid, ['n,p,nv,nu,npts,h,tol,xs,flux,err,errx,erry,errz,' ...
    're_alpha,im_alpha,t_solve_s\n']);
fclose(fid);
fprintf('=== convstell_sweep start %s\n  -> %s\n', datestr(now), outcsv);

for k = 1:size(runs,1)
    n = runs{k,1};
    for nv = runs{k,2}
        nu = 3*nv;
        fprintf('\n\t========\n\tn = %d, nv = %d, tol = %g\n', n, nv, tol);

        [dom, qn{1}, qw{1}] = prepare_stellarator(n, nu, nv, nr, nt);
        [~, qn{2}, qw{2}] = prepare_stellarator(n, nu, nv, nr, nt, xscfix);
        dom = dom{1};
        domparams = [n, nu, nv];
        B0 = reftaylorsurffun(dom, n, ntheta, rmin, rmaj, jmag, zk);

        flux = zeros(1,2);
        for iv = 1:2
            for i = 1:nr*nt
                B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qn{iv}{1}(:,i));
                flux(iv) = flux(iv) + B0eval(2)*qw{iv}{1}(i);
            end
        end
        fprintf('flux orig = %.15e, fix = %.15e\n', flux(1), flux(2));

        xsname = {'orig', 'fix'};
        for iv = 1:2
            t0 = tic;
            tv = RefTaylorState({dom}, domparams, zk, flux(iv), {B0}, ...
                qn{iv}, qw{iv}, tol);
            tv = tv.solve(true);
            B = tv.surface_B();
            t_solve = toc(t0);

            D = B0 - B{1};
            err = vecinfnorm(D)/vecinfnorm(B0);
            % location of the max |B0 - B|; surfacefun subsref needs
            % per-patch indexing, not vals{:}
            emax = -1;
            for ip = 1:length(dom.x)
                e = abs(D.components{1}.vals{ip}).^2 ...
                    + abs(D.components{2}.vals{ip}).^2 ...
                    + abs(D.components{3}.vals{ip}).^2;
                [em, ij] = max(e(:));
                if em > emax
                    emax = em;
                    xyz = [dom.x{ip}(ij) dom.y{ip}(ij) dom.z{ip}(ij)];
                end
            end

            fprintf('n=%d nv=%d xs=%s  err=%.6e at (%.3f,%.3f,%.3f)  alpha=%.12e%+.12ei\n', ...
                n, nv, xsname{iv}, err, xyz, real(tv.alpha), imag(tv.alpha));
            fid = fopen(outcsv, 'a');
            fprintf(fid, '%d,%d,%d,%d,%d,%.10f,%g,%s,%.15e,%.6e,%.4f,%.4f,%.4f,%.15e,%.15e,%.1f\n', ...
                n, n-1, nv, nu, nu*nv*n*n, 1/sqrt(nu*nv), tol, xsname{iv}, ...
                flux(iv), err, xyz, real(tv.alpha), imag(tv.alpha), t_solve);
            fclose(fid);
        end
    end
end

fprintf('\nDONE_SENTINEL %s\n', datestr(now));
end

function N = vecinfnorm(f)
N = max([norm(f.components{1}, inf) ...
        norm(f.components{2}, inf), ...
        norm(f.components{3}, inf)]);
end
