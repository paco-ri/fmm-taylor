function convrotellipse_sweep(outcsv)
%CONVROTELLIPSE_SWEEP Uniform-refinement convergence on the rotating ellipse.
%
%   Geometry: prepare_rotellipse (surfacefun-plus's surfacemesh.torus, max
%   curvature radius 0.108). Reference Taylor state with B0 from the reftaylor
%   current ring centred at (0,4.4,0), radius 2.75, which links the tube with
%   1.16 clearance to the surface.
%
%   Setup: zk = 0, nu = 3*nv, nr = 12, nt = 50, ntheta = 400, tol = 1e-8.
%   Error is relative L^inf of (B - B0) over the surface, max over components.

if nargin < 1
    outcsv = fullfile(fileparts(mfilename('fullpath')), 'convdata_rotellipse.csv');
end

zk     = 0;
nr     = 12;
nt     = 50;
ntheta = 400;
rmin   = 2.75;
rmaj   = 4.4;
jmag   = 1.0;
tol    = 1e-8;
runs   = {5, [4 6 8 12 16]
          7, [4 6 8 12 16]
          9, [4 6 8 12]};

fid = fopen(outcsv, 'w');
fprintf(fid, 'n,p,nv,nu,npts,h,tol,flux,err,errx,erry,errz,re_alpha,im_alpha,sigma_inf,t_solve_s\n');
fclose(fid);
fprintf('=== convrotellipse_sweep start %s\n  -> %s\n', datestr(now), outcsv);

for k = 1:size(runs,1)
    n = runs{k,1};
    for nv = runs{k,2}
        nu = 3*nv;
        fprintf('\n\t========\n\tn = %d, nv = %d, tol = %g\n', n, nv, tol);

        [dom, qn, qw] = prepare_rotellipse(n, nu, nv, nr, nt);
        dom = dom{1};
        B0 = reftaylorsurffun(dom, n, ntheta, rmin, rmaj, jmag, zk);

        flux = 0;
        for i = 1:nr*nt
            B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qn{1}(:,i));
            flux = flux + B0eval(2)*qw{1}(i);
        end
        fprintf('flux = %.15e\n', flux);

        t0 = tic;
        ts = RefTaylorState({dom}, [n, nu, nv], zk, flux, {B0}, qn, qw, tol);
        ts = ts.solve(true);
        B = ts.surface_B();
        t_solve = toc(t0);

        D = B0 - B{1};
        err = vecinfnorm(D)/vecinfnorm(B0);
        % location of the max |B0 - B|; surfacefun subsref needs per-patch
        % indexing, not vals{:}
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

        fprintf('n=%d nv=%d  err=%.6e at (%.3f,%.3f,%.3f)  alpha=%.12e%+.12ei\n', ...
            n, nv, err, xyz, real(ts.alpha), imag(ts.alpha));
        fid = fopen(outcsv, 'a');
        fprintf(fid, '%d,%d,%d,%d,%d,%.10f,%g,%.15e,%.6e,%.4f,%.4f,%.4f,%.15e,%.15e,%.4e,%.1f\n', ...
            n, n-1, nv, nu, nu*nv*n*n, 1/sqrt(nu*nv), tol, flux, err, xyz, ...
            real(ts.alpha), imag(ts.alpha), norm(ts.sigma{1}, inf), t_solve);
        fclose(fid);
    end
end

fprintf('\nDONE_SENTINEL %s\n', datestr(now));
end

function N = vecinfnorm(f)
N = max([norm(f.components{1}, inf) ...
        norm(f.components{2}, inf), ...
        norm(f.components{3}, inf)]);
end
