function convoop_sweep(outcsv)
%CONVOOP_SWEEP One-surface convergence sweep, reproducing convdata_onesurface.csv.
%
%   Corrected form of convoop.m. The differences that matter:
%     - three orders (n = 5, 7, 9), five resolutions (nv = 6, 8, 10, 12, 16)
%     - tolerance is an explicit per-order list, not a *1e-2 ladder: n = 5 runs
%       at 1e-6 only, n = 7 and n = 9 at both 1e-6 and 1e-8
%     - ntheta = 200 and nt = 50 (convoop.m still has 1e3 and 100)
%     - per-run time, GMRES iteration counts and |sigma|_inf are recorded
%
%   Setup: prepare_torus, zk = 0, nu = 3*nv, nr = 12, ntheta = 200,
%   rmin = rmaj = 2, jmag = 1. Error is relative L^inf of (B - B0) over the
%   surface, max over components.
%
%   25 solves, about 59 min total at OMP_NUM_THREADS = 4 (4.2 GB peak at the
%   largest configuration). Rows are appended as they finish, so the run can be
%   interrupted without losing completed work.
%
%   time_s is the solve plus surface_B time, excluding geometry and B0 setup.
%   The original CSV's definition was not recorded; the full breakdown is
%   printed per run if the two need reconciling.

if nargin < 1
    outcsv = fullfile(fileparts(mfilename('fullpath')), 'convdata_onesurface_rerun.csv');
end

zk     = 0;
nr     = 12;
nt     = 50;
ntheta = 200;
rmin   = 2.0;
rmaj   = 2.0;
jmag   = 1.0;
nvs    = [6 8 10 12 16];

% tolerance, orders run at it, and whether those curves are the plotted ones
runs = { 1e-6, [5 7 9], [1 0 0]
         1e-8, [7 9],   [1 1]   };

fid = fopen(outcsv, 'w');
fprintf(fid, 'n,p,nv,nu,npts,h,dof,tol,err,sigma_inf,iterD,iterW,time_s,on_plot\n');
fclose(fid);
fprintf('=== convoop_sweep start %s\n  -> %s\n', datestr(now), outcsv);

for k = 1:size(runs,1)
    tol = runs{k,1};
    for jn = 1:numel(runs{k,2})
        n = runs{k,2}(jn);
        on_plot = runs{k,3}(jn);
        for nv = nvs
            nu = nv*3;
            fprintf('\n\t========\n\tn = %d, nv = %d, tol = %g\n', n, nv, tol);

            t0 = tic;
            [dom, qnodes, qweights] = prepare_torus(n, nu, nv, nr, nt);
            dom = dom{1}; qnodes = qnodes{1}; qweights = qweights{1};
            domparams = [n, nu, nv];
            t_geom = toc(t0);

            t0 = tic;
            B0 = reftaylorsurffun(dom, n, ntheta, rmin, rmaj, jmag, zk);
            flux = 0;
            for i = 1:nr*nt
                B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes(:,i));
                flux = flux + B0eval(2)*qweights(i);
            end
            t_ref = toc(t0);

            ts = RefTaylorState({dom}, domparams, zk, flux, {B0}, ...
                {qnodes}, {qweights}, tol);

            % solve(true) prints the GMRES iteration counts but does not store
            % them, so capture stdout and read them back out.
            t0 = tic;
            solvelog = evalc('ts = ts.solve(true);');
            B = ts.surface_B();
            t_solve = toc(t0);
            fprintf('%s', solvelog);

            iterD = regexp(solvelog, 'A11\*D = A12.*?/\s*(\d+)\s*iter', 'tokens', 'once');
            iterW = regexp(solvelog, 'A11\*W = A12.*?/\s*(\d+)\s*iter', 'tokens', 'once');
            iterD = str2double([iterD{:}]);
            iterW = str2double([iterW{:}]);

            err = vecinfnorm(B0 - B{1})/vecinfnorm(B0);
            npts = nu*nv*n*n;
            h = 1/sqrt(nu*nv);

            fprintf(['n=%d nv=%d tol=%g  err=%.6e  iter=%d/%d\n' ...
                     '  geom %.1f s + ref %.1f s + solve %.1f s = %.1f s\n'], ...
                n, nv, tol, err, iterD, iterW, t_geom, t_ref, t_solve, ...
                t_geom+t_ref+t_solve);

            fid = fopen(outcsv, 'a');
            fprintf(fid, '%d,%d,%d,%d,%d,%.10f,%.4f,%g,%.6e,%.4e,%g,%g,%.0f,%d\n', ...
                n, n-1, nv, nu, npts, h, sqrt(npts), tol, err, ...
                norm(ts.sigma{1}, inf), iterD, iterW, t_solve, on_plot);
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
