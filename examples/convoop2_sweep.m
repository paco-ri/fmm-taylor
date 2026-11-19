function convoop2_sweep(outcsv, nthreads, runs)
%CONVOOP2_SWEEP Two-surface convergence sweep (plot_conv_twosurface.m).
%
%   Setup: prepare_torus(n,nu,nv,n,nu,nv,ao,ai,nr,nt,np), zk = 0.5,
%   nu = 3*nv, ao = 1.0, ai = 0.6, nt = np = 40, ntheta = 200, scalar tol.
%
%   Each row of RUNS is {tol, n, nvs, nr, on_plot}. The default runs the
%   plotted curves, p = 12 last so it can be cut short:
%     p = 6, 8:    tol = 1e-6,  nr = 16  (the published settings; p = 8,
%                  nv = 8 should give a surface error of ~3.848e-05, 43/28)
%     p = 10, 12:  tol = 1e-10, nr = 24  (nr = 16's cross-section rule is
%                  only accurate to ~2e-8)
%
%   TWO ERRORS, both relative L^inf of B - B0, max over components:
%     err      on both surfaces' nodes (surface_B), with err_s1/err_s2 per
%              surface. Dominated at fine meshes by surface_B's reconstruction
%              error at patch edges, which overstates the error of the
%              solution (see diag_edge.m).
%     err_int  at 768 points in the shell: u, v = 0:22.5:337.5 deg on the
%              tori of scale s = 0.7, 0.8, 0.9 (torus_points.m), normalized
%              by the largest |B0| among those points. The outer layers are
%              at scale distance 0.1 from each wall.
%
%   nthreads (optional, [] to skip) sets maxNumCompThreads. On jupiter the
%   MEX's OpenMP runs on MATLAB's runtime, so OMP_NUM_THREADS has no effect.
%
%   Rows are appended as each configuration finishes. A MATLAB error in one
%   configuration writes a NaN row and the sweep continues; a self-quadrature
%   failure stops the process (ier check in fmm3dbie's ggq-quads.f).

here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(outcsv)
    outcsv = fullfile(here, 'convdata_twosurface_interior.csv');
end
if nargin >= 2 && ~isempty(nthreads)
    maxNumCompThreads(nthreads);
end
if nargin < 3 || isempty(runs)
    runs = { 1e-6,  7,  5:8, 16, 1
             1e-6,  9,  5:8, 16, 1
             1e-10, 11, 5:8, 24, 1
             1e-10, 13, 5:8, 24, 1 };
end

zk     = 0.5;
nt     = 40;
np     = 40;
ntheta = 200;
ao     = 1.0;
ai     = 0.6;
rmin   = 2.0;
rmaj   = 2.0;
jmag   = 1.0;

[uu, vv, ss] = ndgrid((0:22.5:337.5)*pi/180, (0:22.5:337.5)*pi/180, [0.7 0.8 0.9]);
intpts = torus_points(uu, vv, ss);
B0int = zeros(3, size(intpts,2));
for i = 1:size(intpts,2)
    B0int(:,i) = reftaylor(ntheta, rmin, rmaj, jmag, zk, intpts(:,i));
end
B0intmax = max(abs(B0int), [], 'all');

stagepats = {'get Laplace/Helmholtz quad\. corr\.'
             'get \+taylor routine quad\. corr\.'
             'get XS Laplace/Helmholtz quad\. corr\.'
             'get XS \+taylor routine quad\. corr\.'
             'compute sigma and alpha'
             'compute m'};

fid = fopen(outcsv, 'w');
fprintf(fid, ['n,p,nv,nu,npts,h,dof,tol,nr,err,err_s1,err_s2,sigma_inf,' ...
    'iterD,iterW,time_s,on_plot,t_quad_lh,t_quad_tay,t_quad_xs_lh,' ...
    't_quad_xs_tay,t_gmres_d,t_gmres_w,t_sigma_alpha,t_m,t_solve_sum,' ...
    't_geom,t_b0,t_flux,t_surfb,err_int,t_int\n']);
fclose(fid);
fprintf('=== convoop2_sweep start %s\n  -> %s\n', datestr(now), outcsv);
fprintf('  maxNumCompThreads = %d, OMP_NUM_THREADS = %s\n', ...
    maxNumCompThreads, getenv('OMP_NUM_THREADS'));

for k = 1:size(runs,1)
    [tol, n, nvs, nr, on_plot] = runs{k,:};
    for nv = nvs
        nu = 3*nv;
        npts = 2*nu*nv*n*n;
        fprintf('\n\t========\n\tn = %d, nv = %d, tol = %g, nr = %d\n', ...
            n, nv, tol, nr);
        try
            t0 = tic;
            [dom, qnodes, qweights] = prepare_torus(n, nu, nv, n, nu, nv, ...
                ao, ai, nr, nt, np);
            t_geom = toc(t0);

            t0 = tic;
            B0 = {reftaylorsurffun(dom{1}, n, ntheta, rmin, rmaj, jmag, zk), ...
                  reftaylorsurffun(dom{2}, n, ntheta, rmin, rmaj, jmag, zk)};
            t_b0 = toc(t0);

            t0 = tic;
            flux = zeros(1,2);
            for i = 1:nr*nt
                B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes{1}(:,i));
                flux(1) = flux(1) + B0eval(2)*qweights{1}(i);
            end
            for i = 1:nr*np
                B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes{2}(:,i));
                flux(2) = flux(2) - B0eval(3)*qweights{2}(i);
            end
            t_flux = toc(t0);

            ts = RefTaylorState(dom, [n, nu, nv], zk, flux, B0, qnodes, ...
                qweights, tol);

            % solve(true) prints but does not return its stage timings
            t0 = tic;
            solvelog = evalc('ts = ts.solve(true);');
            t_solve = toc(t0);
            t0 = tic;
            B = ts.surface_B();
            t_surfb = toc(t0);
            fprintf('%s', solvelog);

            tstage = nan(1,6);
            for ip = 1:6
                tok = regexp(solvelog, ['^' stagepats{ip} ': *([\d.eE+-]+) s'], ...
                    'tokens', 'once', 'lineanchors');
                if ~isempty(tok); tstage(ip) = str2double(tok{1}); end
            end
            gd = regexp(solvelog, 'A11\*D = A12[^:]*: *([\d.eE+-]+) s / *(\d+) iter', ...
                'tokens', 'once');
            gw = regexp(solvelog, 'A11\*W = A12[^:]*: *([\d.eE+-]+) s / *(\d+) iter', ...
                'tokens', 'once');
            if isempty(gd); gd = {'NaN','NaN'}; end
            if isempty(gw); gw = {'NaN','NaN'}; end

            B0max = max(vecinfnorm(B0{1}), vecinfnorm(B0{2}));
            err_s = [vecinfnorm(B0{1}-B{1}), vecinfnorm(B0{2}-B{2})]/B0max;
            sigma_inf = max(norm(ts.sigma{1}, inf), norm(ts.sigma{2}, inf));

            t0 = tic;
            Bint = ts.interior_B(intpts);
            t_int = toc(t0);
            err_int = max(abs(Bint - B0int), [], 'all')/B0intmax;
        catch ME
            fprintf('FAILED n=%d nv=%d: %s\n', n, nv, ME.message);
            [t_geom, t_b0, t_flux, t_solve, t_surfb, sigma_inf, t_int, ...
                err_int] = deal(NaN);
            tstage = nan(1,6); gd = {'NaN','NaN'}; gw = {'NaN','NaN'};
            err_s = [NaN NaN];
        end
        t_gmres_d = str2double(gd{1}); iterD = str2double(gd{2});
        t_gmres_w = str2double(gw{1}); iterW = str2double(gw{2});
        err = max(err_s);

        fprintf(['n=%d nv=%d tol=%g nr=%d  err=%.6e (%.3e, %.3e)  ' ...
                 'err_int=%.6e  iter=%d/%d\n' ...
                 '  geom %.1f + B0 %.1f + flux %.1f + solve %.1f ' ...
                 '+ surfB %.1f + int %.1f s\n'], n, nv, tol, nr, err, err_s, ...
            err_int, iterD, iterW, t_geom, t_b0, t_flux, t_solve, t_surfb, t_int);

        fid = fopen(outcsv, 'a');
        fprintf(fid, ['%d,%d,%d,%d,%d,%.10f,%.4f,%g,%d,%.6e,%.6e,%.6e,%.4e,' ...
            '%g,%g,%.0f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,' ...
            '%.4f,%.4f,%.4f,%.4f,%.6e,%.4f\n'], ...
            n, n-1, nv, nu, npts, 1/sqrt(nu*nv), sqrt(npts), tol, nr, err, ...
            err_s, sigma_inf, iterD, iterW, t_solve+t_surfb, on_plot, ...
            tstage(1:4), t_gmres_d, t_gmres_w, tstage(5:6), sum(tstage), ...
            t_geom, t_b0, t_flux, t_surfb, err_int, t_int);
        fclose(fid);
    end
end

fprintf('\nDONE_SENTINEL %s\n', datestr(now));
end

function N = vecinfnorm(f)
N = max([norm(f.components{1}, inf), ...
         norm(f.components{2}, inf), ...
         norm(f.components{3}, inf)]);
end
