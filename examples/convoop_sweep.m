function convoop_sweep(outcsv, runs, nvs)
%CONVOOP_SWEEP One-surface convergence sweep (plot_conv_onesurface.m).
%
%   CONVOOP_SWEEP(OUTCSV, RUNS, NVS). RUNS rows are {tol, orders n, on_plot};
%   the default is p = 6, 8 at tol = 1e-8 and p = 10, 12 at tol = 1e-10, the
%   plotted curves. NVS defaults to 6, 8, 10, 12, 16.
%
%   Setup: prepare_torus, zk = 0, nu = 3*nv, nr = 12, nt = 50, ntheta = 200,
%   rmin = rmaj = 2, jmag = 1.
%
%   TWO ERRORS, both relative L^inf of B - B0, max over components:
%     err      on the surface nodes (surface_B). Dominated at fine meshes by
%              surface_B's reconstruction error at patch edges, which
%              overstates the error of the solution (see diag_edge.m).
%     err_int  at 768 interior points: u, v = 0:22.5:337.5 deg on the tori of
%              scale s = 0.3, 0.6, 0.9 (torus_points.m), normalized by the
%              largest |B0| among those points. s = 0.9 is the nearest layer
%              to the wall, at scale distance 0.1, the same as the nearest
%              layers in convoop2_sweep.m.
%
%   Rows are appended as they finish, so the run can be interrupted without
%   losing completed work.
%
%   TIMING COLUMNS. time_s keeps its original meaning (solve + surface_B); the
%   original CSV's definition was not recorded, so the breakdown is carried
%   alongside it. The six solve stages are chained inside RefTaylorState.solve,
%   so t_solve_sum should match the wall time of the solve -- a NaN there means
%   a stage line was missed in the log.
%     t_quad_lh, t_quad_tay          surface quadrature corrections
%     t_quad_xs_lh, t_quad_xs_tay    cross-section ditto (small, ~nr*nt points)
%     t_gmres_d, t_gmres_w           the two GMRES solves, pairing with iterD/iterW
%     t_sigma_alpha                  enclosing stage; minus the two GMRES times
%                                    it gives the direct (A21*D + A22)^-1 cost
%     t_m                            compute m (mH only at zk = 0)
%     t_solve_sum                    sum of the six stages
%     t_geom, t_b0, t_flux, t_surfb  outside the solve
%     t_int                          interior_B at the 768 interior points

if nargin < 1 || isempty(outcsv)
    outcsv = fullfile(fileparts(mfilename('fullpath')), 'convdata_onesurface_interior.csv');
end
if nargin < 2 || isempty(runs)
    runs = { 1e-8,  [7 9],   [1 1]
             1e-10, [11 13], [1 1] };
end
if nargin < 3 || isempty(nvs)
    nvs = [6 8 10 12 16];
end

zk     = 0;
nr     = 12;
nt     = 50;
ntheta = 200;
rmin   = 2.0;
rmaj   = 2.0;
jmag   = 1.0;

[uu, vv, ss] = ndgrid((0:22.5:337.5)*pi/180, (0:22.5:337.5)*pi/180, [0.3 0.6 0.9]);
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
fprintf(fid, ['n,p,nv,nu,npts,h,dof,tol,err,sigma_inf,iterD,iterW,time_s,on_plot,' ...
    't_quad_lh,t_quad_tay,t_quad_xs_lh,t_quad_xs_tay,t_gmres_d,t_gmres_w,' ...
    't_sigma_alpha,t_m,t_solve_sum,t_geom,t_b0,t_flux,t_surfb,err_int,t_int\n']);
fclose(fid);
fprintf('=== convoop_sweep start %s\n  -> %s\n', datestr(now), outcsv);
fprintf('  maxNumCompThreads = %d\n', maxNumCompThreads);

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
            t_b0 = toc(t0);

            t0 = tic;
            flux = 0;
            for i = 1:nr*nt
                B0eval = reftaylor(ntheta, rmin, rmaj, jmag, zk, qnodes(:,i));
                flux = flux + B0eval(2)*qweights(i);
            end
            t_flux = toc(t0);

            ts = RefTaylorState({dom}, domparams, zk, flux, {B0}, ...
                {qnodes}, {qweights}, tol);

            % solve(true) prints its stage timings and GMRES iteration counts
            % but does not store them, so capture stdout and read them back out.
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
            t_gmres_d = str2double(gd{1}); iterD = str2double(gd{2});
            t_gmres_w = str2double(gw{1}); iterW = str2double(gw{2});

            err = vecinfnorm(B0 - B{1})/vecinfnorm(B0);
            npts = nu*nv*n*n;
            h = 1/sqrt(nu*nv);

            t0 = tic;
            Bint = ts.interior_B(intpts);
            t_int = toc(t0);
            err_int = max(abs(Bint - B0int), [], 'all')/B0intmax;

            fprintf(['n=%d nv=%d tol=%g  err=%.6e  err_int=%.6e  iter=%d/%d\n' ...
                     '  geom %.1f + B0 %.1f + flux %.1f + solve %.1f ' ...
                     '+ surfB %.1f + int %.1f s  (stages sum %.1f s)\n'], ...
                n, nv, tol, err, err_int, iterD, iterW, t_geom, t_b0, t_flux, ...
                t_solve, t_surfb, t_int, sum(tstage));

            fid = fopen(outcsv, 'a');
            fprintf(fid, ['%d,%d,%d,%d,%d,%.10f,%.4f,%g,%.6e,%.4e,%g,%g,%.0f,%d,' ...
                '%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,' ...
                '%.6e,%.4f\n'], ...
                n, n-1, nv, nu, npts, h, sqrt(npts), tol, err, ...
                norm(ts.sigma{1}, inf), iterD, iterW, t_solve+t_surfb, on_plot, ...
                tstage(1), tstage(2), tstage(3), tstage(4), t_gmres_d, t_gmres_w, ...
                tstage(5), tstage(6), sum(tstage), t_geom, t_b0, t_flux, t_surfb, ...
                err_int, t_int);
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
