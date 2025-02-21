%% self-convergence test for Balpha, TaylorState objects

ns = [7 8 9 10];
nfine = 12;
nv = 10;
nu = nv*3;

nr = 16;
nt = 40;
np = 40;

tols = 1e-9;
zk = 0;

ao = 1.0; % outer minor radius
ai = 0.6; % inner minor radius
flux = [1.0 0.7];
nptspersurf_fine = nfine^2*nu*nv;

% fine solutions
dom = prepare_torus(nfine,nu,nv,nfine,nu,nv,ao,ai,nr,nt,np);
domparams = [nfine nu nv];

do_load = true;
fname = '../saved_data/Balpha_arr_crit0/Balpha_arr_nfine12_nv10.mat';
% load/compute Balpha
if do_load
    domain_fine = Domain(dom,domparams,nfine);
    b_fine = load(fname, 'b_fine');
    b_fine = b_fine.b_fine;
else
    tol_fine = 1e-12;
    ts = TaylorState(dom,domparams,zk,flux,tols);
    disp('created fine TS')
    ts = ts.get_quad_corr_laphelm();
    disp('created fine TS Laplace/Helmholtz quad. corr.')
    ts = ts.get_quad_corr_taylor();
    disp('created fine TS grad/curl Sk quad. corr.')
    Balpha_fine = TaylorState.mtxBalpha(ts.domain,ts.zk, ...
        ts.eps_taylor,ts.eps_laphelm,ts.domain.surf, ...
        ts.quad_opts_taylor,ts.quad_opts_laphelm);
    disp('created Balpha_fine')
    b_fine = zeros(2*ts.domain.nptspersurf,2);
    for i = 1:2
        for j = 1:2
            inds = ts.domain.nptspersurf*(i-1)+1:ts.domain.nptspersurf*i;
            b_fine(inds,j) = surfacefun_to_array(Balpha_fine{i,j},...
                ts.domain.dom{i},ts.domain.surf{i});
        end
    end
    domain_fine = ts.domain;

    save(fname, 'b_fine');
end

errs = zeros(4,size(ns,2));
ind = 1;

disp('critical mode set to 0')

for n = ns

    fprintf('n = %d, nv = %d\n', n, nv)

    dom = prepare_torus(n,nu,nv,n,nu,nv,ao,ai,nr,nt,np);

    domparams = [n nu nv];
    ts = TaylorState(dom,domparams,zk,flux,tols);
    disp('created TS')
    ts = ts.get_quad_corr_laphelm();
    disp('created Lap/Helm quad corr')
    ts = ts.get_quad_corr_taylor();
    disp('created grad/curl Sk quad corr')
    Balpha_coarse = TaylorState.mtxBalpha(ts.domain,ts.zk, ...
        ts.eps_taylor,ts.eps_laphelm,ts.domain.surf, ...
        ts.quad_opts_taylor,ts.quad_opts_laphelm);
    disp('created Balpha_coarse')
    for i = 1:2
        for j = 1:2
            Balpha_coarse{i,j} = resample(Balpha_coarse{i,j},nfine);
        end
    end
    b_coarse = zeros(2*nptspersurf_fine,2);
    for i = 1:2
        for j = 1:2
            inds = nptspersurf_fine*(i-1)+1:nptspersurf_fine*i;
            b_coarse(inds,j) = surfacefun_to_array(Balpha_coarse{i,j},...
                domain_fine.dom{i},domain_fine.surf{i});
        end
    end

    errs(1,ind) = n;
    errs(2,ind) = nv;
    errs(3,ind) = nu;
    errs(4,ind) = max(abs(b_coarse-b_fine),[],'all')/max(abs(b_fine),[],'all');

    ind = ind + 1;
    tols = tols*1e-1;
end

loglog(sqrt(errs(1,1:size(nv,2)).^2.*errs(2,1:size(nv,2)).* ...
    errs(3,1:size(nv,2))), abs(errs(4,1:size(nv,2))), 'o-')
hold on
for i = 1:size(ns,2)-1
    ind = i*size(nv,2)+1:(i+1)*size(nv,2);
    loglog(sqrt(errs(1,ind).^2.*errs(2,ind).* ...
        errs(3,ind)), abs(errs(4,ind)), 'o-')
end