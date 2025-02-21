%% self-convergence test for curl S0[mH] on torus

ns = [7 8 9 10];
nfine = 13;
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

do_load = false;
fname = '../saved_data/curlS0mH_arr_crit0/curlS0mH_arr_nfine13_nv10.mat';
% load/compute Balpha
if do_load
    domain_fine = Domain(dom,domparams,nfine);
    curl_fine = load(fname, 'curl_fine');
    curl_fine = curl_fine.curl_fine;
else
    tol_fine = 1e-12;
    ts = TaylorState(dom,domparams,zk,flux,tols);
    disp('created fine TS')
    ts = ts.get_quad_corr_laphelm();
    disp('created fine TS Laplace/Helmholtz quad. corr.')
    ts = ts.get_quad_corr_taylor();
    disp('created fine TS grad/curl Sk quad. corr.')
    curl_fine = zeros(3,4*ts.domain.nptspersurf);
    inds = 1:ts.domain.nptspersurf;
    for i = 1:2
        for j = 1:2
            mHvals = surfacefun_to_array(ts.domain.mH{i},...
                ts.domain.dom{i},ts.domain.surf{i});
            curlSmH = taylor.static.eval_curlS0(ts.domain.surf{i},...
                mHvals.',tol_fine,ts.domain.surf{j},...
                ts.quad_opts_taylor{i,j});
            curl_fine(:,inds) = curlSmH;
            inds = inds + ts.domain.nptspersurf;
        end
    end
    disp('computed curl_fine')
    domain_fine = ts.domain;

    save(fname, 'curl_fine');
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
    curl_coarse = zeros(3,4*nptspersurf_fine);
    inds = 1:nptspersurf_fine;
    for i = 1:2
        for j = 1:2
            mHvals = surfacefun_to_array(ts.domain.mH{i},...
                ts.domain.dom{i},ts.domain.surf{i});
            curlSmH = taylor.static.eval_curlS0(ts.domain.surf{i},...
                mHvals.',tol_fine,ts.domain.surf{j},...
                ts.quad_opts_taylor{i,j});
            curlSmH = array_to_surfacefun(curlSmH.',ts.domain.dom{j},...
                ts.domain.surf{j});
            curlSmHcoarse = surfacefunv(domain_fine.dom{j});
            for k = 1:3
                curlSmHcoarse.components{k} = resample(...
                    curlSmH.components{k},nfine);
            end
            curlSmHcoarse = surfacefun_to_array(curlSmHcoarse,...
                domain_fine.dom{j},domain_fine.surf{j});
            curl_coarse(:,inds) = curlSmHcoarse.';
            inds = inds + nptspersurf_fine;
        end
    end

    errs(1,ind) = n;
    errs(2,ind) = nv;
    errs(3,ind) = nu;
    errs(4,ind) = max(abs(curl_coarse-curl_fine),[],'all');

    ind = ind + 1;
    % tols = tols*1e-1;
end

loglog(sqrt(errs(1,1:size(nv,2)).^2.*errs(2,1:size(nv,2)).* ...
    errs(3,1:size(nv,2))), abs(errs(4,1:size(nv,2))), 'o-')
hold on
for i = 1:size(ns,2)-1
    ind = i*size(nv,2)+1:(i+1)*size(nv,2);
    loglog(sqrt(errs(1,ind).^2.*errs(2,ind).* ...
        errs(3,ind)), abs(errs(4,ind)), 'o-')
end