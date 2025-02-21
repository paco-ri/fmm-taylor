%% self-convergence test for sigma, TaylorState objects

ns = [7 8 9 10];
nfine = 12;
nv = 10;
nu = nv*3;

nr = 16;
nt = 40;
np = 40;

tols = 1e-6;
zk = 0;

ao = 1.0; % outer minor radius
ai = 0.6; % inner minor radius
flux = [1.0 0.7];

% do_load = false;
% % load/compute sigma
% if do_load
%     % Add code to load sigma here
% else
%     [dom, qnodes, qweights] = prepare_torus(nfine,nu,nv,nfine,nu,nv,ao,ai,nr,nt,np);
% 
%     domparams = [nfine nu nv];
%     ts = TaylorState(dom,domparams,zk,flux,tols);
%     ts = ts.solve(true);
% 
%     sigma_o_fine = ts.sigma{1};
%     sigma_i_fine = ts.sigma{2};
% 
%     save('../saved_data/sigma_o_fine.mat', 'sigma_o_fine');
%     save('../saved_data/sigma_i_fine.mat', 'sigma_i_fine');
% end

sigma_errs = zeros(4,size(ns,2));
sind = 1;

for n = ns

    fprintf('n = %d, nv = %d\n', n, nv)

    [dom, qnodes, qweights] = prepare_torus(n,nu,nv,n,nu,nv,ao,ai,nr,nt,np);

    domparams = [n nu nv];
    ts = TaylorState(dom,domparams,zk,flux,tols);
    ts = ts.solve(true);
    sigma_o_coarse = resample(ts.sigma{1},nfine);
    sigma_i_coarse = resample(ts.sigma{2},nfine);

    sigma_errs(1,sind) = n;
    sigma_errs(2,sind) = nv;
    sigma_errs(3,sind) = nu;
    sigma_errs(4,sind) = max(norm(sigma_o_coarse-sigma_o_fine,'inf'), ...
        norm(sigma_i_coarse-sigma_i_fine,'inf'))/ ...
        max(norm(sigma_o_fine,'inf'), norm(sigma_i_fine,'inf'));

    sind = sind + 1;
    tols = tols*1e-1;
end

loglog(sqrt(sigma_errs(1,1:size(nv,2)).^2.*sigma_errs(2,1:size(nv,2)).* ...
    sigma_errs(3,1:size(nv,2))), sigma_errs(4,1:size(nv,2)), 'o-')
hold on
for i = 1:size(ns,2)-1
    ind = i*size(nv,2)+1:(i+1)*size(nv,2);
    loglog(sqrt(sigma_errs(1,ind).^2.*sigma_errs(2,ind).* ...
        sigma_errs(3,ind)), sigma_errs(4,ind), 'o-')
end