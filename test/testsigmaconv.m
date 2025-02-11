%% self-convegence test for sigma

ns = 11;%[7 8 9 10];
nfine = 12;
nv = 7;
nu = nv*3;

tols = 1e-6;
zk = 0;

% load sigma and alpha
sigma_o_fine = load('../saved_data/sigma_o_n12_nv7.mat', 'sigma_o');
sigma_o_fine = sigma_o_fine.sigma_o;
sigma_i_fine = load('../saved_data/sigma_i_n12_nv7.mat', 'sigma_i');
sigma_i_fine = sigma_i_fine(1).sigma_i;

sigma_errs = zeros(4,size(ns,2));
sind = 1;

for n = ns

    fprintf('n = %d, nv = %d\n', n, nv)

    ao = 1.0; % outer minor radius
    ai = 0.6; % inner minor radius

    nr = 16;
    nt = 40;
    np = 40;
    [dom, qnodes, qweights] = prepare_torus(n,nu,nv,n,nu,nv,ao,ai,nr,nt,np);

    ntheta = 1e3;
    rmin = 2.0;
    rmaj = 2.0;
    jmag = 1.0;
    B0o = reftaylorsurffun(dom{1},n,nu,nv,ntheta,rmin,rmaj,jmag,zk);
    B0i = reftaylorsurffun(dom{2},n,nu,nv,ntheta,rmin,rmaj,jmag,zk);
        
    flux = zeros(1,2);
    for i = 1:nr*nt
        B0eval = reftaylor(ntheta,rmin,rmaj,jmag,zk,qnodes{1}(:,i));
        flux(1) = flux(1) + B0eval(2)*qweights{1}(i);
    end
    for i = 1:nr*np
        [B0eval, B0srcpts] = reftaylor(ntheta,rmin,rmaj,jmag,zk,qnodes{2}(:,i));
        flux(2) = flux(2) - B0eval(3)*qweights{2}(i);
    end
        
    B0 = {B0o,B0i};
    domparams = [n nu nv];
    ts = RefTaylorState(dom,domparams,zk,flux,B0,qnodes,qweights,tols);
    ts = ts.solve(true);

    sigma_o_coarse = ts.sigma{1};
    sigma_i_coarse = ts.sigma{2};
    sigma_o_coarse = resample(sigma_o_coarse,nfine);
    sigma_i_coarse = resample(sigma_i_coarse,nfine);

    sigma_errs(1,sind) = n;
    sigma_errs(2,sind) = nv;
    sigma_errs(3,sind) = nu;
    % sigma_errs(4,sind) = max(vecinfnorm(sigma_o_coarse-sigma_o_fine), ...
    %     vecinfnorm(sigma_i_coarse-sigma_i_fine))/ ...
    %     max(vecinfnorm(sigma_o_fine), vecinfnorm(sigma_i_fine));
    sigma_errs(4,sind) = max(norm(sigma_o_coarse-sigma_o_fine,'inf'), ...
        norm(sigma_i_coarse-sigma_i_fine,'inf'))/ ...
        max(norm(sigma_o_fine,'inf'), norm(sigma_i_fine,'inf'));

    sind = sind + 1;
end

loglog(sqrt(sigma_errs(1,1:size(nv,2)).^2.*sigma_errs(2,1:size(nv,2)).* ...
    sigma_errs(3,1:size(nv,2))), sigma_errs(4,1:size(nv,2)), 'o-')
hold on
for i = 1:size(ns,2)-1
    ind = i*size(nv,2)+1:(i+1)*size(nv,2);
    loglog(sqrt(sigma_errs(1,ind).^2.*sigma_errs(2,ind).* ...
        sigma_errs(3,ind)), sigma_errs(4,ind), 'o-')
end