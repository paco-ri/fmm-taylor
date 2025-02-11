%% self-convergence test for mH

ns = [7 9 11];
nfine = 16;
nvs = [7 10 13];

% these don't matter for this test
nr = 16;
nt = 40;
np = 40;
ao = 1.0;
ai = 0.6;

% high-order mH
% dom = prepare_torus(16,30,10,16,30,10,ao,ai,nr,nt,np);
% domparams = [16 30 10];
% domain = Domain(dom,domparams);
% mHfine = domain.mH;

mHerr = zeros(4,size(ns,2)*size(nvs,2));
mHind = 1;

for nv = nvs
    % high-order mH
    nu = nv*3;
    % domfine = prepare_torus(nfine,nu,nv,nfine,nu,nv,ao,ai,nr,nt,np);
    % domfineparams = [nfine nu nv];
    % domainfine = Domain(domfine,domfineparams,nfine);
    % mHfine = domainfine.mH;
    for n = ns
        fprintf('n = %d, nv = %d\n', n, nv)
        disp('---------------')
        nu = nv*3;
        dom = prepare_torus(n,nu,nv,n,nu,nv,ao,ai,nr,nt,np);
        domparams = [n nu nv];
        domain = Domain(dom,domparams,nfine);
        % mHcoarse = cell(1,2);
        % for i = 1:2
        %     mHcoarse{i} = surfacefunv(domfine{i});
        %     for j = 1:3
        %         mHcoarse{i}.components{j} = resample(domain.mH{i}.components{j},nfine);
        %     end
        % end
        mHfine = domain.mH;
        mHcoarse = domain.mHcoarse;
        mHerr(1,mHind) = n; % number of points on each patch
        mHerr(2,mHind) = nv; % geom. param. 
        mHerr(3,mHind) = nu; % geom. param. 
        mHerr(4,mHind) = max(vecinfnorm(mHcoarse{1}-mHfine{1}),vecinfnorm(mHcoarse{2}-mHfine{2}))...
            /max(vecinfnorm(mHfine{1}),vecinfnorm(mHfine{2})); 
        mHind = mHind + 1;
    end
end

figure(1)
loglog(sqrt(mHerr(1,1:size(nvs,2)).^2.*mHerr(2,1:size(nvs,2)).* ...
    mHerr(3,1:size(nvs,2))), mHerr(4,1:size(nvs,2)), 'o-')
hold on
for i = 1:size(ns,2)-1
    ind = i*size(nvs,2)+1:(i+1)*size(nvs,2);
    loglog(sqrt(mHerr(1,ind).^2.*mHerr(2,ind).* ...
        mHerr(3,ind)), mHerr(4,ind), 'o-')
end

