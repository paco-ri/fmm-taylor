%% self-convergence test for surface area

ns = [7 9 11];
nfine = 16;
nvs = [7 10 13];

% these don't matter for this test
nr = 16;
nt = 40;
np = 40;
ao = 1.0;
ai = 0.6;

err = zeros(4,size(ns,2)*size(nvs,2));
ind = 1;

nvfine = 16;
nufine = 3*nvfine;
torus_fine = prepare_torus(nfine,nufine,nvfine,nfine,nufine,nvfine,...
    ao,ai,nr,nt,np);
area_fine = surfacearea(torus_fine{1});

for nv = nvs
    nu = nv*3;
    for n = ns
        fprintf('n = %d, nv = %d\n', n, nv)
        disp('---------------')
        nu = nv*3;
        dom = prepare_torus(n,nu,nv,n,nu,nv,ao,ai,nr,nt,np);
        area_coarse = surfacearea(dom{1});
        
        err(1,ind) = n; % number of points on each patch
        err(2,ind) = nv; % geom. param. 
        err(3,ind) = nu; % geom. param. 
        err(4,ind) = abs(area_fine-area_coarse);
        ind = ind + 1;
    end
end

figure(1)
loglog(sqrt(err(1,1:size(nvs,2)).^2.*err(2,1:size(nvs,2)).* ...
    err(3,1:size(nvs,2))), err(4,1:size(nvs,2)), 'o-')
hold on
for i = 1:size(ns,2)-1
    ind = i*size(nvs,2)+1:(i+1)*size(nvs,2);
    loglog(sqrt(err(1,ind).^2.*err(2,ind).* ...
        err(3,ind)), err(4,ind), 'o-')
end