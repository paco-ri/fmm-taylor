% ===== copied from lap3d matlab test ==============
% For dan meshes

nsph = 5;
nref = 3;
dom = surfacemesh.sphere(nsph,nref);
S = surfer.surfacemesh_to_surfer(dom);

[srcvals,~,~,~,~,wts] = extract_arrays(S);
rr2 = sum(wts);
rr3 = norm(S.r-S.n);

% tic, plot(S); toc;

ndeg = 1;
mdeg = 0;

r2 = S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2;
rho2 = S.r(1,:).^2 + S.r(2,:).^2;
npts = size(r2,2);
rhs = zeros([3 npts]);
rhs(1,:) = S.r(1,:).*S.r(3,:)./r2;
rhs(2,:) = S.r(2,:).*S.r(3,:)./r2;
rhs(3,:) = -rho2./r2;
% rhs = surface grad cos(theta) = (sin(theta)/r) theta^

eps = 1e-7;

% ===================================================

rjvec = array_to_surfacefun(rhs.',dom,S);
ncurlS0j = mtxBalpha(S,dom,rjvec,eps);
vn = normal(dom);
ncurlS0j = surfacefun_to_array(ncurlS0j,dom,S);
ncurlS0j = ncurlS0j.';
err = norm(sqrt(sum(ncurlS0j.^2.*wts.')));
fprintf('sphere(%d,%d)\n',nsph,nref)
fprintf('error = %f\n',err)
