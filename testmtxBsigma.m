% ===== copied from lap3d matlab test ==============
nsph = 5;
nref = 3;
dom = surfacemesh.sphere(nsph,nref);
S = surfer.surfacemesh_to_surfer(dom);
[srcvals,~,~,~,~,wts] = extract_arrays(S);

ndeg = 1;
mdeg = 0;

rr = sqrt(S.r(1,:).^2 + S.r(2,:).^2 + S.r(3,:).^2);
rhs = S.r(3,:)./rr;
rhs = rhs(:);

sigma = array_to_surfacefun(rhs,dom,S);
eps = 1e-7;

% ===================================================

% without precomputed quadrature
Bsigma = mtxBsigma(S,dom,sigma,eps);
Bsigma_ex = (ndeg+1)/(2*ndeg+1).*rhs - rhs;
Bsigma = surfacefun_to_array(Bsigma,dom,S);
err = sqrt(sum((Bsigma_ex - Bsigma).^2.*wts));
fprintf('error w/o Q = %f\n',err)

% with precomputed quadrature
Q = taylor.static.get_quadrature_correction(S,eps,S);
Bsigma = mtxBsigma(S,dom,sigma,eps,S,Q);
Bsigma = surfacefun_to_array(Bsigma,dom,S);
err = sqrt(sum((Bsigma_ex - Bsigma).^2.*wts));
fprintf('error w/ Q = %f\n',err)


