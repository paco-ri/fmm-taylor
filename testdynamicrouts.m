nsph = 6;
nref = 3;
dom = surfacemesh.sphere(nsph,nref);
S = surfer.surfacemesh_to_surfer(dom);
[srcvals,~,~,~,~,wts] = extract_arrays(S);
vn = normal(dom);

ndeg = 3; 
mdeg = 0;

f = spherefun.sphharm(ndeg,mdeg);
sigma = surfacefun(@(x,y,z) f(x,y,z),dom);
sigvals = surfacefun_to_array(sigma,dom,S);

zk = 1.0+0.0i;
eps = 1e-7;

% =========================================

gradSksigma = taylor.dynamic.eval_gradSk(S,zk,sigvals.',eps);
ngradSksigma_ex = -(ndeg+1)/(2*ndeg+1).*sigvals;
gradSksigma = array_to_surfacefun(gradSksigma.',dom,S);
ngradSksigma = dot(vn,gradSksigma);
ngradSksigma = surfacefun_to_array(ngradSksigma,dom,S);
err = sqrt(sum(abs(ngradSksigma_ex - ngradSksigma).^2.*wts));
fprintf('error w/o Q = %f\n',err)