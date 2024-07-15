nsph = 5;
nref = 4;
dom = surfacemesh.sphere(nsph,nref);
S = surfer.surfacemesh_to_surfer(dom);
[srcvals,~,~,~,~,wts] = extract_arrays(S);
vn = normal(dom);

ndeg = 4; 
mdeg = 1;

f = spherefun.sphharm(ndeg,mdeg);
sigma = surfacefun(@(x,y,z) f(x,y,z),dom);
sigvals = surfacefun_to_array(sigma,dom,S);

unm = grad(sigma)./sqrt(ndeg*(ndeg+1));
xnm = cross(vn,unm);
cunm = 1.0;
cxnm = 0;
rjvec = cunm.*unm + cxnm.*xnm;
rjvals = surfacefun_to_array(rjvec,dom,S);

zk = 1.0+0.0i;
eps = 1e-9;

% =========================================

if false
fprintf('nsph = %d, nref = %d\n',nsph,nref)
gradSksigma = taylor.dynamic.eval_gradSk(S,zk,sigvals.',eps);
gradSksigma = array_to_surfacefun(gradSksigma.',dom,S);
ngradSksigma = dot(vn,gradSksigma);
ngradSksigma = surfacefun_to_array(ngradSksigma,dom,S);

% see eqn. 29 in Vico et al. 
ngradSksigma_ex = 1i * zk^2 * sphj(ndeg,zk) * sphhp(ndeg,zk) .* sigvals;
ngradSksigma_ex = ngradSksigma_ex + sigvals./2;

err = sqrt(sum(abs(ngradSksigma_ex - ngradSksigma).^2.*wts));
fprintf('error w/o Q = %f\n',err)
end

% =========================================

curlSksigma = taylor.dynamic.eval_curlSk(S,zk,rjvals.',eps);
curlSksigma = array_to_surfacefun(curlSksigma.',dom,S);
nxcurlSksigma = cross(vn,curlSksigma);
% nxcurlSksigma = surfacefun_to_array(nxcurlSksigma,dom,S);

% see eqn. 54 in Vico et al.
nxcurlSksigma_ex = cunm * 1i * ricjp(ndeg,zk) * rich(ndeg,zk) .* unm;
nxcurlSksigma_ex = nxcurlSksigma_ex - cxnm * 1i * ricj(ndeg,zk) * richp(ndeg,zk) .* xnm;
nxcurlSksigma_ex = nxcurlSksigma_ex - rjvec./2;

% err = sqrt(sum(norm(ngradSksigma_ex - ngradSksigma,2).^2.*wts));
err = norm(norm(nxcurlSksigma - nxcurlSksigma_ex));
fprintf('error w/o Q = %f\n',err)

% =========================================

function J = sphj(n, z)
J = sqrt(pi/2./z) .* besselj(n+1/2,z);
end

function Jp = sphjp(n, z)
if n == 0
    Jp = -sphj(1, z);
else
    Jp = sphj(n-1,z) - (n+1)./z.*sphj(n,z);
end
end

function Y = sphy(n, z)
Y = sqrt(pi/2./z) .* bessely(n+1/2,z);
end

function Yp = sphyp(n, z)
if n == 0
    Yp = -sphy(1, z);
else
    Yp = sphy(n-1,z) - (n+1)./z.*sphy(n,z);
end
end

function H = sphh(n, z)
H = sphj(n,z) + 1i.*sphy(n,z);
end

function Hp = sphhp(n, z)
Hp = sphjp(n,z) + 1i.*sphyp(n,z);
end

function rJ = ricj(n, z)
rJ = sqrt(pi.*z./2) .* besselj(n+1/2,z);
end

function rJp = ricjp(n, z)
rJp = z./(2*n+1) .* (n.*sphj(n-1,z) - (n+1).*sphj(n+1,z));
end

function rH = rich(n, z)
rH = sqrt(pi.*z./2) .* (besselj(n+1/2,z) + 1i*bessely(n+1/2,z));
end

function rHp = richp(n, z)
rHp = z./(2*n+1) .* (n.*sphh(n-1,z) - (n+1).*sphh(n+1,z));
end