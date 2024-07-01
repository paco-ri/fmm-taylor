n = 12;
nu = 6;
nv = nu*3;

dom = circulartorus(n, nu, nv);

pdo = [];
pdo.lap = 1;
L = surfaceop(dom, pdo);
L.rankdef = true;

tic
[V, D] = eigs(L);
toc

% diag(D)
% plot(V)

S = surfer.surfacemesh_to_surfer(dom);
eps = 1e-10;
Q = taylor.static.get_quadrature_correction(S,eps,S);
sigma = V(3);
sigvals = surfacefun_to_array(sigma,dom,S);
S0sigvals = taylor.static.eval_gradS0(S,sigvals,eps,S,Q);
S0sig = array_to_surfacefun(S0sigvals.',dom,S);

vn = normal(dom);
nS0sig = dot(vn,S0sig);

figure(1)
plot(sigma)
colorbar

figure(2)
plot(nS0sig)
colorbar

nS0sigvals = surfacefun_to_array(nS0sig,dom,S);
figure(3)
plot(real(sigvals./nS0sigvals), 'o')
hold on
plot(imag(sigvals./nS0sigvals), 'o')