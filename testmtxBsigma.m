% small problem
n = 4;
nu = 2;
nv = nu*3;

dom = surfacemesh.torus(n, nu, nv);
[S,domsurfer] = surfer.surfacemesh_to_surfer(dom);

x0 = 0.2; y0 = 0.2; z0 = 0.2;
denom = @(x,y,z) ((x-x0).^2 + (y-y0).^2 + (z-z0).^2).^(3/2);
sigma = surfacefun(denom,dom);

eps = 1e-8;
opts_quad = [];
opts_quad.format='rsc';
Q = lap3d.dirichlet.get_quadrature_correction(S,[1.0 0],eps,S,opts_quad);

Bsigma = mtxBsigma(S,dom,sigma,eps,Q);

