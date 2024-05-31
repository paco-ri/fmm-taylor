S = surfer.sphere(6,1,2,11);

rhs = ones([S.npts 1]);
dpars = [1.0, 0.0];
zpars = complex([1.0, 1.0, 0.0]);
opts_quad = [];
opts_quad.format='rsc';

Q = lap3d.dirichlet.get_quadrature_correction(S,dpars,eps,S);
