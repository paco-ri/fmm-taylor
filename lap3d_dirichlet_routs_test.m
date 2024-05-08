S = surfer.sphere(6,1,2,11);
rhs = ones([S.npts 1]);
zpars = [1.0+0j; 0];

p = lap3d.dirichlet.eval(S,zpars,rhs,eps);
