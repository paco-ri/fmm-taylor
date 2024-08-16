% --- Geometry parameters ---
n = 8; % polynomial order + 1
nv = 10; % number of patches in poloidal direction
nu = nv*3; % number of patches in toroidal direction
a = 2.0; % minor radius, horiz. axis
a0 = 5.0; % major radius
b = 3.0; % minor radius vert. axis
dom = twisted_ellipse_torus(a,a0,b,n,nu,nv); % surfacemesh
domparams = [n, nu, nv];

% --- Specify flux ---
flux = 1.0;

% --- Beltrami parameter ---
zk = 0.5;

% --- Tolerances ---
tol = 1e-12;

ts = TaylorState(dom,domparams,zk,flux,tol);
ts = ts.solve(true);
B = ts.surface_B();

rmaj = 5.0;
phi = 4*pi/7;
h = 1e-6;
center = [rmaj*cos(phi)+.2 rmaj*sin(phi)-.1 .5];

[errB, curlB, kB] = ts.fd_test(center,h);
disp(errB)
