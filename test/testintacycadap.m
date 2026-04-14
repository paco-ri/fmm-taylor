clear
addpath('../examples')

n = 7; 
nv = 10;
nu = nv*3;
rmaj = 4.0;
rmin = 2.0;
dom = circulartorus(n,nu,nv,rmin,rmaj);

% mark some patches for refinement
marked = 1:8;

% initialize amr parameters
amr_tol = 1e-10;
rmax = 4;
mode = 1;

% refine the mesh
do_adap_ref = true;
if do_adap_ref
    [dom, qf, p2q] = surfacemesh.adap_ref(dom, amr_tol, rmax, mode, marked);
end
% plot(dom)

% function to integrate
% cos phi = x/(rmaj + rmin cos theta), sin phi = y/(rmaj + rmin cos theta)
costheta = @(x,y,z) sign(x.^2+y.^2-rmaj^2).*sqrt(1 - z.^2./(rmin^2));
sintheta = @(x,y,z) z/rmin;
cosphi = @(x,y,z) x./(rmaj + rmin*costheta(x,y,z));
sinphi = @(x,y,z) y./(rmaj + rmin*costheta(x,y,z));
f = surfacefunv(@(x,y,z) -sintheta(x,y,z).*cosphi(x,y,z), ...
    @(x,y,z) -sintheta(x,y,z).*sinphi(x,y,z), ...
    @(x,y,z) costheta(x,y,z), dom);

% if all the entries in p2q(:,3) are zero
if ~do_adap_ref || all(p2q(:,3) == 0)
    [integrala, quadpts, quadwts] = TaylorState.intacyc(f, n, nu, nv);
else
    [integrala, quadpts, quadwts] = TaylorState.intacyc(f, n, nu, nv, qf, p2q);
end
% hold on
% plot3(quadpts(:,1),quadpts(:,2),quadpts(:,3),'ro')
surfacemesh_to_vtk(dom, 'testintacycadap_quadpts.vtk', [], Points=quadpts, ...
    Title="adaptive intacyc mesh and quadrature points");

if (abs(integrala - 4 * pi) < 10*amr_tol)
    fprintf("intacyc test passed! \terror = %.6e\n", abs(integrala - 4 * pi))
else
    fprintf("intacyc test failed! \terror = %.6e\n", abs(integrala - 4 * pi))
end