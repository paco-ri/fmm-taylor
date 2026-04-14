clear
addpath('../examples')

n = 8; 
nv = 10;
nu = nv*3;
rmaj = 4.0;
rmin = 2.0;
dom = circulartorus(n,nu,nv,rmin,rmaj);

% mark some patches for refinement
marked = 1:10;

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

% function to integrate: azimuthal 1/r field around z-axis
% integral along one toroidal loop is +/-2*pi depending on orientation
f = surfacefunv(@(x,y,z) -y./(x.^2 + y.^2), ...
    @(x,y,z)  x./(x.^2 + y.^2), ...
    @(x,y,z)  0*x, dom);

if ~do_adap_ref || all(p2q(:,3) == 0)
    [integralb, quadpts, quadwts] = TaylorState.intbcyc(f, n, nu, nv);
else
    [integralb, quadpts, quadwts] = TaylorState.intbcyc(f, n, nu, nv, qf, p2q);
end
% hold on
% plot3(quadpts(:,1),quadpts(:,2),quadpts(:,3),'ro')
surfacemesh_to_vtk(dom, 'testintbcycadap_quadpts.vtk', [], Points=quadpts, ...
    Title="adaptive intbcyc mesh and quadrature points");

error = abs(integralb - 2*pi);
if (error < 10*amr_tol)
    fprintf("intbcyc test passed! \terror = %.6e\n", error)
else
    fprintf("intbcyc test failed! \terror = %.6e, integralb = %.6e, expected = %.6e", error, integralb, 2*pi*(rmaj + rmin)*(rmaj + rmin))
end