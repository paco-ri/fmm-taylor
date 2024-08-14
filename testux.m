% test curl of B resulting of ux 

rmaj = 5.0;
phi = 4*pi/7;
h = 1e-5;
center = [rmaj*cos(phi)+.2 rmaj*sin(phi)-.1 .5];
targs = [center;
    center + [h 0 0];
    center - [h 0 0];
    center + [0 h 0];
    center - [0 h 0];
    center + [0 0 h];
    center - [0 0 h]];
targinfo = [];
targinfo.r = targs.';

if abs(zk) < 1e-16
    Q = taylor.static.get_quadrature_correction(S,eps,targinfo);
else
    Q = taylor.dynamic.get_quadrature_correction(S,zk,eps,targinfo);
end

Bint = interiorB(sigma,m,S,dom,zk,targinfo,eps,Q,optsqc,optsqclh);
Bint = Bint.';
curlB = [Bint(4,3)-Bint(5,3)-Bint(6,2)+Bint(7,2); % DyBz-DzBy
    Bint(6,1)-Bint(7,1)-Bint(2,3)+Bint(3,3); % DzBx-DxBz
    Bint(2,2)-Bint(3,2)-Bint(4,1)+Bint(5,1)]; % DxBy-DyBx
curlB = curlB./(2*h);
lambdaB = zk.*Bint(1,:).';
disp(norm(curlB-lambdaB))

% Compute B at some interior points
function B = interiorB(sigma,m,S,dom,zk,targinfo,eps,Q,optsqc,optsqclh)
sigmavals = surfacefun_to_array(sigma,dom,S);
mvals = surfacefun_to_array(m,dom,S);
if abs(zk) < 1e-16
    gradSksigma = taylor.static.eval_gradS0(S,sigmavals.',eps, ...
        targinfo,Q,optsqc);
    curlSkm = taylor.static.eval_curlS0(S,mvals.',eps, ...
        targinfo,Q,optsqc);

    B = -gradSksigma + 1i.*curlSkm;
else
    gradSksigma = taylor.dynamic.eval_gradSk(S,zk,sigmavals.',eps, ...
        targinfo,Q,optsqc);
    curlSkm = taylor.dynamic.eval_curlSk(S,zk,mvals.',eps, ...
        targinfo,Q,optsqc);

    dpars = [1.0, 0.0];
    Qlh = helm3d.dirichlet.get_quadrature_correction(S,eps,zk,dpars, ...
        targinfo,optsqclh);

    optslh = [];
    optslh.format = 'rsc';
    optslh.precomp_quadrature = Qlh;
    Skm1 = helm3d.dirichlet.eval(S,mvals(:,1),targinfo,eps, ...
        zk,dpars,optslh);
    Skm2 = helm3d.dirichlet.eval(S,mvals(:,2),targinfo,eps, ...
        zk,dpars,optslh);
    Skm3 = helm3d.dirichlet.eval(S,mvals(:,3),targinfo,eps, ...
        zk,dpars,optslh);
    Skm = [Skm1 Skm2 Skm3];
    Skm = Skm.';

    B = 1i.*zk.*Skm - gradSksigma + 1i.*curlSkm;
end
end