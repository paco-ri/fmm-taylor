r = 2.0;
ao = 1.0;
ai = 0.6;

n = 5;
nv = 5;
nu = 3*nv;

domo = circulartorus(n,nu,nv,ao,r);
domi = circulartorus(n,nu,nv,ai,r);
doms = {domo, domi};

zk = 0;%.5;

% Reference Taylor State
ntheta = 1e3;
rmin = 2.0;
rmaj = r;
jmag = 1.0;
B0o = reftaylorsurffun(domo,n,nu,nv,ntheta,rmin,rmaj,jmag,zk);
B0i = reftaylorsurffun(domi,n,nu,nv,ntheta,rmin,rmaj,jmag,zk);

% plot(norm(B0o))

nr = 16;
nt = 40;
np = 40;
[tornodes, torweights] = toroidalfluxquad(nr,nt,r,ao,r,ai);
[polnodes, polweights] = poloidalfluxquad(nr,np,r,ao,r,ai);

torflux = 0;
for i = 1:nr*nt
    B0eval = reftaylor(ntheta,rmin,rmaj,jmag,zk,tornodes(:,i));
    torflux = torflux + B0eval(2)*torweights(i);
end
polflux = 0;
for i = 1:nr*np
    B0eval = reftaylor(ntheta,rmin,rmaj,jmag,zk,polnodes(:,i));
    polflux = polflux + B0eval(3)*polweights(i);
end

% plot3(tornodes(1,:),tornodes(2,:),tornodes(3,:),'.')
% hold on
% plot3(polnodes(1,:),polnodes(2,:),polnodes(3,:),'.')
% plot(domo)
% alpha .5
% plot(domi)

flux = [1.0,1.0];
tol = 1e-7;
ts = TaylorState(doms,[n,nu,nv],zk,flux,tol);
% ts = TaylorState(domo,[n,nu,nv],zk,torflux,1e-6);
% ts = TaylorState(doms,[n,nu,nv],zk,[torflux,polflux],tol);

% rts = RefTaylorState(doms,[n,nu,nv],zk,[torflux,polflux],{B0o,B0i}, ...
    % {tornodes,polnodes},{torweights,polweights},tol);
% rts = RefTaylorState(domo,[n,nu,nv],zk,torflux,{B0o,B0i}, ...
%     {tornodes},{torweights},tol);

ts = ts.solve(true);
% rts = rts.solve(true);

% B = rts.surface_B();

% plot(norm(B0o-B{1}));
% colorbar

% interior point
intpt = [2.8*cos(pi/12) 2.8*sin(pi/12) 0.01];

h = 1e-4;
[errB, curlB, kB] = ts.fd_test(intpt,h);
disp(errB)

% [errB2, curlB2, kB2] = rts.fd_test(intpt,h);
% disp(errB2)

fprintf('true flux = [%f,%f]\n', torflux, polflux)
intBtor = ts.interior_B(tornodes);
torfluxcheck = sum(intBtor(2,:).*torweights);
intBpol = ts.interior_B(polnodes);
polfluxcheck = sum(intBpol(3,:).*polweights);
fprintf('TS flux = [%f+i%f,%f+i%f]\n', real(torfluxcheck), ...
    imag(torfluxcheck), real(polfluxcheck), imag(polfluxcheck))

% intBtor = rts.interior_B(tornodes);
% torfluxcheck = sum(intBtor(2,:).*torweights);
% intBpol = rts.interior_B(polnodes);
% polfluxcheck = sum(intBpol(3,:).*polweights);
% fprintf('RTS flux = [%f+i%f,%f+i%f]\n', real(torfluxcheck), ...
%     imag(torfluxcheck), real(polfluxcheck), imag(polfluxcheck))

function [qnodes, qweights] = toroidalfluxquad(nr,nt,ro,ao,ri,ai)
%TOROIDALFLUXQUAD Computes quadrature for toroidal cross-section
%   Gauss-Legendre in r, periodic trapezoidal in theta
% 
%   Arguments:
%     nr: [int] number of r nodes
%     nt: [int] number of theta nodes
%   Returns:
%     qnodes: [double(3,nr*nt)] quadrature nodes
%     qweights: [double(1,nr*nt)] quadrature weights

[rnodes, rweights] = chebpts(nr,[0 1],1);
qnodes = zeros([3 nr*nt]);
qweights = zeros([1 nr*nt]);
for i = 1:nr
    rr = rnodes(i);
    wr = rweights(i);
    for j = 1:nt
        ij = (i-1)*nt+j;
        tt = 2*pi*(j-1)/nt;
        [go1, ~, go2] = evalTorus(0,tt,ro,ao);
        [gi1, ~, gi2] = evalTorus(0,tt,ri,ai);
        [dgo1, ~, dgo2] = dvEvalTorus(0,tt,ao);
        [dgi1, ~, dgi2] = dvEvalTorus(0,tt,ai);
        qnodes(:,ij) = (1-rr)*[gi1; 0; gi2] + rr*[go1; 0; go2];
        qweights(1,ij) = (2*pi/nt) ...
            * wr*((-gi1+go1)*((1-rr)*dgi2+rr*dgo2) ...
            - (-gi2+go2)*((1-rr)*dgi1+rr*dgo1));
    end
end

end

function [qnodes, qweights] = poloidalfluxquad(nr,np,ro,ao,ri,ai)
%POLOIDALFLUXQUAD Computes quadrature for poloidal cross-section of
%   toroidal_shell
%   Gauss-Legendre in r, periodic trapezoidal in phi
% 
%   Arguments:
%     nr: [int] number of r nodes
%     np: [int] number of phi nodes
%   Returns:
%     qnodes: [double(3,nr*nt)] quadrature nodes
%     qweights: [double(1,nr*nt)] quadrature weights

[rnodes, rweights] = chebpts(nr,[0 1],1);
qnodes = zeros([3 nr*np]);
qweights = zeros([1 nr*np]);
for i = 1:nr
    rr = rnodes(i);
    wr = rweights(i);
    for j = 1:np
        ij = (i-1)*np+j;
        pp = 2*pi*(j-1)/np;
        [go1, go2] = evalTorus(pp,0,ro,ao);
        [gi1, gi2] = evalTorus(pp,0,ri,ai);
        [dgo1, dgo2] = duEvalTorus(pp,0,ro,ao);
        [dgi1, dgi2] = duEvalTorus(pp,0,ri,ai);
        qnodes(:,ij) = (1-rr)*[gi1; gi2; 0] + rr*[go1; go2; 0];
        qweights(1,ij) = (2*pi/np) ...
            * wr*((-gi1+go1)*((1-rr)*dgi2+rr*dgo2) ...
            - (-gi2+go2)*((1-rr)*dgi1+rr*dgo1));
    end
end
end

function [x, y, z] = evalTorus(u, v, rmaj, rmin)

x = (rmaj + rmin*cos(v)).*cos(u);
y = (rmaj + rmin*cos(v)).*sin(u);
z = rmin*sin(v);
 
end

function [x, y, z] = duEvalTorus(u, v, rmaj, rmin)

x = -(rmaj + rmin*cos(v)).*sin(u);
y = (rmaj + rmin*cos(v)).*cos(u);
z = 0.*u;

end

function [x, y, z] = dvEvalTorus(u, v, rmin)

x = -rmin*sin(v).*cos(u);
y = -rmin*sin(v).*sin(u);
z = rmin*cos(v);

end
