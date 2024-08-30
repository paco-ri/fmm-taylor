ro = 2.0;
ao = 1.0;%0.7;
bo = 1.0;
ri = 2.0;
ai = 0.55;%0.3;
bi = 0.55;

n = 6;
nv = 8;
nu = 3*nv;
npat = 2*nu*nv; % number of patches

dom = toroidal_shell(ro,ao,bo,ri,ai,bi,n,nu,nv);
domo = twisted_ellipse_torus(ao,ro,bo,n,nu,nv);
domi = twisted_ellipse_torus(ai,ri,bi,n,nu,nv);
doms = {domo,domi};

mesh(dom)
alpha 0.5

ntheta = 1e3;
rmin = 2.0;
rmaj = 2.0;
jmag = 1.0;
zk = 0;
B0 = reftaylorshell(dom,n,npat,ntheta,rmin,rmaj,jmag,zk);
% plot(norm(B0))
% alpha 0.5
% colorbar

nr = 6;
nt = 8;
np = 40;
[tornodes, torweights] = toroidalfluxquad(nr,nt,ro,ao,bo,ri,ai,bi);
[polnodes, polweights] = poloidalfluxquad(nr,np,ro,ao,bo,ri,ai,bi);
% plot3(tornodes(1,:),tornodes(2,:),tornodes(3,:),'.')
% hold on
% plot3(polnodes(1,:),polnodes(2,:),polnodes(3,:),'.')
% wireframe(dom)
% alpha 0.5

% compute fluxes
flux = zeros(1,2);
for j = 1:length(torweights)
    B0int = reftaylor(nt,rmin,rmaj,jmag,zk,tornodes(:,j));
    flux(1) = flux(1) + B0int(2)*torweights(j);
end
for j = 1:length(polweights)
    B0int = reftaylor(nt,rmin,rmaj,jmag,zk,polnodes(:,j));
    flux(2) = flux(2) + B0int(2)*polweights(j);
end

rts = TaylorState({domo,domi},[n nu nv],zk,flux,1e-6);
rts = rts.solve();

% interior point
intpt = 1.15*[1.51171 1.51171 .137886];

h = 1e-4;
[errB, curlB, kB] = rts.fd_test(intpt,h);
disp(errB)

function B0 = reftaylorshell(dom,n,npat,ntheta,rmin,rmaj,jmag,lambda)

B0 = surfacefunv(dom);
B0x = cell(npat,1);
B0y = cell(npat,1);
B0z = cell(npat,1);
for i = 1:npat
    B0x{i} = zeros(n);
    B0y{i} = zeros(n);
    B0z{i} = zeros(n);
    for j = 1:n
        for k = 1:n
            B0eval = reftaylor(ntheta,rmin,rmaj,jmag,lambda,...
                [dom.x{i}(k,j); dom.y{i}(k,j); dom.z{i}(k,j)]);
            B0x{i}(k,j) = B0eval(1);
            B0y{i}(k,j) = B0eval(2);
            B0z{i}(k,j) = B0eval(3);
        end
    end
end
B0.components{1} = surfacefun(B0x,dom);
B0.components{2} = surfacefun(B0y,dom);
B0.components{3} = surfacefun(B0z,dom);

end

function [qnodes, qweights] = toroidalfluxquad(nr,nt,ro,ao,bo,ri,ai,bi)
%TOROIDALFLUXQUAD Computes quadrature for toroidal cross-section of
%   toroidal_shell
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
        [go1, ~, go2] = evalTorus(0,tt,ro,ao,bo);
        [gi1, ~, gi2] = evalTorus(0,tt,ri,ai,bi);
        [dgo1, ~, dgo2] = dvEvalTorus(0,tt,ao,bo);
        [dgi1, ~, dgi2] = dvEvalTorus(0,tt,ai,bi);
        qnodes(:,ij) = (1-rr)*[gi1; 0; gi2] + rr*[go1; 0; go2];
        qweights(1,ij) = (2*pi/nt) ...
            * wr*((-gi1+go1)*((1-rr)*dgi2+rr*dgo2) ...
            - (-gi2+go2)*((1-rr)*dgi1+rr*dgo1));
    end
end

end

function [qnodes, qweights] = poloidalfluxquad(nr,np,ro,ao,bo,ri,ai,bi)
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
options = optimset('Display','off');
for i = 1:nr
    rr = rnodes(i);
    wr = rweights(i);
    for j = 1:np
        ij = (i-1)*np+j;
        pp = 2*pi*(j-1)/np;
        % For this phi pp, find angle v corresp. with z = 0 
        v0 = fsolve(@(v) evalTorusZ(pp,v,ao,bo),-pp,options);
        [go1, go2] = evalTorus(pp,v0,ro,ao,bo);
        [gi1, gi2] = evalTorus(pp,v0,ri,ai,bi);
        [dgo1, dgo2] = duEvalTorus(pp,v0,ro,ao,bo);
        [dgi1, dgi2] = duEvalTorus(pp,v0,ri,ai,bi);
        qnodes(:,ij) = (1-rr)*[gi1; gi2; 0] + rr*[go1; go2; 0];
        qweights(1,ij) = (2*pi/np) ...
            * wr*((-gi1+go1)*((1-rr)*dgi2+rr*dgo2) ...
            - (-gi2+go2)*((1-rr)*dgi1+rr*dgo1));
    end
end
end

function [x, y, z] = evalTorus(u, v, r, a, b)

R = a*cos(v).*cos(u) - b*sin(v).*sin(u) + r;
x = R.*cos(u);
y = R.*sin(u);
z = a*cos(v).*sin(u) + b*cos(u).*sin(v);

end

function z = evalTorusZ(u, v, a, b)

z = a*cos(v).*sin(u) + b*cos(u).*sin(v);

end

function [x, y, z] = duEvalTorus(u, v, r, a, b)

R = a*cos(v).*cos(u) - b*sin(v).*sin(u) + r;
drdu = -a*cos(v).*sin(u) - b*sin(v).*cos(u);
x = drdu.*cos(u) - R.*sin(u);
y = drdu.*sin(u) + R.*cos(u);
z = a*cos(v).*cos(u) - b*sin(u).*sin(v);

end

function [x, y, z] = dvEvalTorus(u, v, a, b)

drdv = -a*sin(v).*cos(u) - b*cos(v).*sin(u);

x = drdv.*cos(u);
y = drdv.*sin(u);
z = -a*sin(v).*sin(u) + b*cos(u).*cos(v);
 
end