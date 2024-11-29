% --- Tolerances ---
tol = 1e-4;

% --- Beltrami parameter ---
zk = 0;

% --- Geometry parameters ---
ns = [5];% 7 9]; % polynomial order + 1
nvs = [4 6];% 8 10]; % number of patches in poloidal direction

lerr = zeros(4,size(ns,2)*size(nvs,2));
lind = 1;

for n = ns
    for nv = nvs
        fprintf('\t========\n\tn = %d, nv = %d \n', n, nv)
        nu = nv*3; % number of patches in toroidal direction
        r = 2.0; % major radius
        ao = 1.0; % outer minor radius
        ai = 0.6; % inner minor radius
        domo = circulartorus(n,nu,nv,ao,r);
        domi = circulartorus(n,nu,nv,ai,r);
        dom = {domo,domi};
        domparams = [n, nu, nv];

        % --- Compute B0 ---
        ntheta = 1e3;
        rmin = 2.0;
        rmaj = 2.0;
        jmag = 1.0;
        B0o = reftaylorsurffun(domo,n,nu,nv,ntheta,rmin,rmaj,jmag,zk);
        B0i = reftaylorsurffun(domi,n,nu,nv,ntheta,rmin,rmaj,jmag,zk);

        % --- Compute XS quad. and flux ---
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
            polflux = polflux + B0eval(3)*polweights(i); % note sign
        end
        
        B0 = {B0o,B0i};
        qnodes = {tornodes,polnodes};
        qweights = {torweights,polweights};
        flux = [torflux,polflux];
        ts = RefTaylorState(dom,domparams,zk,flux,B0,qnodes,qweights,tol);
        ts = ts.solve(true);
        B = ts.surface_B();

        lerr(1,lind) = n; % number of points on each patch
        lerr(2,lind) = nv; % geom. param. 
        lerr(3,lind) = nu; % geom. param. 
        lerr(4,lind) = max(vecinfnorm(B0{1}-B{1}),vecinfnorm(B0{2}-B{2}))...
            /max(vecinfnorm(B0{1}),vecinfnorm(B0{2})); 

        lind = lind + 1;
    end
    tol = tol*1e-2;
end

loglog(sqrt(lerr(1,1:size(nvs,2)).^2.*lerr(2,1:size(nvs,2)).* ...
    lerr(3,1:size(nvs,2))), lerr(4,1:size(nvs,2)), 'o-')
hold on
for i = 1:size(ns,2)-1
    ind = i*size(nvs,2)+1:(i+1)*size(nvs,2);
    loglog(sqrt(lerr(1,ind).^2.*lerr(2,ind).* ...
        lerr(3,ind)), lerr(4,ind), 'o-')
end

h1 = 50:90;
h2 = 80:120;
h3 = 100:150;
loglog(h1,5e1.*h1.^(-2),'--')
loglog(h2,1e4.*h2.^(-4),'--')
loglog(h3,1e7.*h3.^(-6),'--')

legend('p=poly. order=4','p=6','p=8','O(h^{-2})','O(h^{-4})','O(h^{-6})','location','southwest')
xlabel('h, panel size')
ylabel('||B-B_0||_\infty/||B_0||_\infty')

function N = vecinfnorm(f)
N = max([norm(f.components{1}, inf) ...
        norm(f.components{2}, inf), ...
        norm(f.components{3}, inf)]);
end

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
%     qnodes: [double(3,nr*np)] quadrature nodes
%     qweights: [double(1,nr*np)] quadrature weights

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

