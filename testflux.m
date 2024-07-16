% 1 Feb 24 ref Taylor state flux test
% should S0j0 technique for the flux even work when B0.n is not zero?
% test "manual" computation on different cross-sections

% define surface 
n = 3; % polynomial order 
nu = 8; 
nv = nu*3;
dom = circulartorus(n,nu,nv);

docirctorus = true;

if docirctorus
    lambda = 1.0;
    % compute reference Taylor state
    rmaj = 5.0;
    rmin = 3.0; % radius of current ring
    ntheta = 1e2;
    jmag = 1.0;
    B0 = reftaylorsurffun(dom,n,nu,nv,ntheta,rmin,rmaj,jmag,lambda);
    
    % integrate on a disc
    nr = 4*6; % number of disc pts in radial dir
    nt = 2*6; % '' in angular dir
    surfint = 0;
    domrmaj = 5.0;
    domrmin = 2.0;
    plot(dom)
    hold on
    flux = computefluxB0(domrmaj,domrmin,nr,nt,ntheta,rmin,rmaj,jmag,lambda);
    
    fprintf('nu = %d\n', nu)
    fprintf('lambda = %f\n', lambda)
    fprintf('flux, circ. torus = %f\n', flux)
end

% dom = surfacemesh.torus(n, nu, nv);
% % get A-cycle points
% xonapatches = cell(nu,3);
% whichcycle = 1;
% for i = 1:nu
%     xonapatches{i,1} = dom.x{(i-1)*nv+whichcycle};
%     xonapatches{i,2} = dom.y{(i-1)*nv+whichcycle};
%     xonapatches{i,3} = dom.z{(i-1)*nv+whichcycle};
% end
% 
% x = zeros(n*nu,3); % x on A-cycle
% for i = 1:nu
%     x((i-1)*n+1:i*n,:) = [xonapatches{i,1}(1,:);
%         xonapatches{i,2}(1,:); 
%         xonapatches{i,3}(1,:)].';
% end

% compute reference Taylor state
% rmaj = 4.5;
% rmin = 2.0; % radius of current ring
% ntheta = 1e2;
% jmag = 1.0;
% B0 = reftaylorsurffun(dom,n,nu,nv,ntheta,rmin,rmaj,jmag,lambda);
% 
% % integrate on a disc
% nr = 6; % number of disc pts in radial dir
% nt = 30; % '' in angular dir
% surfint = 0;
% domrmaj = 4.5; % d_10 in surfacemesh.torus
% B0evals = zeros([nr nt]);
% tt = 2*pi.*(0:nt-1)./nt;
% wt = 2*pi/nt;
% % plot(dom)
% % hold on
% [rr,wr] = chebpts(nr,[0,1]);
% for j = 1:nt
%     % disp(evalTorus(tt(j),0))
%     [edgeptx, edgepty, edgeptz] = evalTorus(tt(j),0);
%     edgept = [edgeptx edgepty edgeptz];
%     for i = 1:nr
%         % targpt = [domrmaj+rr(i)*sin(tt(j)) 0 rr(i)*cos(tt(j))];
%         targpt = rr(i).*(edgept - [domrmaj 0 0]) + [domrmaj 0 0];
%         plot3(targpt(1),targpt(2),targpt(3),'ob')
%         hold on
%         rmax = norm(edgept - [domrmaj 0 0]);
%         wr = rmax.*wr;
%         B0eval = reftaylor(ntheta,rmin,rmaj,jmag,lambda,targpt);
%         B0evals(i,j) = -B0eval(2);
%         surfint = surfint - B0eval(2)*rr(i)*wt*wr(i);
%         %surfint = surfint - rr(i)*wt*wr(i);
%     end
% end
% 
% plot3(x(:,1),x(:,2),x(:,3),'ro-')
% 
% fprintf('nu = %d\n', nu)
% fprintf('lambda = %f\n', lambda)
% fprintf('flux, circ. torus = %f\n', surfint)

function [x, y, z] = evalTorus(u, v)

d = [0.17 0.11  0    0;
     0    1     0.01 0;
     0    4.5   0    0;
     0   -0.25 -0.45 0];

x = zeros(size(u));
y = zeros(size(u));
z = zeros(size(u));

for i = -1:2
    for j = -1:2
        ii = i+2;
        jj = j+2;
        x = x + d(ii,jj)*cos(v).*cos((1-i)*u+j*v);
        y = y + d(ii,jj)*sin(v).*cos((1-i)*u+j*v);
        z = z + d(ii,jj)*sin((1-i)*u+j*v);
    end
end
 
end
