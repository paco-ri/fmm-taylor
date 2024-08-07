function m0 = debyem0(sigma,lambda)
%DEBYEM0 compute part of the Debye current associated with sigma
%   Detailed explanation goes here

pdo = [];
pdo.lap = 1;

% resample
n = size(sigma.vals{1,1},1);
sigma = resample(sigma,n+4);
dom = sigma.domain;
S = surfer.surfacemesh_to_surfer(dom);

% solve lap(u) = sigma
L = surfaceop(dom, pdo, sigma);
L.rankdef = true;
u = L.solve();

vn = normal(dom);
% ==== vn computation ====
% sinphi = @(x,y,z) y./sqrt(x.^2 + y.^2);
% cosphi = @(x,y,z) x./sqrt(x.^2 + y.^2);
% domrmin = 1.0;
% domrmaj = 2.0;
% sintheta = @(x,y,z) z./domrmin; 
% costheta = @(x,y,z) (sqrt(x.^2 + y.^2) - domrmaj)./domrmin;
% vn = surfacefunv(@(x,y,z) costheta(x,y,z).*cosphi(x,y,z), ...
%                  @(x,y,z) costheta(x,y,z).*sinphi(x,y,z), ...
%                  @(x,y,z) sintheta(x,y,z), dom); 
% ========================
% m0 = 1i.*lambda.*(grad(u) + 1i.*cross(vn, grad(u)));
uvals = surfacefun_to_array(u,dom,S);
graduvals = get_surface_grad(S,uvals);
gradu = array_to_surfacefun(graduvals.',dom,S);
m0 = 1i.*lambda.*(gradu + 1i.*cross(vn, gradu));

% m0err = div(m0) - 1i*lambda.*sigma;
% fprintf('upsampled m0err = %f\n', norm(m0err))

% resample
m0 = resample(m0,n);

end