function m0 = debyem0(sigma,lambda,L,vn)
%DEBYEM0 compute part of the Debye current associated with sigma
%   
%    Arguments
%      sigma [surfacefun] 
%      lambda [double] 
%      L [surfaceop]
%      vn [surfacefun] normal

% moving oversampling outside of this, orders of sigma and L must match

% oversample = false;
% if nargin > 2
%     oversample = true;
%     np = varargin{1};
% end

% pdo = [];
% pdo.lap = 1;

% % resample
% n = size(sigma.vals{1,1},1);
% if oversample
%     sigma = resample(sigma,n+np);
% end
% 
% dom = sigma.domain;

% % solve lap(u) = sigma
% L = surfaceop(dom, pdo, sigma);
% L.rankdef = true;
L.rhs = sigma;
u = L.solve();
m0 = 1i.*lambda.*(grad(u) + 1i.*cross(vn, grad(u)));

% if oversample
%     % resample
%     m0 = resample(m0,n);
% end

end