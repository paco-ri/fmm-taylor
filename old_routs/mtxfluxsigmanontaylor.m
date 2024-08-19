function fluxsigma = mtxfluxsigmanontaylor(S,dom,nodes,weights,sigma, ...
    zk,eps,varargin)
%MTXFLUXSIGMANONTAYLOR compute sigma-dep. terms of flux condition for
%                      a non-Taylor-state magnetic field
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * nodes: [double(3,*)] quadrature nodes for computing flux 
%              Gauss-Legendre in r, trapezoidal in theta
%     * weights: [double(*)] corresponding quadrature weights
%     * sigma: [surfacefun] density for which 
%                  \oint (grad S_0[sigma]) . da
%              is computed
%     * eps: [double] precision requested
% 
%   Optional arguments:
%    * Q: precomputed quadrature corrections struct (optional)
%           currently only supports quadrature corrections
%           computed in rsc format 

% Get near quadrature correction 
targinfo = [];
targinfo.r = nodes;
opts_quad = [];
opts_quad.format = 'rsc';

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals.';

if (nargin > 6)
    Q = varargin{1};
elseif zk == 0
    Q = taylor.static.get_quadrature_correction(S,eps,targinfo,opts_quad);
else
    Q = taylor.dynamic.get_quadrature_correction(S,zk,eps,targinfo,opts_quad);
end
opts = [];
opts.format = 'rsc';
opts.precomp_quadrature = Q;

% Evaluate layer potential
if zk == 0
    sigmaterms = taylor.static.eval_gradS0(S,sigmavals,eps,targinfo,opts);
else
    gradSsigma = taylor.dynamic.eval_gradSk(S,zk,sigmavals,eps, ...
        targinfo,opts);
    m0 = debyem0(sigma,zk);
    m0vals = surfacefun_to_array(m0,dom,S);
    m0vals = m0vals.';
    curlSm0 = taylor.dynamic.eval_curlSk(S,zk,m0vals,eps,targinfo,opts);
    Qhelm = helm3d.dirichlet.get_quadrature_correction(S,eps,zk, ...
        [1.0 0],targinfo,opts_quad); 
    opts_helm = [];
    opts_helm.precomp_quadrature = Qhelm;
    Sm0 = complex(zeros(size(nodes)));
    for j = 1:3
        Sm0(j,:) = helm3d.dirichlet.eval(S,m0vals(j,:),targinfo,eps,zk, ...
            [1.0 0],opts_helm);
    end
    sigmaterms = gradSsigma - 1i.*(zk.*Sm0 + curlSm0);
end

% Compute flux
% ASSUMES y^ IS THE NORMAL TO THE CROSS-SECTION
% fluxsigma = dot(sigmaterms(2,:),weights);
fluxsigma = sum(sigmaterms(2,:).*weights);

end