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
sigmavals = sigmavals';

if (nargin > 6)
    Q = varargin{1};
elseif zk == 0
    Q = taylor.static.get_quadrature_correction(S,eps,targinfo,opts_quad);
else
    Q = taylor.dynamic.get_quadrature_correction(S,zk,eps,targinfo,opts_quad);
end

% Evaluate layer potential
if zk == 0
    gradSsigma = taylor.static.eval_gradS0(S,sigmavals,eps,targinfo,Q, ...
        opts_quad);
else
    gradSsigma = taylor.dynamic.eval_gradSk(S,zk,sigmavals,eps, ...
        targinfo,Q,opts_quad);
end

% Compute flux
% ASSUMES y^ IS THE NORMAL TO THE CROSS-SECTION
fluxsigma = dot(gradSsigma(2,:),weights);

end