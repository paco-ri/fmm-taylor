function fluxsigma = mtxfluxsigmanontaylor(S,dom,nodes,weights,sigma, ...
    eps,varargin)
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
if (nargin > 6)
    Q = varargin{1};
else
    Q = taylor.static.get_quadrature_correction(S,eps,targinfo,opts_quad);
end

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals';

% Evaluate layer potential
gradS0sigma = taylor.static.eval_gradS0(S,sigmavals,eps,targinfo,Q, ...
    opts_quad);

% Compute flux
% ASSUMES y^ IS THE NORMAL TO THE CROSS-SECTION
fluxsigma = dot(gradS0sigma(2,:),weights);

end