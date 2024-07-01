function fluxalpha = mtxfluxalphanontaylor(S,dom,nodes,weights,mH, ...
    eps,varargin)
%MTXFLUXALPHA compute alpha coefficient in flux condition
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * nodes: [double(3,*)] quadrature nodes for computing flux 
%              Gauss-Legendre in r, trapezoidal in theta
%     * weights: [double(*)] corresponding quadrature weights
%     * mH: [surfacefunv] density for which 
%                  \oint curl S0[mH] . da
%              is computed
%     * eps: [double] precision requested

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

mHvals = surfacefun_to_array(mH,dom,S);
mHvals = mHvals';

% Evaluate layer potential
curlS0mH = taylor.static.eval_curlS0(S,mHvals,eps,targinfo,Q,opts_quad);

% Compute flux
% ASSUMES y^ IS THE NORMAL TO THE CROSS-SECTION
fluxalpha = dot(curlS0mH(2,:),weights);

end