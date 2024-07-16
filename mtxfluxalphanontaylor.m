function fluxalpha = mtxfluxalphanontaylor(S,dom,nodes,weights,mH, ...
    zk,eps,varargin)
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

mHvals = surfacefun_to_array(mH,dom,S);
mHvals = mHvals';

if (nargin > 6)
    Q = varargin{1};
elseif zk == 0
    Q = taylor.static.get_quadrature_correction(S,eps,targinfo,opts_quad);
else
    Q = taylor.dynamic.get_quadrature_correction(S,zk,eps,targinfo,opts_quad);
end

% Evaluate layer potential
if zk == 0
    mHterms = taylor.static.eval_curlS0(S,mHvals,eps,targinfo,Q,opts_quad);
else
    curlSmH = taylor.dynamic.eval_curlSk(S,zk,mHvals,eps,targinfo,Q,opts_quad);
    % compute Sk[mH]
    Qhelm = helm3d.dirichlet.get_quadrature_correction(S,eps,zk, ...
        [1.0 0],targinfo,opts_quad);
    SmH = zeros(size(nodes));
    opts_eval = [];
    opts_eval.precomp_quadrature = Qhelm;
    for j=1:3
        SmH(j,:) = helm3d.dirichlet.eval(S,mHvals(j,:),targinfo,eps, ...
            zk,[1.0 0],opts_eval);
    end
    mHterms = curlSmH + SmH;
end

% Compute flux
% ASSUMES y^ IS THE NORMAL TO THE CROSS-SECTION
fluxalpha = dot(mHterms(2,:),weights);

end