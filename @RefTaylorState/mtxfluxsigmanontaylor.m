function fluxsigma = mtxfluxsigmanontaylor(S,dom,params,nodes, ...
    weights,sigma,zk,epstaylor,epslh,varargin)
%MTXFLUXSIGMANONTAYLOR compute sigma-dep. terms of flux condition for
%                      a non-Taylor-state magnetic field
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * params: parameters describing cross-sectional surface [pol]
%         pol: [int] set to 1 if computing poloidal flux
%     * nodes: [double(3,*)] quadrature nodes for computing flux 
%              Gauss-Legendre in r, trapezoidal in theta
%     * weights: [double(*)] corresponding quadrature weights
%     * sigma: [surfacefun] density for which 
%                  \oint (grad S_0[sigma]) . da
%              is computed
%     * zk: [double complex] wavenumber
%     * epstaylor: [double] precision requested for taylor routines
%     * epslh: [double] precision requested for lap3d/helm3d routines
% 
%   Optional arguments:
%    * opts: options struct for taylor routines
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%        opts.precomp_quadrature - precomputed quadrature corrections 
%           struct for taylor routine evaluation
%             currently only supports quadrature corrections 
%             computed in rsc format 
%    * optslh: options struct for S_k evaluation
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%        opts.precomp_quadrature - precomputed quadrature corrections 
%           struct for S_k evaluation
%             currently only supports quadrature corrections 
%             computed in rsc format 

pol = params;
dpars = [1.0, 0];
targinfo = [];
targinfo.r = nodes;
if nargin < 9
    opts = [];
    opts.format = 'rsc';
    if abs(zk) < eps
        Q = taylor.static.get_quadrature_correction(S,epstaylor, ...
            targinfo,opts);
    else
        Q = taylor.dynamic.get_quadrature_correction(S,zk,epstaylor, ...
            targinfo,opts);
    end
    opts.precomp_quadrature = Q;
else
    opts = varargin{1};
end

if abs(zk) >= eps
    if nargin < 10
        if nargin == 9
            targinfo = [];
            targinfo.r = nodes;
        end
        optslh = [];
        optslh.format = 'rsc';
        Qlh = helm3d.dirichlet.get_quadrature_correction(S,epslh,zk, ...
            dpars,targinfo,optslh);
        optslh.precomp_quadrature = Qlh;
    else
        optslh = varargin{2};
    end
end

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals.';

% Evaluate layer potential
if abs(zk) < eps
    sigmaterms = taylor.static.eval_gradS0(S,sigmavals,epstaylor, ...
        targinfo,opts);
else
    gradSsigma = taylor.dynamic.eval_gradSk(S,zk,sigmavals,epstaylor, ...
        targinfo,opts);
    m0 = TaylorState.debyem0(sigma,zk);
    m0vals = surfacefun_to_array(m0,dom,S);
    m0vals = m0vals.';
    curlSm0 = taylor.dynamic.eval_curlSk(S,zk,m0vals,epstaylor, ...
        targinfo,opts);
    Sm0 = complex(zeros(size(nodes)));
    for j = 1:3
        Sm0(j,:) = helm3d.dirichlet.eval(S,m0vals(j,:),targinfo,epslh, ...
            zk,dpars,optslh);
    end
    sigmaterms = gradSsigma - 1i.*(zk.*Sm0 + curlSm0);
end

% Compute flux
% ASSUMES y^ IS THE NORMAL TO THE CROSS-SECTION
% fluxsigma = dot(sigmaterms(2,:),weights);
fluxsigma = sum(sigmaterms(2+pol,:).*weights);

end