function fluxalpha = mtxfluxalphanontaylor(S,dom,params,nodes,weights, ...
    mH,zk,epstaylor,epslh,varargin)
%MTXFLUXALPHANONTAYLOR compute alpha coefficient in flux condition for
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
%     * mH: [surfacefunv] density for which 
%                  \oint curl S0[mH] . da
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
    targinfo = [];
    targinfo.r = nodes;
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

mHvals = surfacefun_to_array(mH,dom,S);
mHvals = mHvals.';

% Evaluate layer potential
if abs(zk) < eps
    mHterms = taylor.static.eval_curlS0(S,mHvals,epstaylor,targinfo,opts);
else
    curlSmH = taylor.dynamic.eval_curlSk(S,zk,mHvals,epstaylor, ...
        targinfo,opts);
    % compute Sk[mH]
    SmH = complex(zeros(size(nodes)));
    for j=1:3
        SmH(j,:) = helm3d.dirichlet.eval(S,mHvals(j,:),targinfo,epslh, ...
            zk,[1.0 0],optslh);
    end
    mHterms = curlSmH + zk.*SmH;
end

% Compute flux
% ASSUMES y^ IS THE NORMAL TO THE CROSS-SECTION
% fluxalpha = dot(mHterms(2,:),weights);
fluxalpha = sum(mHterms(2+pol,:).*weights);

end