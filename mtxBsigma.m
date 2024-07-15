function Bsigma = mtxBsigma(S,dom,sigma,zk,eps,varargin)
%MTXBSIGMA compute sigma-dep. terms of surface magnetic field
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * sigma: [surfacefun] density for which 
%                  sigma/2 + n . grad S0[sigma]
%              is computed
%     * zk: [dcomplex] wavenumber
%     * eps: [double] precision requested
% 
%   Optional arguments:
%     * targinfo: target info (optional)
%         targinfo.r = (3,nt) target locations
%         targinfo.du = u tangential derivative info
%         targinfo.dv = v tangential derivative info
%         targinfo.n = normal info
%         targinfo.patch_id (nt,) patch id of target, = -1, if target
%           is off-surface (optional)
%         targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%           patch if on-surface (optional)
%    * Q: precomputed quadrature corrections struct (optional)
%           currently only supports quadrature corrections
%           computed in rsc format 
%    * opts: options struct
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)

if nargin < 8
    opts = [];
    opts.format = 'rsc';
else
    opts = varargin{3};
end

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals.';

% evaulate layer potential
if zk == 0
    gradSsigma = taylor.static.eval_gradS0(S,sigmavals,eps,varargin{:});
else
    gradSsigma = taylor.dynamic.eval_gradSk(S,zk,sigmavals,eps,varargin{:});
end
gradSsigma = array_to_surfacefun(gradSsigma.',dom,S); % note transpose 
n = normal(dom);
ngradSsigma = dot(n,gradSsigma);

% construct Bsigma
if zk ~= 0
    % compute m0
    m0 = debyem0(sigma,zk);
    m0vals = surfacefun_to_array(m0,dom,S);

    % compute n . Sk[m0]
    zpars = [zk 1.0 0];
    Qhelm = helm3d.dirichlet.get_quadrature_correction(S,zpars,eps, ...
        varargin{1},opts);
    Sm0 = zeros(size(m0vals));
    for j=1:3
        Sm0(:,j) = helm3d.dirichlet.eval(S,zpars,m0vals(:,j),eps, ...
            varargin{1},Qhelm,opts);
    end
    Sm0 = array_to_surfacefun(Sm0,dom,S);

    % compute n . curl Sk[m0]
    Qhelm = taylor.dynamic.get_quadrature_correction(S,zk,eps, ...
        varargin{1},opts);
    curlSm0 = taylor.dynamic.eval_curlSk(S,zk,m0vals.',eps,varargin{1}, ...
        Qhelm,opts);
    curlSm0 = array_to_surfacefun(curlSm0.',dom,S);

    % combine
    m0terms = 1i.*dot(n,zk.*Sm0+curlSm0);
    Bsigma = sigma./2 + ngradSsigma - m0terms;
else
    Bsigma = sigma./2 + ngradSsigma;
end

end
