function Bsigma = mtxBsigma(S,dom,domparams,sigma,zk,epstaylor,epslh,varargin)
%MTXBSIGMA compute sigma-dep. terms of surface magnetic field
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * domparams: parameter describing dom and circulation [io]
%         io: [int] if 1, negate vn because inner torus
%     * sigma: [surfacefun] density for which 
%                  sigma/2 + n . grad S0[sigma]
%              is computed
%     * zk: [dcomplex] wavenumber
%     * epstaylor: [double] precision requested for taylor routines
%     * epslh: [double] precision requested for lap3d/helm3d routines
% 
%   Optional arguments:
%     * targinfo: target info
%         targinfo.r = (3,nt) target locations
%         targinfo.du = u tangential derivative info
%         targinfo.dv = v tangential derivative info
%         targinfo.n = normal info
%         targinfo.patch_id (nt,) patch id of target, = -1, if target
%           is off-surface (optional)
%         targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%           patch if on-surface (optional)
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

io = domparams;
dpars = [1.0,0];

if nargin < 8
    targinfo = S;
else
    targinfo = varargin{1};
end

if nargin < 9
    if abs(zk) < eps
        Q = taylor.static.get_quadrature_correction(S,epstaylor,targinfo);
    else
        Q = taylor.dynamic.get_quadrature_correction(S,zk, ...
            epstaylor,targinfo);
    end
    opts = [];
    opts.format = 'rsc';
    opts.precomp_quadrature = Q;
else
    opts = varargin{2};
end

if nargin < 10
    if abs(zk) < eps
        Qlh = lap3d.dirichlet.get_quadrature_correction(S, ...
            epslh,dpars,targinfo,opts);
    else
        Qlh = helm3d.dirichlet.get_quadrature_correction(S, ...
            epslh,zk,dpars,targinfo,opts);
    end
    optslh = [];
    optslh.format = 'rsc';
    optslh.precomp_quadrature = Qlh;
else
    optslh = varargin{3};
end

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals.';

% evaulate layer potential
if abs(zk) < eps
    gradSsigma = taylor.static.eval_gradS0(S,sigmavals,epstaylor, ...
        targinfo,opts);
else
    gradSsigma = taylor.dynamic.eval_gradSk(S,zk,sigmavals,epstaylor, ...
        targinfo,opts);
end
gradSsigma = array_to_surfacefun(gradSsigma.',dom,S); % note transpose 
vn = (-1)^io.*normal(dom);
ngradSsigma = dot(vn,gradSsigma);

% construct Bsigma
if abs(zk) > eps
    % compute m0
    m0 = TaylorState.debyem0(sigma,zk);
    m0vals = surfacefun_to_array(m0,dom,S);

    % compute n . Sk[m0]
    Sm0 = complex(zeros(size(m0vals)));
    for j=1:3
        Sm0(:,j) = helm3d.dirichlet.eval(S,m0vals(:,j),targinfo,epslh, ...
            zk,[1.0 0],optslh);
    end
    Sm0 = array_to_surfacefun(Sm0,dom,S);

    % compute n . curl Sk[m0]
    curlSm0 = taylor.dynamic.eval_curlSk(S,zk,m0vals.',epstaylor, ...
        targinfo,opts);
    curlSm0 = array_to_surfacefun(curlSm0.',dom,S);

    % combine
    % m0terms = 1i.*dot(n,zk.*Sm0+curlSm0);
    m0terms = 1i*zk.*dot(vn,Sm0) + 1i.*dot(vn,curlSm0);
    Bsigma = sigma./2 + ngradSsigma - m0terms;
else
    Bsigma = sigma./2 + ngradSsigma;
end

end
