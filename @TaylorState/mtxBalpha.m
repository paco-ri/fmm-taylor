function Balpha = mtxBalpha(S,dom,domparams,mH,zk,epstaylor,epslh,varargin)
%MTXBALPHA compute alpha coefficient for surface magnetic field
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * domparams: parameter describing dom and circulation [io]
%         io: [int] if 1, negate vn because inner torus
%     * mH: [surfacefunv] density for which 
%                  n . curl S0[mH]
%           is computed. Note that there is no factor of i. 
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
        Q = taylor.static.get_quadrature_correction(S,epstaylor, ...
            targinfo);
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

mHvals = surfacefun_to_array(mH,dom,S);
mHvals = mHvals.';

% evaluate layer potential
vn = (-1)^io.*normal(dom);

if abs(zk) > eps
    % compute curl Sk[mH]
    curlSmH = taylor.dynamic.eval_curlSk(S,zk,mHvals,epstaylor, ...
        targinfo,opts);
    curlSmH = array_to_surfacefun(curlSmH.',dom,S);

    % compute Sk[mH]
    SmH = complex(zeros(size(mHvals)));
    mHvals = surfacefun_to_array(mH,dom,S); 
    mHvals = mHvals.';
    for j=1:3
        SmH(j,:) = helm3d.dirichlet.eval(S,mHvals(j,:),targinfo, ...
            epslh,zk,dpars,optslh);
    end
    SmH = array_to_surfacefun(SmH.',dom,S);

    % Balpha = dot(n,zk.*SmH + curlSmH);
    Balpha = zk.*dot(vn,SmH) + dot(vn,curlSmH);
else
    curlSmH = taylor.static.eval_curlS0(S,mHvals,epstaylor,varargin{:});
    curlSmH = array_to_surfacefun(curlSmH.',dom,S);
    Balpha = dot(vn,curlSmH);
end

end