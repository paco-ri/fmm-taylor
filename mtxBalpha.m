function Balpha = mtxBalpha(S,dom,mH,zk,eps,varargin)
%MTXBALPHA compute alpha coefficient for surface magnetic field
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * mH: [surfacefunv] density for which 
%                  n . curl S0[mH]
%           is computed. Note that there is no factor of i. 
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

mHvals = surfacefun_to_array(mH,dom,S);
mHvals = mHvals.';

% evaluate layer potential
n = normal(dom);

if zk ~= 0
    % compute curl Sk[mH]
    curlSmH = taylor.dynamic.eval_curlSk(S,zk,mHvals,eps,varargin{:});
    curlSmH = array_to_surfacefun(curlSmH.',dom,S);

    % compute Sk[mH]
    Qhelm = helm3d.dirichlet.get_quadrature_correction(S,eps,zk, ...
        [1.0 0],varargin{1},opts);
    SmH = zeros(size(mHvals));
    opts_eval = [];
    opts_eval.precomp_quadrature = Qhelm;
    opts_eval.format = 'rsc';
    for j=1:3
        SmH(j,:) = helm3d.dirichlet.eval(S,mHvals(j,:),varargin{1},eps, ...
            zk,[1.0 0],opts_eval);
    end
    SmH = array_to_surfacefun(SmH.',dom,S);

    Balpha = dot(n,zk.*SmH + curlSmH);
else
    curlSmH = taylor.static.eval_curlS0(S,mHvals,eps,varargin{:});
    curlSmH = array_to_surfacefun(curlSmH.',dom,S);
    Balpha = dot(n,curlSmH);
end

end