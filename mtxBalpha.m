function Balpha = mtxBalpha(S,dom,mH,eps,varargin)
%MTXBALPHA compute alpha coefficient for surface magnetic field
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * mH: [surfacefunv] density for which 
%                  n . curl S0[mH]
%           is computed. Note that there is no factor of i. 
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

mHvals = surfacefun_to_array(mH,dom,S);
mHvals = mHvals.';

% evaluate layer potential
n = normal(dom);
curlS0mH = taylor.static.eval_curlS0(S,mHvals,eps,varargin{:});
curlS0mH = array_to_surfacefun(curlS0mH.',dom,S);
Balpha = dot(n,curlS0mH);

end