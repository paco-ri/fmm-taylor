function Bsigma = mtxBsigma(S,dom,sigma,eps,varargin)
%MTXBSIGMA compute sigma-dep. terms of surface magnetic field
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * sigma: [surfacefun] density for which 
%                  sigma/2 + n . grad S0[sigma]
%              is computed
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

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals.';

% evaulate layer potential
% tic
gradS0sigma = taylor.static.eval_gradS0(S,sigmavals,eps,varargin{:});
% toc
gradS0sigma = array_to_surfacefun(gradS0sigma.',dom,S); % note transpose 
n = normal(dom);
ngradS0sigma = dot(n,gradS0sigma);

% construct Bsigma
Bsigma = sigma./2 + ngradS0sigma;
% if eval_gradS0 is average of interior and exterior limit
% Bsigma = -sigma./2 + ngradS0sigma; 
% if eval_gradS0 is exterior limit
% Bsigma = sigma + ngradS0sigma;
% no identity term
% Bsigma = ngradS0sigma;

end
