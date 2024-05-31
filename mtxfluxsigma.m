function fluxsigma = mtxfluxsigma(S,dom,n,nu,nv,sigma,eps,varargin)
%MTXFLUXSIGMA compute sigma-dep. terms of flux condition
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * n: [int] polynomial order on each surface patch
%     * nu: [int] number of patches in toroidal direction
%     * nv: [int] number of patches in poloidal direction
%     * sigma: [surfacefun] density for which 
%                  \oint S0[n times grad S0[sigma]]
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
sigmavals = sigmavals';

% Evalaute layer potentials
gradS0sigma = taylor.static.eval_gradS0(S,sigmavals,eps,varargin{:});
gradS0sigma = array_to_surfacefun(gradS0sigma',dom,S); % note transpose
vn = normal(dom);
nxgradS0sigma = cross(vn,gradS0sigma);
nxvals = surfacefun_to_array(nxgradS0sigma,dom,S);

% single layer potential evals
dpars = [1.0, 0.0];
% TODO: allow user to input quadrature 
S0nx1 = lap3d.dirichlet.eval(S,dpars,nxvals(:,1),eps);
S0nx2 = lap3d.dirichlet.eval(S,dpars,nxvals(:,2),eps);
S0nx3 = lap3d.dirichlet.eval(S,dpars,nxvals(:,3),eps);
S0nx = [S0nx1 S0nx2 S0nx3];
S0nx = array_to_surfacefun(S0nx,dom,S);

fluxsigma = intacyc(S0nx,n,nu,nv);

end