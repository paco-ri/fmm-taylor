function fluxalpha = mtxfluxalpha(S,dom,n,nu,nv,mH,eps,varargin)
%MTXFLUXALPHA compute alpha coefficient in flux condition
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * n: [int] polynomial order on each surface patch
%     * nu: [int] number of patches in toroidal direction
%     * nv: [int] number of patches in poloidal direction
%     * mH: [surfacefunv] density for which 
%                  \oint -mH/2 + n x curl S0[mH]
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

mHvals = surfacefun_to_array(mH,dom,S);
mHvals = mHvals';

% derivative of layer potential
curlS0mH = taylor.static.eval_curlS0(S,mHvals,eps,varargin{:});
curlS0mH = array_to_surfacefun(curlS0mH',dom,S);
vn = normal(dom);
nxcurlS0mH = cross(vn,curlS0mH);
nxminusmH = nxcurlS0mH - mH./2;
nxminusmHvals = surfacefun_to_array(nxminusmH,dom,S);

% single layer potential evals
dpars = [1.0, 0.0];
% TODO: allow user to input Laplace quadrature 
Qlap = lap3d.dirichlet.get_quadrature_correction(S,dpars,eps,S);
S0nx1 = lap3d.dirichlet.eval(S,dpars,nxminusmHvals(:,1),eps,S,Qlap);
S0nx2 = lap3d.dirichlet.eval(S,dpars,nxminusmHvals(:,2),eps,S,Qlap);
S0nx3 = lap3d.dirichlet.eval(S,dpars,nxminusmHvals(:,3),eps,S,Qlap);
S0nx = [S0nx1 S0nx2 S0nx3];
S0nx = array_to_surfacefun(S0nx,dom,S);

fluxalpha = intacyc(S0nx,n,nu,nv);

end