function fluxalpha = mtxfluxalpha(S,dom,domparams,mH,zk,epstaylor,epslh,varargin)
%MTXFLUXALPHA compute alpha coefficient in flux condition
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * domparams: parameters describing dom [n, nu, nv]
%         n: [int] polynomial order on each surface patch
%         nu: [int] number of patches in toroidal direction
%         nv: [int] number of patches in poloidal direction
%     * mH: [surfacefunv] density for which 
%                  \oint -mH/2 + n x curl S0[mH]
%              is computed
%     * zk: [double complex] wavenumber
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

n = domparams(1);
nu = domparams(2);
nv = domparams(3);
dpars = [1.0, 0.0];
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
    if abs(zk) < 1e-16
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
mHvals = mHvals';
vn = normal(dom);

if abs(zk) < eps
    % n x curl S0[mH] - mH/2
    curlS0mH = taylor.static.eval_curlS0(S,mHvals,epstaylor,targinfo,opts);
    curlS0mH = array_to_surfacefun(curlS0mH',dom,S);
    nxcurlS0mH = cross(vn,curlS0mH);
    nxminusmH = nxcurlS0mH - mH./2;
    nxminusmHvals = surfacefun_to_array(nxminusmH,dom,S);
    
    % S0[n x curl S0[mH] - mH/2]
    S0nx1 = lap3d.dirichlet.eval(S,nxminusmHvals(:,1),targinfo,epslh, ...
        dpars,optslh);
    S0nx2 = lap3d.dirichlet.eval(S,nxminusmHvals(:,2),targinfo,epslh, ...
        dpars,optslh);
    S0nx3 = lap3d.dirichlet.eval(S,nxminusmHvals(:,3),targinfo,epslh, ...
        dpars,optslh);
    S0nx = [S0nx1 S0nx2 S0nx3];
    S0nx = array_to_surfacefun(S0nx,dom,S);

    fluxalpha = 1i*TaylorState.intacyc(S0nx,n,nu,nv);
else
    % Sk[mH]
    SkmH1 = helm3d.dirichlet.eval(S,mHvals(1,:),targinfo,epslh,zk, ...
        dpars,optslh);
    SkmH2 = helm3d.dirichlet.eval(S,mHvals(2,:),targinfo,epslh,zk, ...
        dpars,optslh);
    SkmH3 = helm3d.dirichlet.eval(S,mHvals(3,:),targinfo,epslh,zk, ...
        dpars,optslh);
    SkmH = [SkmH1 SkmH2 SkmH3];
    SkmH = array_to_surfacefun(SkmH,dom,S);

    % curl Sk[mH]
    curlSkmH = taylor.dynamic.eval_curlSk(S,zk,mHvals,epstaylor, ...
        targinfo,opts);
    curlSkmH = array_to_surfacefun(curlSkmH.',dom,S);

    % curl S0[mH]
    curlS0mH = taylor.static.eval_curlS0(S,mHvals,eps);
    curlS0mH = array_to_surfacefun(curlS0mH.',dom,S);

    % A-cycle integral
    integrand = 1i.*( SkmH + (curlSkmH - curlS0mH)./zk );
    fluxalpha = TaylorState.intacyc(integrand,n,nu,nv);
end

end