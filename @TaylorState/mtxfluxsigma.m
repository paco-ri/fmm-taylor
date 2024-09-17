function fluxsigma = mtxfluxsigma(S,dom,L,domparams,sigma,zk,epstaylor,epslh,varargin)
%MTXFLUXSIGMA compute sigma-dep. terms of flux condition
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * domparams: parameters describing dom and circulation [n, nu, nv, io, aint]
%         n: [int] polynomial order on each surface patch
%         nu: [int] number of patches in toroidal direction
%         nv: [int] number of patches in poloidal direction
%         io: [int] if 1, negate vn because inner torus
%         aint: [int] if 1, do A-cyc. integral; otherwise, B-cyc.
%     * sigma: [surfacefun] density for which 
%                  \oint S0[n times grad S0[sigma]]
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
io = domparams(4);
aint = domparams(5);
dpars = [1.0, 0.0];
nreqarg = 8;
if nargin < nreqarg + 1
    targinfo = S;
else
    targinfo = varargin{1};
end

if nargin < nreqarg + 2
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

if nargin < nreqarg + 3
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
vn = (-1)^io.*normal(dom);

if abs(zk) < eps
    % n x grad S_0[sigma]
    gradS0sigma = taylor.static.eval_gradS0(S,sigmavals,epstaylor, ...
        targinfo,opts);
    gradS0sigma = array_to_surfacefun(gradS0sigma.',dom,S);
    nxgradS0sigma = cross(vn,gradS0sigma);
    nxvals = surfacefun_to_array(nxgradS0sigma,dom,S);
    
    % S_0[n x grad S_0[sigma]]
    S0nx1 = lap3d.dirichlet.eval(S,real(nxvals(:,1)),targinfo,epslh, ...
        dpars,optslh);
    S0nx1 = S0nx1 + ...
        1i*lap3d.dirichlet.eval(S,imag(nxvals(:,1)),targinfo,epslh, ...
        dpars,optslh);
    S0nx2 = lap3d.dirichlet.eval(S,real(nxvals(:,2)),targinfo,epslh, ...
        dpars,optslh);
    S0nx2 = S0nx2 + ...
        1i*lap3d.dirichlet.eval(S,imag(nxvals(:,2)),targinfo,epslh, ...
        dpars,optslh);
    S0nx3 = lap3d.dirichlet.eval(S,real(nxvals(:,3)),targinfo,epslh, ...
        dpars,optslh);
    S0nx3 = S0nx3 + ...
        1i*lap3d.dirichlet.eval(S,imag(nxvals(:,3)),targinfo,epslh, ...
        dpars,optslh);
    S0nx = [S0nx1 S0nx2 S0nx3];
    S0nx = array_to_surfacefun(S0nx,dom,S);
    
    if aint
        fluxsigma = -TaylorState.intacyc(S0nx,n,nv);
    else
        fluxsigma = -TaylorState.intbcyc(S0nx,n,nu);
    end
else
    m0 = TaylorState.debyem0(sigma,zk,L,vn);
    m0vals = surfacefun_to_array(m0,dom,S);

    % S_k[m_0]
    Skm01 = helm3d.dirichlet.eval(S,m0vals(:,1),targinfo,epslh,zk, ...
        dpars,optslh);
    Skm02 = helm3d.dirichlet.eval(S,m0vals(:,2),targinfo,epslh,zk, ...
        dpars,optslh);
    Skm03 = helm3d.dirichlet.eval(S,m0vals(:,3),targinfo,epslh,zk, ...
        dpars,optslh);
    Skm0 = [Skm01 Skm02 Skm03];
    Skm0 = array_to_surfacefun(Skm0,dom,S);

    % curl S_k[m_0]
    m0vals = m0vals.';
    curlSkm0 = taylor.dynamic.eval_curlSk(S,zk,m0vals,epstaylor, ...
        targinfo,opts);
    curlSkm0 = array_to_surfacefun(curlSkm0.',dom,S);
    
    % A-cycle integral
    integrand = 1i.*Skm0 + (m0./2 + 1i.*curlSkm0)./zk;
    if aint
        fluxsigma = TaylorState.intacyc(integrand,n,nv);
    else
        fluxsigma = TaylorState.intbcyc(integrand,n,nu);
    end
end

end