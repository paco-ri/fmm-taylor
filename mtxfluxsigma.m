function fluxsigma = mtxfluxsigma(S,dom,n,nu,nv,sigma,zk,eps,varargin)
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
%     * zk: [double complex] wavenumber
%     * eps: [double] precision requested
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

dpars = [1.0, 0.0];
if nargin < 9
    targinfo = S;
else
    targinfo = varargin{1};
    if nargin < 10
        if abs(zk) < 1e-16
            Q = taylor.static.get_quadrature_correction(S,eps, ...
                targinfo,opts);
        else
            Q = taylor.dynamic.get_quadrature_correction(S,zk,eps, ...
                targinfo,opts);
        end
        opts = [];
        opts.format = 'rsc';
        opts.precomp_quadrature = Q;
    else
        opts = varargin{2};
        if nargin < 11
            if abs(zk) < 1e-16
                Qlh = lap3d.dirichlet.get_quadrature_correction(S,eps, ...
                    dpars,targinfo,opts);
            else
                Qlh = helm3d.dirichlet.get_quadrature_correction(S, ...
                    eps,zk,dpars,targinfo,opts);
            end
            optslh = [];
            optslh.format = 'rsc';
            optslh.precomp_quadrature = Qlh;
        else
            optslh = varargin{3};
        end
    end
end

sigmavals = surfacefun_to_array(sigma,dom,S);
sigmavals = sigmavals.';
vn = normal(dom);

if abs(zk) < 1e-16
    % n x grad S_0[sigma]
    gradS0sigma = taylor.static.eval_gradS0(S,sigmavals,eps,targinfo,opts);
    gradS0sigma = array_to_surfacefun(gradS0sigma.',dom,S);
    nxgradS0sigma = cross(vn,gradS0sigma);
    nxvals = surfacefun_to_array(nxgradS0sigma,dom,S);
    
    % S_0[n x grad S_0[sigma]]
    S0nx1 = lap3d.dirichlet.eval(S,nxvals(:,1),targinfo,eps,dpars,optslh);
    S0nx2 = lap3d.dirichlet.eval(S,nxvals(:,2),targinfo,eps,dpars,optslh);
    S0nx3 = lap3d.dirichlet.eval(S,nxvals(:,3),targinfo,eps,dpars,optslh);
    S0nx = [S0nx1 S0nx2 S0nx3];
    S0nx = array_to_surfacefun(S0nx,dom,S);
    
    fluxsigma = -intacyc(S0nx,n,nu,nv);
else
    m0 = debyem0(sigma,zk);
    m0vals = surfacefun_to_array(m0,dom,S);

    % S_k[m_0]
    Skm01 = helm3d.dirichlet.eval(S,m0vals(:,1),targinfo,eps,zk, ...
        dpars,optslh);
    Skm02 = helm3d.dirichlet.eval(S,m0vals(:,2),targinfo,eps,zk, ...
        dpars,optslh);
    Skm03 = helm3d.dirichlet.eval(S,m0vals(:,3),targinfo,eps,zk, ...
        dpars,optslh);
    Skm0 = [Skm01 Skm02 Skm03];
    Skm0 = array_to_surfacefun(Skm0,dom,S);

    % curl S_k[m_0]
    m0vals = m0vals.';
    curlSkm0 = taylor.dynamic.eval_curlSk(S,zk,m0vals,eps,targinfo,opts);
    curlSkm0 = array_to_surfacefun(curlSkm0.',dom,S);
    
    % A-cycle integral
    integrand = 1i.*Skm0 + (m0./2 + 1i.*curlSkm0)./zk;
    fluxsigma = intacyc(integrand,n,nu,nv);
end

end