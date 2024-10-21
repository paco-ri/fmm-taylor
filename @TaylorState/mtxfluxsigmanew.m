function fluxsigma = mtxfluxsigmanew(S,dom,vn,L,domparams,sigma,zk,epstaylor,epslh,varargin)
%MTXFLUXSIGMA compute sigma-dep. terms of flux condition
% 
%   Required arguments:
%     * S: cell of surfer objects 
%     * dom: surfacemesh version of S
%     * vn: normal vectors on S
%     * L: surfaceops on S
%     * domparams: parameters describing dom and circulation [n, nu, nv, io, aint]
%         n: [int] polynomial order on each surface patch
%         nu: [int] number of patches in toroidal direction
%         nv: [int] number of patches in poloidal direction
%     * sigma: [cell of surfacefuns] densities for which 
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

dpars = [1.0, 0.0];
nreqarg = 9;

if length(dom) == 1

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
        
        fluxsigma = -TaylorState.intacyc(S0nx,n,nv);
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
        fluxsigma = TaylorState.intacyc(integrand,n,nv);
    end

% toroidal shell case
else

    So = S{1};
    Si = S{2};
    domo = dom{1};
    domi = dom{2};

    if nargin < nreqarg + 1
        targinfoo = S{1};
        targinfoi = S{2};
    else
        targinfoo = varargin{1}{1};
        targinfoi = varargin{1}{2};
    end
    
    if nargin < nreqarg + 2
        if abs(zk) < eps
            Qo = taylor.static.get_quadrature_correction(So, ...
                epstaylor,targinfoo);
            Qi = taylor.static.get_quadrature_correction(Si, ...
                epstaylor,targinfoi);
        else
            Qo = taylor.dynamic.get_quadrature_correction(So,zk, ...
                epstaylor,targinfoo);
            Qi = taylor.dynamic.get_quadrature_correction(Si,zk, ...
                epstaylor,targinfoi);
        end
        optso = [];
        optso.format = 'rsc';
        optso.precomp_quadrature = Qo;
        optsi = [];
        optsi.format = 'rsc';
        optsi.precomp_quadrature = Qi;
    else
        optso = varargin{2}{1};
        optsi = varargin{2}{2};
    end
    
    if nargin < nreqarg + 3
        if abs(zk) < eps
            Qlho = lap3d.dirichlet.get_quadrature_correction(So, ...
                epslh,dpars,targinfoo);
            Qlhi = lap3d.dirichlet.get_quadrature_correction(Si, ...
                epslh,dpars,targinfoi);
        else
            Qlho = helm3d.dirichlet.get_quadrature_correction(So, ...
                epslh,zk,dpars,targinfoo);
            Qlhi = helm3d.dirichlet.get_quadrature_correction(Si, ...
                epslh,zk,dpars,targinfoi);
        end
        optslho = [];
        optslho.format = 'rsc';
        optslho.precomp_quadrature = Qlho;
        optslhi = [];
        optslhi.format = 'rsc';
        optslhi.precomp_quadrature = Qlhi;
    else
        optslho = varargin{3}{1};
        optslhi = varargin{3}{2};
    end

    % ====
    
    sigmavalso = surfacefun_to_array(sigma{1},domo,So);
    sigmavalso = sigmavalso.';
    sigmavalsi = surfacefun_to_array(sigma{2},domi,Si);
    sigmavalsi = sigmavalsi.';
    vno = vn{1};
    vni = vn{2};

    if abs(zk) < eps
        % n x grad S_0[sigma]
        gradS0sigmao = taylor.static.eval_gradS0(So,sigmavalso, ...
            epstaylor,targinfoo,optso);
        gradS0sigmao = array_to_surfacefun(gradS0sigmao.',domo,So);
        nxgradS0sigmao = cross(vno,gradS0sigmao);
        nxvalso = surfacefun_to_array(nxgradS0sigmao,domo,So);

        gradS0sigmai = taylor.static.eval_gradS0(Si,sigmavalsi, ...
            epstaylor,targinfoi,optsi);
        gradS0sigmai = array_to_surfacefun(gradS0sigmai.',domi,Si);
        nxgradS0sigmai = cross(vni,gradS0sigmai);
        nxvalsi = surfacefun_to_array(nxgradS0sigmai,domi,Si);

        gradS0sigmao2i = taylor.static.eval_gradS0(So,sigmavalso, ...
            epstaylor,targinfoi);
        gradS0sigmao2i = array_to_surfacefun(gradS0sigmao2i.',domi,Si);
        nxgradS0sigmao2i = cross(vni,gradS0sigmao2i);
        nxvalso2i = surfacefun_to_array(nxgradS0sigmao2i,domi,Si);

        gradS0sigmai2o = taylor.static.eval_gradS0(Si,sigmavalsi, ...
            epstaylor,targinfoo);
        gradS0sigmai2o = array_to_surfacefun(gradS0sigmai2o.',domo,So);
        nxgradS0sigmai2o = cross(vno,gradS0sigmai2o);
        nxvalsi2o = surfacefun_to_array(nxgradS0sigmai2o,domo,So);
        
        % S_0[n x grad S_0[sigma]]
        S0nxo = taylor.helper.lap_dir_vec_eval(So,nxvalso.',targinfoo, ...
            epslh,dpars,optslho);
        S0nxo = array_to_surfacefun(S0nxo.',domo,So);
        S0nxi = taylor.helper.lap_dir_vec_eval(Si,nxvalsi.',targinfoi, ...
            epslh,dpars,optslhi);
        S0nxi = array_to_surfacefun(S0nxi.',domi,Si);

        % 18 Oct 2024 trying "double coupling"
        S0nxo2i2i = taylor.helper.lap_dir_vec_eval(Si,nxvalso2i.',targinfoi, ...
            epslh,dpars,optslhi);
        S0nxo2i2i = array_to_surfacefun(S0nxo2i2i.',domi,Si);
        S0nxi2o2o = taylor.helper.lap_dir_vec_eval(So,nxvalsi2o.',targinfoo, ...
            epslh,dpars,optslho);
        S0nxi2o2o = array_to_surfacefun(S0nxi2o2o.',domo,So);

        S0nxo2o2i = taylor.helper.lap_dir_vec_eval(So,nxvalso.',targinfoi, ...
            epslh,dpars);
        S0nxo2o2i = array_to_surfacefun(S0nxo2o2i.',domi,Si);
        S0nxi2i2o = taylor.helper.lap_dir_vec_eval(Si,nxvalsi.',targinfoo, ...
            epslh,dpars);
        S0nxi2i2o = array_to_surfacefun(S0nxi2i2o.',domo,So);

        S0nxo2i2o = taylor.helper.lap_dir_vec_eval(Si,nxvalso2i.',targinfoo, ...
            epslh,dpars);
        S0nxo2i2o = array_to_surfacefun(S0nxo2i2o.',domo,So);
        S0nxi2o2i = taylor.helper.lap_dir_vec_eval(So,nxvalsi2o.',targinfoi, ...
            epslh,dpars);
        S0nxi2o2i = array_to_surfacefun(S0nxi2o2i.',domi,Si);

        fluxsigma = zeros(2);
        % note negative sign
        fluxsigma(1,1) = -TaylorState.intacyc(S0nxo,n,nv) ...
            + TaylorState.intacyc(S0nxo2i2i,n,nv) ...
            - TaylorState.intacyc(S0nxo2i2o,n,nv) ...
            + TaylorState.intacyc(S0nxo2o2i,n,nv);
        fluxsigma(1,2) = -TaylorState.intacyc(S0nxi2o2o,n,nv) ...
            + TaylorState.intacyc(S0nxi,n,nv) ...
            - TaylorState.intacyc(S0nxi2i2o,n,nv) ...
            + TaylorState.intacyc(S0nxi2o2i,n,nv);
        fluxsigma(2,1) = -TaylorState.intbcyc(S0nxo,n,nu) ...
            + TaylorState.intbcyc(S0nxo2i2i,n,nu) ...
            - TaylorState.intbcyc(S0nxo2i2o,n,nu) ...
            + TaylorState.intbcyc(S0nxo2o2i,n,nu);
        fluxsigma(2,2) = -TaylorState.intbcyc(S0nxi2o2o,n,nu) ...
            + TaylorState.intbcyc(S0nxi,n,nu) ...
            - TaylorState.intbcyc(S0nxi2i2o,n,nu) ...
            + TaylorState.intbcyc(S0nxi2o2i,n,nu);
        % 17 Oct 2024 testing alternate splitting
        % fluxsigma(1,1) = -TaylorState.intacyc(S0nxo,n,nv) ...
        %     - TaylorState.intacyc(S0nxi2o,n,nv);
        % fluxsigma(1,2) = TaylorState.intacyc(S0nxo2i,n,nv) ...
        %     + TaylorState.intacyc(S0nxi,n,nv);
        % fluxsigma(2,1) = -TaylorState.intbcyc(S0nxo,n,nu) ...
        %     - TaylorState.intbcyc(S0nxi2o,n,nu);
        % fluxsigma(2,2) = TaylorState.intbcyc(S0nxo2i,n,nu) ...
        %     + TaylorState.intbcyc(S0nxi,n,nu);
    else
        m0o = TaylorState.debyem0(sigma{1},zk,L{1},vn{1});
        m0i = TaylorState.debyem0(sigma{2},zk,L{2},vn{2});
        m0ovals = surfacefun_to_array(m0o,dom{1},S{1});
        m0ivals = surfacefun_to_array(m0i,dom{2},S{2});

        % S_k[m_0]
        Skm0o = taylor.helper.helm_dir_vec_eval(S{1},m0ovals.',targinfoo, ...
            epslh,zk,dpars,optslho);
        Skm0i = taylor.helper.helm_dir_vec_eval(S{2},m0ivals.',targinfoi, ...
            epslh,zk,dpars,optslhi); 
        Skm0o2i = taylor.helper.helm_dir_vec_eval(S{1},m0ovals.', ...
            targinfoi,epslh,zk,dpars);
        Skm0i2o = taylor.helper.helm_dir_vec_eval(S{2},m0ivals.', ...
            targinfoo,epslh,zk,dpars);

        Skm0o = array_to_surfacefun(Skm0o.',dom{1},S{1});
        Skm0i = array_to_surfacefun(Skm0i.',dom{2},S{2});
        Skm0o2i = array_to_surfacefun(Skm0o2i.',dom{2},S{2});
        Skm0i2o = array_to_surfacefun(Skm0i2o.',dom{1},S{1});

        % curl S_k[m_0]
        curlSkm0o = taylor.dynamic.eval_curlSk(S{1},zk,m0ovals.', ...
            epstaylor,targinfoo,optso);
        curlSkm0i = taylor.dynamic.eval_curlSk(S{2},zk,m0ivals.', ...
            epstaylor,targinfoi,optsi);
        curlSkm0o2i = taylor.dynamic.eval_curlSk(S{1},zk,m0ovals.', ...
            epstaylor,targinfoi);
        curlSkm0i2o = taylor.dynamic.eval_curlSk(S{2},zk,m0ivals.', ...
            epstaylor,targinfoo);

        curlSkm0o = array_to_surfacefun(curlSkm0o.',dom{1},S{1});
        curlSkm0i = array_to_surfacefun(curlSkm0i.',dom{2},S{2});
        curlSkm0o2i = array_to_surfacefun(curlSkm0o2i.',dom{2},S{2});
        curlSkm0i2o = array_to_surfacefun(curlSkm0i2o.',dom{1},S{1});

        fluxsigma = zeros(2);
        % note negative sign
        % fluxsigma(1,1) = TaylorState.intacyc(...
        %     1i.*Skm0o + (m0o./2 + 1i.*curlSkm0o)./zk,n,nv) ...
        %     - TaylorState.intacyc(1i.*Skm0o2i + 1i.*curlSkm0o2i./zk,n,nv);
        % fluxsigma(1,2) = TaylorState.intacyc(...
        %     1i.*Skm0i2o + 1i.*curlSkm0i2o./zk,n,nv) ...
        %     - TaylorState.intacyc(...
        %     1i.*Skm0i + (m0i./2 + 1i.*curlSkm0i)./zk,n,nv);
        % fluxsigma(2,1) = TaylorState.intbcyc(...
        %     1i.*Skm0o + (m0o./2 + 1i.*curlSkm0o)./zk,n,nu) ...
        %     - TaylorState.intbcyc(1i.*Skm0o2i + 1i.*curlSkm0o2i./zk,n,nu);
        % fluxsigma(2,2) = TaylorState.intbcyc(...
        %     1i.*Skm0i2o + 1i.*curlSkm0i2o./zk,n,nu) ...
        %     - TaylorState.intbcyc(...
        %     1i.*Skm0i + (m0i./2 + 1i.*curlSkm0i)./zk,n,nu);
        fluxsigma(1,1) = TaylorState.intacyc(...
            1i.*Skm0o + (m0o./2 + 1i.*curlSkm0o)./zk,n,nv) ...
            + TaylorState.intacyc(...
            1i.*Skm0i2o + 1i.*curlSkm0i2o./zk,n,nv);
        fluxsigma(1,2) = -TaylorState.intacyc(...
            1i.*Skm0o2i + 1i.*curlSkm0o2i./zk,n,nv)...
            - TaylorState.intacyc(...
            1i.*Skm0i + (m0i./2 + 1i.*curlSkm0i)./zk,n,nv);
        fluxsigma(2,1) = TaylorState.intbcyc(...
            1i.*Skm0o + (m0o./2 + 1i.*curlSkm0o)./zk,n,nu) ...
            + TaylorState.intbcyc(...
            1i.*Skm0i2o + 1i.*curlSkm0i2o./zk,n,nu);
        fluxsigma(2,2) = - TaylorState.intbcyc(...
            1i.*Skm0o2i + 1i.*curlSkm0o2i./zk,n,nu) ...
            - TaylorState.intbcyc(...
            1i.*Skm0i + (m0i./2 + 1i.*curlSkm0i)./zk,n,nu);
    end
    
end

end