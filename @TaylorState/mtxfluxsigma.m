function fluxsigma = mtxfluxsigma(domain,sigma,zk,epstaylor,epslh,varargin)
%MTXFLUXSIGMA compute sigma-dep. terms of flux condition
% 
%   Required arguments:
%     * domain: Domain object
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

S = domain.surf;
dom = domain.dom;
vn = domain.vn;
L = domain.L;
domparams = domain.domparams;

n = domparams(1);
nu = domparams(2);
nv = domparams(3);

dpars = [1.0, 0.0];
nreqarg = 5;

if isa(dom,'surfacemesh')

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

    if nargin < nreqarg + 1
        targinfoo = S{1};
        targinfoi = S{2};
        targinfo = S;
    else
        targinfoo = varargin{1}{1};
        targinfoi = varargin{1}{2};
        targinfo = varargin{1};
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
        opts = {optso,optsi};
    else
        opts = varargin{2};
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
        optslh = {optslho,optslhi};
    else
        optslh = varargin{3};
    end

    % ====
    
    sig = cell(2);
    for i = 1:2
        for j = 1:2
            sig{i,j} = surfacefun_to_array(sigma{i,j},dom{i},S{i});
        end
    end
    % i: inner or outer surface
    % j: RHS was mH1 or mH2

    if abs(zk) < eps
        % n x grad S_0[sigma]
        nxgradS0 = cell(2,2,2);
        for i = 1:2
            for j = 1:2
                for k = 1:2
                    gradS0sig = taylor.static.eval_gradS0(S{i}, ...
                        sig{i,j},epstaylor,S{k},opts{i,k});
                    gradS0sig = array_to_surfacefun(gradS0sig.', ...
                        dom{k},S{k});
                    nxgradS0{i,j,k} = cross(vn{k},gradS0sig);
                    nxgradS0{i,j,k} = surfacefun_to_array(...
                        nxgradS0{i,j,k},dom{k},S{k});
                end
            end
        end
        % k: n x grad S0[sigma] target inner or outer surface
        
        % S_0[n x grad S_0[sigma]]
        S0nx = cell(2,2,2,2);
        for i = 1:2
            for j = 1:2
                for k = 1:2
                    for l = 1:2
                        S0nx{i,j,k,l} = ...
                            taylor.helper.lap_dir_vec_eval( ...
                            S{k},nxgradS0{i,j,k}.',S{l},epslh, ...
                            dpars,optslh{k,l});
                        S0nx{i,j,k,l} = array_to_surfacefun( ...
                            S0nx{i,j,k,l}.',dom{l},S{l});
                    end
                end
            end
        end
        % l: S0[n x grad S_0[sigma]] target points outer or inner surf

        fluxsigma = zeros(2);
        % rows: toroidal or poloidal flux
        % j: as before
        for j = 1:2
            for i = 1:2
                for k = 1:2
                    for l = 1:2
                        fluxsigma(1,j) = fluxsigma(1,j) ...
                            + (-1)^l.*TaylorState.intacyc( ...
                            S0nx{i,j,k,l},n,nv);
                        fluxsigma(2,j) = fluxsigma(2,j) ...
                            + (-1)^l.*TaylorState.intbcyc( ...
                            S0nx{i,j,k,l},n,nu);
                    end
                end
            end
        end
    else
        % i: inner or outer surface
        % j: RHS was mH1 or mH2
        m0 = cell(2);
        m0vals = cell(2);
        for i = 1:2
            for j = 1:2
                m0{i,j} = TaylorState.debyem0(sigma{i,j},zk,L{i},vn{i});
                m0vals{i,j} = surfacefun_to_array(m0{i,j},dom{i},S{i});
            end
        end

        % S_k[m_0]
        Skm0 = cell(2,2,2);
        for i = 1:2
            for j = 1:2
                for k = 1:2
                    temp = taylor.helper.helm_dir_vec_eval( ...
                        S{i},m0vals{i,j}.',targinfo{k},epslh,zk, ...
                        dpars,optslh{i,k});
                    Skm0{i,j,k} = array_to_surfacefun(temp.', ...
                        dom{k},S{k});
                end
            end
        end
        % k: Sk[m0] target inner or outer surface

        % curl S_k[m_0]
        curlSkm0 = cell(2,2,2);
        for i = 1:2
            for j = 1:2
                for k = 1:2
                    temp = taylor.dynamic.eval_curlSk(S{i},zk, ...
                        m0vals{i,j}.',epstaylor,targinfo{k},opts{i,k});
                    curlSkm0{i,j,k} = array_to_surfacefun(temp.', ...
                        dom{k},S{k});
                end
            end
        end
        % k: curl Sk[m0] target inner or outer surface

        fluxsigma = zeros(2);
        % rows: toroidal or poloidal flux
        % j: as before
        for j = 1:2
            for i = 1:2
                for k = 1:2
                    if i == k
                        fluxsigma(1,j) = fluxsigma(1,j) ...
                            + (-1)^k.*TaylorState.intacyc( ...
                            1i.*Skm0{i,j,k} + (m0{i,j}./2 ...
                            + 1i.*curlSkm0{i,j,k})./zk,n,nv);
                        fluxsigma(2,j) = fluxsigma(2,j) ...
                            + (-1)^k.*TaylorState.intbcyc( ...
                            1i.*Skm0{i,j,k} + (m0{i,j}./2 ...
                            + 1i.*curlSkm0{i,j,k})./zk,n,nu);
                    else
                       fluxsigma(1,j) = fluxsigma(1,j) ...
                            + (-1)^k.*TaylorState.intacyc( ...
                            1i.*(Skm0{i,j,k} + curlSkm0{i,j,k}./zk),n,nv);
                       fluxsigma(2,j) = fluxsigma(2,j) ...
                            + (-1)^k.*TaylorState.intbcyc( ...
                            1i.*(Skm0{i,j,k} + curlSkm0{i,j,k}./zk),n,nu); 
                    end
                end
            end
        end
    end
    
end

end