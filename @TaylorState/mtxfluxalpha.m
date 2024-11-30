function fluxalpha = mtxfluxalpha(S,dom,vn,domparams,mH,zk,epstaylor,epslh,varargin)
%MTXFLUXALPHA compute alpha coefficient in flux condition
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * vn: normal to surface
%     * domparams: parameters describing dom and circulation [n, nu, nv, io, aint]
%         n: [int] polynomial order on each surface patch
%         nu: [int] number of patches in toroidal direction
%         nv: [int] number of patches in poloidal direction
%         io: [int] if 1, negate vn because inner torus
%         aint: [int] if 1, do A-cyc. integral; otherwise, B-cyc.
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
nreqarg = 8;

if isa(dom,'surfacemesh')
    if nargin < nreqarg + 1
        targinfo = S;
    else
        targinfo = varargin{1};
    end
    
    if nargin < nreqarg + 2
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
    
    mHvals = surfacefun_to_array(mH,dom,S);
    mHvals = mHvals.';
    
    if abs(zk) < eps
        % n x curl S0[mH] - mH/2
        curlS0mH = taylor.static.eval_curlS0(S,mHvals,epstaylor,targinfo,opts);
        curlS0mH = array_to_surfacefun(curlS0mH.',dom,S);
        nxcurlS0mH = cross(vn,curlS0mH);
        nxminusmH = nxcurlS0mH - mH./2;
        nxminusmHvals = surfacefun_to_array(nxminusmH,dom,S);
        
        % S0[n x curl S0[mH] - mH/2]
        S0nx1 = lap3d.dirichlet.eval(S,real(nxminusmHvals(:,1)),targinfo,epslh, ...
            dpars,optslh);
        S0nx1 = S0nx1 + ...
            1i*lap3d.dirichlet.eval(S,imag(nxminusmHvals(:,1)),targinfo,epslh, ...
            dpars,optslh);
        S0nx2 = lap3d.dirichlet.eval(S,real(nxminusmHvals(:,2)),targinfo,epslh, ...
            dpars,optslh);
        S0nx2 = S0nx2 + ...
            1i*lap3d.dirichlet.eval(S,imag(nxminusmHvals(:,2)),targinfo,epslh, ...
            dpars,optslh);
        S0nx3 = lap3d.dirichlet.eval(S,real(nxminusmHvals(:,3)),targinfo,epslh, ...
            dpars,optslh);
        S0nx3 = S0nx3 + ...
            1i*lap3d.dirichlet.eval(S,imag(nxminusmHvals(:,3)),targinfo,epslh, ...
            dpars,optslh);
        S0nx = [S0nx1 S0nx2 S0nx3];
        S0nx = array_to_surfacefun(S0nx,dom,S);
    
        fluxalpha = 1i*TaylorState.intacyc(S0nx,n,nv);
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
    
        % A/B-cycle integral
        integrand = 1i.*( SkmH + (curlSkmH - curlS0mH)./zk );
        fluxalpha = TaylorState.intacyc(integrand,n,nv);
    end

% toroidal shell case
else

    So = S{1};
    Si = S{2};
    domo = dom{1};
    domi = dom{2};

    if nargin < nreqarg + 1
        targinfo = S;
        targinfoo = S{1};
        targinfoi = S{2};
    else
        targinfo = varargin{1};
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
        opts = varargin{2};
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
        optslh = varargin{3};
        optslho = varargin{3}{1};
        optslhi = varargin{3}{2};
    end

    % ====
    
    mHvals = cell(1,2);
    for i = 1:2
        mHvals{i} = surfacefun_to_array(mH{i},dom{i},S{i});
        mHvals{i} = mHvals{i}.';
    end
    
    if abs(zk) < eps
        nxcurlS0mH = cell(2);
        % i: source is mHo or i
        % j: target is o or i
        for i = 1:2
            for j = 1:2
                if i == j
                    curlS0mH = taylor.static.eval_curlS0(S{i}, ...
                        mHvals{i},epstaylor,targinfo{i},opts{i});
                else
                    curlS0mH = taylor.static.eval_curlS0(S{i}, ...
                        mHvals{i},epstaylor,targinfo{j});
                end
                curlS0mH = array_to_surfacefun(curlS0mH.',dom{j},S{j});
                nxcurlS0mH{i,j} = cross(vn{j},curlS0mH);
                % diagonal terms for on-surface potentials
                if i == j
                    nxcurlS0mH{i,i} = nxcurlS0mH{i,i} - mH{i}./2;
                end
                nxcurlS0mH{i,j} = surfacefun_to_array(nxcurlS0mH{i,j}, ...
                    dom{j},S{j});
            end
        end

        % k: targ o or i
        S0nx = cell(2,2,2);
        for i = 1:2
            for j = 1:2
                for k = 1:2
                    if j == k
                        S0nx{i,j,k} = taylor.helper.lap_dir_vec_eval(S{j}, ...
                            nxcurlS0mH{i,j}.',targinfo{j},epslh,dpars,optslh{j});
                    else
                        S0nx{i,j,k} = taylor.helper.lap_dir_vec_eval(S{j}, ...
                            nxcurlS0mH{i,j}.',targinfo{k},epslh,dpars);
                    end
                    S0nx{i,j,k} = array_to_surfacefun(S0nx{i,j,k}.', ...
                        dom{k},S{k});
                end
            end
        end
        
        fluxalpha = zeros(2);
        for i = 1:2
            for j = 1:2
                for k = 1:2
                    fluxalpha(1,i) = fluxalpha(1,i) ...
                        + (-1)^(k-1)*1i.*TaylorState.intacyc(S0nx{i,j,k},n,nv);
                    fluxalpha(2,i) = fluxalpha(2,i) ...
                        + (-1)^(k-1)*1i.*TaylorState.intbcyc(S0nx{i,j,k},n,nu);
                end
            end
        end
    else
        % Sk[mH]
        % SkmHo = taylor.helper.helm_dir_vec_eval(So,mHvalso, ...
        %     targinfoo,epslh,zk,dpars,optslho);
        % SkmHi = taylor.helper.helm_dir_vec_eval(Si,mHvalsi, ...
        %     targinfoi,epslh,zk,dpars,optslhi);
        % SkmHo2i = taylor.helper.helm_dir_vec_eval(So,mHvalso, ...
        %     targinfoi,epslh,zk,dpars);
        % SkmHi2o = taylor.helper.helm_dir_vec_eval(Si,mHvalsi, ...
        %     targinfoo,epslh,zk,dpars);
        % 
        % SkmHo = array_to_surfacefun(SkmHo.',domo,So);
        % SkmHi = array_to_surfacefun(SkmHi.',domi,Si);
        % SkmHo2i = array_to_surfacefun(SkmHo2i.',domi,Si);
        % SkmHi2o = array_to_surfacefun(SkmHi2o.',domo,So);

        SkmH = cell(2);
        curlSkmH = cell(2);
        curlS0mH = cell(2);
        % i: source is mHo or i
        % j: target is o or i
        for i = 1:2
            for j = 1:2
                if i == j
                    SkmH{i,j} = taylor.helper.helm_dir_vec_eval(S{i}, ...
                        mHvals{i},targinfo{i},epslh,zk,dpars,optslh{i});
                    curlSkmH{i,j} = taylor.dynamic.eval_curlSk(S{i}, ...
                        zk,mHvals{i},epstaylor,targinfo{i},opts{i});
                else
                    SkmH{i,j} = taylor.helper.helm_dir_vec_eval(S{i}, ...
                        mHvals{i},targinfo{j},epslh,zk,dpars);
                    curlSkmH{i,j} = taylor.dynamic.eval_curlSk(S{i}, ...
                        zk,mHvals{i},epstaylor,targinfo{j},opts{j});
                end
                curlS0mH{i,j} = taylor.static.eval_curlS0(S{i}, ...
                    mHvals{i},epstaylor,targinfo{j});
                SkmH{i,j} = array_to_surfacefun(SkmH{i,j}.',dom{j},S{j});
                curlSkmH{i,j} = array_to_surfacefun(curlSkmH{i,j}.', ...
                    dom{j},S{j});
                curlS0mH{i,j} = array_to_surfacefun(curlS0mH{i,j}.', ...
                    dom{j},S{j});
            end
        end

        % curl Sk[mH]
        % curlSkmHo = taylor.dynamic.eval_curlSk(So,zk,mHvalso, ...
        %     epstaylor,targinfoo,optso);
        % curlSkmHi = taylor.dynamic.eval_curlSk(Si,zk,mHvalsi, ...
        %     epstaylor,targinfoi,optsi);
        % curlSkmHo2i = taylor.dynamic.eval_curlSk(So,zk,mHvalso, ...
        %     epstaylor,targinfoi);
        % curlSkmHi2o = taylor.dynamic.eval_curlSk(Si,zk,mHvalsi, ...
        %     epstaylor,targinfoo);
        % 
        % curlSkmHo = array_to_surfacefun(curlSkmHo.',domo,So);
        % curlSkmHi = array_to_surfacefun(curlSkmHi.',domi,Si);
        % curlSkmHo2i = array_to_surfacefun(curlSkmHo2i.',domi,Si);
        % curlSkmHi2o = array_to_surfacefun(curlSkmHi2o.',domo,So);

        % curl S0[mH]
        % TODO request quadrature correction?
        % curlS0mHo = taylor.static.eval_curlS0(So,mHvalso, ...
        %     epstaylor);
        % curlS0mHi = taylor.static.eval_curlS0(Si,mHvalsi, ...
        %     epstaylor);
        % curlS0mHo2i = taylor.static.eval_curlS0(So,mHvalso, ...
        %     epstaylor,targinfoi);
        % curlS0mHi2o = taylor.static.eval_curlS0(So,mHvalsi, ...
        %     epstaylor,targinfoo);
        % 
        % curlS0mHo = array_to_surfacefun(curlS0mHo.',domo,So);
        % curlS0mHi = array_to_surfacefun(curlS0mHi.',domi,Si);
        % curlS0mHo2i = array_to_surfacefun(curlS0mHo2i.',domi,Si);
        % curlS0mHi2o = array_to_surfacefun(curlS0mHi2o.',domo,So);

        fluxalpha = zeros(2);
        % fluxalpha(1,1) = TaylorState.intacyc( ...
        %     1i.*(SkmHo + (curlSkmHo - curlS0mHo)./zk), n, nv) ...
        %     + TaylorState.intacyc( ...
        %     1i.*(SkmHi2o + (curlSkmHi2o - curlS0mHi2o)./zk), n, nv);
        % fluxalpha(1,2) = -TaylorState.intacyc( ...
        %     1i.*(SkmHo2i + (curlSkmHo2i - curlS0mHo2i)./zk), n, nv) ...
        %     - TaylorState.intacyc( ...
        %     1i.*(SkmHi + (curlSkmHi - curlS0mHi)./zk), n, nv);
        % fluxalpha(2,1) = TaylorState.intbcyc( ...
        %     1i.*(SkmHo + (curlSkmHo - curlS0mHo)./zk), n, nu) ...
        %     + TaylorState.intbcyc( ...
        %     1i.*(SkmHi2o + (curlSkmHi2o - curlS0mHi2o)./zk), n, nu);
        % fluxalpha(2,2) = -TaylorState.intbcyc( ...
        %     1i.*(SkmHo2i + (curlSkmHo2i - curlS0mHo2i)./zk), n, nu) ...
        %     - TaylorState.intbcyc( ...
        %     1i.*(SkmHi + (curlSkmHi - curlS0mHi)./zk), n, nu);
        for i = 1:2
            for j = 1:2
                fluxalpha(1,i) = fluxalpha(1,i) ...
                    + (-1)^(j-1).*TaylorState.intacyc( ...
                    1i.*(SkmH{i} + (curlSkmH{i}-curlS0mH{i})./zk), n, nv);
                fluxalpha(2,i) = fluxalpha(2,i) ...
                    + (-1)^(j-1).*TaylorState.intbcyc( ...
                    1i.*(SkmH{i} + (curlSkmH{i}-curlS0mH{i})./zk), n, nu);
            end
        end
    end    

end

end