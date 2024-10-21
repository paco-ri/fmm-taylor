function fluxalpha = mtxfluxalphanew(S,dom,vn,domparams,mH,zk,epstaylor,epslh,varargin)
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

if length(dom) == 1
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
        nxminusmH = nxcurlS0mH - bd.*mH./2;
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
    
    mHo = mH{1};
    mHi = mH{2};
    mHvalso = surfacefun_to_array(mHo,domo,So);
    mHvalso = mHvalso.';
    mHvalsi = surfacefun_to_array(mHi,domi,Si);
    mHvalsi = mHvalsi.';
    vno = vn{1};
    vni = vn{2};
    
    if abs(zk) < eps
        % n x curl S0[mH] - mH/2
        curlS0mHo = taylor.static.eval_curlS0(So,mHvalso,epstaylor, ...
            targinfoo,optso);
        curlS0mHo = array_to_surfacefun(curlS0mHo.',domo,So);
        nxcurlS0mHo = cross(vno,curlS0mHo);

        curlS0mHi = taylor.static.eval_curlS0(Si,mHvalsi,epstaylor, ...
            targinfoi,optsi);
        curlS0mHi = array_to_surfacefun(curlS0mHi.',domi,Si);
        nxcurlS0mHi = cross(vni,curlS0mHi);

        curlS0mHo2i = taylor.static.eval_curlS0(So,mHvalso,epstaylor, ...
            targinfoi);
        curlS0mHo2i = array_to_surfacefun(curlS0mHo2i.',domi,Si);
        nxcurlS0mHo2i = cross(vni,curlS0mHo2i);

        curlS0mHi2o = taylor.static.eval_curlS0(Si,mHvalsi,epstaylor, ...
            targinfoo);
        curlS0mHi2o = array_to_surfacefun(curlS0mHi2o.',domo,So);
        nxcurlS0mHi2o = cross(vno,curlS0mHi2o);

        % boundary terms for on-surface potentials
        nxminusmHo = nxcurlS0mHo - mHo./2;
        nxminusmHi = nxcurlS0mHi - mHi./2;

        % surfacefuns to arrays
        nxo = surfacefun_to_array(nxminusmHo,domo,So);
        nxi = surfacefun_to_array(nxminusmHi,domi,Si);
        nxo2i = surfacefun_to_array(nxcurlS0mHo2i,domi,Si);
        nxi2o = surfacefun_to_array(nxcurlS0mHi2o,domo,So);
        
        % S0[n x curl S0[mH] - mH/2]
        S0nxo = taylor.helper.lap_dir_vec_eval(So,nxo.',targinfoo,epslh, ...
            dpars,optslho);
        S0nxi = taylor.helper.lap_dir_vec_eval(Si,nxi.',targinfoi,epslh, ...
            dpars,optslhi);

        % 18 Oct 2024 trying "double coupling"
        S0nxo2i2i = taylor.helper.lap_dir_vec_eval(Si,nxo2i.',targinfoi,epslh, ...
            dpars,optslhi);
        S0nxi2o2o = taylor.helper.lap_dir_vec_eval(So,nxi2o.',targinfoo,epslh, ...
            dpars,optslho);

        S0nxo2o2i = taylor.helper.lap_dir_vec_eval(So,nxo.',targinfoi,epslh, ...
            dpars);
        S0nxi2i2o = taylor.helper.lap_dir_vec_eval(Si,nxi.',targinfoo,epslh, ...
            dpars);

        S0nxo2i2o = taylor.helper.lap_dir_vec_eval(Si,nxo2i.',targinfoo,epslh, ...
            dpars);
        S0nxi2o2i = taylor.helper.lap_dir_vec_eval(So,nxi2o.',targinfoi,epslh, ...
            dpars);

        S0nxo = array_to_surfacefun(S0nxo.',domo,So);
        S0nxi = array_to_surfacefun(S0nxi.',domi,Si);
        S0nxo2i2i = array_to_surfacefun(S0nxo2i2i.',domi,Si);
        S0nxi2o2o = array_to_surfacefun(S0nxi2o2o.',domo,So);
        S0nxo2o2i = array_to_surfacefun(S0nxo2o2i.',domi,Si);
        S0nxi2i2o = array_to_surfacefun(S0nxi2i2o.',domo,So);
        S0nxo2i2o = array_to_surfacefun(S0nxo2i2o.',domo,So);
        S0nxi2o2i = array_to_surfacefun(S0nxi2o2i.',domi,Si);
    
        fluxalpha = zeros(2);
        fluxalpha(1,1) = 1i*(TaylorState.intacyc(S0nxo,n,nv) ...
            - TaylorState.intacyc(S0nxo2i2i,n,nv) ...
            + TaylorState.intacyc(S0nxo2i2o,n,nv) ...
            - TaylorState.intacyc(S0nxo2o2i,n,nv));
        fluxalpha(1,2) = 1i*(TaylorState.intacyc(S0nxi2o2o,n,nv) ...
            - TaylorState.intacyc(S0nxi,n,nv) ...
            + TaylorState.intacyc(S0nxi2i2o,n,nv) ...
            - TaylorState.intacyc(S0nxi2o2i,n,nv));
        fluxalpha(2,1) = 1i*(TaylorState.intbcyc(S0nxo,n,nu) ...
            - TaylorState.intbcyc(S0nxo2i2i,n,nu) ...
            + TaylorState.intbcyc(S0nxo2i2o,n,nu) ...
            - TaylorState.intbcyc(S0nxo2o2i,n,nu));
        fluxalpha(2,2) = 1i*(TaylorState.intbcyc(S0nxi2o2o,n,nu) ...
            - TaylorState.intbcyc(S0nxi,n,nu) ...
            + TaylorState.intbcyc(S0nxi2i2o,n,nu) ...
            - TaylorState.intbcyc(S0nxi2o2i,n,nu));
        % 18 Oct 2024 trying "double coupling"
        % fluxalpha(1,1) = 1i*(TaylorState.intacyc(S0nxo,n,nv) ...
        %     + TaylorState.intacyc(S0nxi2o2o,n,nv));
        % fluxalpha(1,2) = 1i*(-TaylorState.intacyc(S0nxo2i2i,n,nv) ...
        %     - TaylorState.intacyc(S0nxi,n,nv));
        % fluxalpha(2,1) = 1i*(TaylorState.intbcyc(S0nxo,n,nu) ...
        %     + TaylorState.intbcyc(S0nxi2o2o,n,nu));
        % fluxalpha(2,2) = 1i*(-TaylorState.intbcyc(S0nxo2i2i,n,nu) ...
        %     - TaylorState.intbcyc(S0nxi,n,nu));
    else
        % Sk[mH]
        SkmHo = taylor.helper.helm_dir_vec_eval(So,mHvalso, ...
            targinfoo,epslh,zk,dpars,optslho);
        SkmHi = taylor.helper.helm_dir_vec_eval(Si,mHvalsi, ...
            targinfoi,epslh,zk,dpars,optslhi);
        SkmHo2i = taylor.helper.helm_dir_vec_eval(So,mHvalso, ...
            targinfoi,epslh,zk,dpars);
        SkmHi2o = taylor.helper.helm_dir_vec_eval(Si,mHvalsi, ...
            targinfoo,epslh,zk,dpars);

        SkmHo = array_to_surfacefun(SkmHo.',domo,So);
        SkmHi = array_to_surfacefun(SkmHi.',domi,Si);
        SkmHo2i = array_to_surfacefun(SkmHo2i.',domi,Si);
        SkmHi2o = array_to_surfacefun(SkmHi2o.',domo,So);

        % curl Sk[mH]
        curlSkmHo = taylor.dynamic.eval_curlSk(So,zk,mHvalso, ...
            epstaylor,targinfoo,optso);
        curlSkmHi = taylor.dynamic.eval_curlSk(Si,zk,mHvalsi, ...
            epstaylor,targinfoi,optsi);
        curlSkmHo2i = taylor.dynamic.eval_curlSk(So,zk,mHvalso, ...
            epstaylor,targinfoi);
        curlSkmHi2o = taylor.dynamic.eval_curlSk(Si,zk,mHvalsi, ...
            epstaylor,targinfoo);

        curlSkmHo = array_to_surfacefun(curlSkmHo.',domo,So);
        curlSkmHi = array_to_surfacefun(curlSkmHi.',domi,Si);
        curlSkmHo2i = array_to_surfacefun(curlSkmHo2i.',domi,Si);
        curlSkmHi2o = array_to_surfacefun(curlSkmHi2o.',domo,So);

        % curl S0[mH]
        % TODO request quadrature correction?
        curlS0mHo = taylor.static.eval_curlS0(So,mHvalso, ...
            epstaylor);
        curlS0mHi = taylor.static.eval_curlS0(Si,mHvalsi, ...
            epstaylor);
        curlS0mHo2i = taylor.static.eval_curlS0(So,mHvalso, ...
            epstaylor,targinfoi);
        curlS0mHi2o = taylor.static.eval_curlS0(So,mHvalsi, ...
            epstaylor,targinfoo);

        curlS0mHo = array_to_surfacefun(curlS0mHo.',domo,So);
        curlS0mHi = array_to_surfacefun(curlS0mHi.',domi,Si);
        curlS0mHo2i = array_to_surfacefun(curlS0mHo2i.',domi,Si);
        curlS0mHi2o = array_to_surfacefun(curlS0mHi2o.',domo,So);

        fluxalpha = zeros(2);
        % fluxalpha(1,1) = TaylorState.intacyc( ...
        %     1i.*(SkmHo + (curlSkmHo - curlS0mHo)./zk), n, nv) ...
        %     - TaylorState.intacyc( ...
        %     1i.*(SkmHo2i + (curlSkmHo2i - curlS0mHo2i)./zk), n, nv);
        % fluxalpha(1,2) = TaylorState.intacyc( ...
        %     1i.*(SkmHi2o + (curlSkmHi2o - curlS0mHi2o)./zk), n, nv) ...
        %     - TaylorState.intacyc( ...
        %     1i.*(SkmHi + (curlSkmHi - curlS0mHi)./zk), n, nv);
        % fluxalpha(2,1) = TaylorState.intbcyc( ...
        %     1i.*(SkmHo + (curlSkmHo - curlS0mHo)./zk), n, nu) ...
        %     - TaylorState.intbcyc( ...
        %     1i.*(SkmHo2i + (curlSkmHo2i - curlS0mHo2i)./zk), n, nu);
        % fluxalpha(2,2) = TaylorState.intbcyc( ...
        %     1i.*(SkmHi2o + (curlSkmHi2o - curlS0mHi2o)./zk), n, nu) ...
        %     - TaylorState.intbcyc( ...
        %     1i.*(SkmHi + (curlSkmHi - curlS0mHi)./zk), n, nu);
        fluxalpha(1,1) = TaylorState.intacyc( ...
            1i.*(SkmHo + (curlSkmHo - curlS0mHo)./zk), n, nv) ...
            + TaylorState.intacyc( ...
            1i.*(SkmHi2o + (curlSkmHi2o - curlS0mHi2o)./zk), n, nv);
        fluxalpha(1,2) = -TaylorState.intacyc( ...
            1i.*(SkmHo2i + (curlSkmHo2i - curlS0mHo2i)./zk), n, nv) ...
            - TaylorState.intacyc( ...
            1i.*(SkmHi + (curlSkmHi - curlS0mHi)./zk), n, nv);
        fluxalpha(2,1) = TaylorState.intbcyc( ...
            1i.*(SkmHo + (curlSkmHo - curlS0mHo)./zk), n, nu) ...
            + TaylorState.intbcyc( ...
            1i.*(SkmHi2o + (curlSkmHi2o - curlS0mHi2o)./zk), n, nu);
        fluxalpha(2,2) = -TaylorState.intbcyc( ...
            1i.*(SkmHo2i + (curlSkmHo2i - curlS0mHo2i)./zk), n, nu) ...
            - TaylorState.intbcyc( ...
            1i.*(SkmHi + (curlSkmHi - curlS0mHi)./zk), n, nu);
    end    

end

end