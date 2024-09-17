function Bsigma = mtxBsigma(S,dom,L,sigmavals,zk,epstaylor,epslh, ...
    varargin)
%MTXBSIGMA compute sigma-dep. terms of surface magnetic field
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * sigmavals: [double complex(*)] density for which 
%                  sigma/2 + n . grad S0[sigma]
%              is computed
%     * zk: [dcomplex] wavenumber
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

dpars = [1.0,0];

% torus case
if length(dom) == 1

    nreqarg = 7;
    if nargin < nreqarg + 1
        targinfo = S{1};
    else
        targinfo = varargin{1}{1};
    end

    if nargin < nreqarg + 2
        if abs(zk) < eps
            Q = taylor.static.get_quadrature_correction(S{1},epstaylor, ...
                targinfo);
        else
            Q = taylor.dynamic.get_quadrature_correction(S{1},zk, ...
                epstaylor,targinfo);
        end
        opts = [];
        opts.format = 'rsc';
        opts.precomp_quadrature = Q;
    else
        opts = varargin{2}{1};
    end

    if nargin < nreqarg + 3
        if abs(zk) < eps
            Qlh = lap3d.dirichlet.get_quadrature_correction(S{1}, ...
            epslh,dpars,targinfo);
            optslh = [];
            optslh.format = 'rsc';
            optslh.precomp_quadrature = Qlh;
        end
    else
        optslh = varargin{3}{1};
    end

    vn = normal(dom{1});
    sigma = array_to_surfacefun(sigmavals,dom{1},S{1});

    if abs(zk) < eps
        % compute n . grad Sk[sigma]
        gradSsigma = taylor.static.eval_gradS0(S{1},sigmavals, ...
            epstaylor,targinfo,opts);
        gradSsigma = array_to_surfacefun(gradSsigma.',dom{1},S{1}); 
        ngradSsigma = dot(vn,gradSsigma);

        Bsigma = {sigma./2 + ngradSsigma};
    else
        % compute m0
        m0 = TaylorState.debyem0(sigma,zk,L{1},vn);
        m0vals = surfacefun_to_array(m0,dom{1},S{1});

        % compute n . Sk[m0]
        Sm0 = complex(zeros(size(m0vals)));
        for j=1:3
            Sm0(:,j) = helm3d.dirichlet.eval(S{1},m0vals(:,j),targinfo, ...
                epslh,zk,[1.0 0],optslh);
        end
        Sm0 = array_to_surfacefun(Sm0,dom{1},S{1});

        % compute n . grad Sk[sigma] and n . curl Sk[m0]
        [gradSsigma, curlSm0] = taylor.dynamic.eval_gradcurlSk(S{1},zk, ...
            sigmavals,m0vals.',epstaylor,targinfo,opts);
        gradSsigma = array_to_surfacefun(gradSsigma.',dom{1},S{1}); 
        ngradSsigma = dot(vn,gradSsigma);
        curlSm0 = array_to_surfacefun(curlSm0.',dom{1},S{1});

        % combine
        % m0terms = 1i.*dot(n,zk.*Sm0+curlSm0);
        m0terms = 1i*zk.*dot(vn,Sm0) + 1i.*dot(vn,curlSm0);
        Bsigma = {sigma./2 + ngradSsigma - m0terms};
    end

% toroidal shell case    
else
    
    % targinfoo = outer surface, targinfoi = inner surface
    if nargin < 7
        targinfoo = S{1};
        targinfoi = S{2};
    else
        targinfoo = varargin{1}{1};
        targinfoi = varargin{1}{2};
    end

    % near quadrature corrections for +taylor routines
    if nargin < 8
        if abs(zk) < eps
            Qo = taylor.static.get_quadrature_correction(S{1}, ...
                epstaylor,targinfoo);
            Qi = taylor.static.get_quadrature_correction(S{2}, ...
                epstaylor,targinfoi);
        else
            Qo = taylor.dynamic.get_quadrature_correction(S{1},zk, ...
                epstaylor,targinfoo);
            Qi = taylor.dynamic.get_quadrature_correction(S{2},zk, ...
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

    % near quadrature corrections for Helmholtz layer potential
    % only needed if zk != 0 
    if nargin < 9
        if abs(zk) < eps
            Qlho = helm3d.dirichlet.get_quadrature_correction(S{1}, ...
                epslh,zk,dpars,targinfoo);
            Qlhi = helm3d.dirichlet.get_quadrature_correction(S{2}, ...
                epslh,zk,dpars,targinfoi);
            optslho = [];
            optslho.format = 'rsc';
            optslho.precomp_quadrature = Qlho;
            optslhi = [];
            optslhi.format = 'rsc';
            optslhi.precomp_quadrature = Qlhi;
        end
    else
        optslho = varargin{3}{1};
        optslhi = varargin{3}{2};
    end

    % split sigmavals into halves. top half is on outer surface, etc.
    npts = length(sigmavals)/2;
    sigmavalso = sigmavals(1:npts);
    sigmavalsi = sigmavals(npts+1:end);

    % evaluate layer potential
    if abs(zk) < eps
        gradSsigmao = taylor.static.eval_gradS0(S{1},sigmavalso, ...
            epstaylor,targinfoo,optso);
        gradSsigmai = taylor.static.eval_gradS0(S{2},sigmavalsi, ...
            epstaylor,targinfoi,optsi);
        gradSsigmai2o = taylor.static.eval_gradS0(S{2},sigmavalsi, ...
            epstaylor,targinfoo);
        gradSsigmao2i = taylor.static.eval_gradS0(S{1},sigmavalso, ...
            epstaylor,targinfoi);
    else
        gradSsigmao = taylor.dynamic.eval_gradSk(S{1},zk,sigmavals, ...
            epstaylor,targinfoo,optso);
        gradSsigmai = taylor.dynamic.eval_gradSk(S{2},zk,sigmavalsi, ...
            epstaylor,targinfoi,optsi);
        gradSsigmai2o = taylor.dynamic.eval_gradSk(S{2},zk,sigmavalsi, ...
            epstaylor,targinfoo);
        gradSsigmao2i = taylor.dynamic.eval_gradSk(S{1},zk,sigmavalso, ...
            epstaylor,targinfoi);
    end
    gradSsigmao = array_to_surfacefun(gradSsigmao.',dom{1},S{1}); % note transpose 
    gradSsigmai = array_to_surfacefun(gradSsigmai.',dom{2},S{2});
    gradSsigmai2o = array_to_surfacefun(gradSsigmai2o.',dom{1},S{1});
    gradSsigmao2i = array_to_surfacefun(gradSsigmao2i.',dom{2},S{2});
    vno = normal(dom{1});
    vni = -1.*normal(dom{2}); % note negative sign 
    nogradSsigmao = dot(vno,gradSsigmao);
    nigradSsigmai = dot(vni,gradSsigmai);
    nogradSsigmai2o = dot(vno,gradSsigmai2o);
    nigradSsigmao2i = dot(vni,gradSsigmao2i);

    sigmao = array_to_surfacefun(sigmavalso,dom{1},S{1});
    sigmai = array_to_surfacefun(sigmavalsi,dom{2},S{2});

    % construct Bsigma
    if abs(zk) > eps
        % compute m0
        m0o = TaylorState.debyem0(sigmao,zk,L{1},vno);
        m0i = TaylorState.debyem0(sigmai,zk,L{2},vni);
        m0ovals = surfacefun_to_array(m0o,dom{1},S{1});
        m0ivals = surfacefun_to_array(m0i,dom{2},S{2});

        % compute n . Sk[m0]
        Sm0o = complex(zeros(size(m0ovals)));
        Sm0i = complex(zeros(size(m0ivals)));
        Sm0i2o = complex(zeros(size(m0ovals)));
        Sm0o2i = complex(zeros(size(m0ivals)));
        for j=1:3
            Sm0o(:,j) = helm3d.dirichlet.eval(S{1},m0ovals(:,j), ...
                targinfoo,epslh,zk,dpars,optslho);
            Sm0i(:,j) = helm3d.dirichlet.eval(S{2},m0ivals(:,j), ...
                targinfoi,epslh,zk,dpars,optslhi);
            Sm0i2o(:,j) = helm3d.dirichlet.eval(S{2},m0ovals(:,j), ...
                targinfoo,epslh,zk,dpars);
            Sm0o2i(:,j) = helm3d.dirichlet.eval(S{2},m0ivals(:,j), ...
                targinfoi,epslh,zk,dpars);
        end
        Sm0o = array_to_surfacefun(Sm0o,dom{1},S{1});
        Sm0i = array_to_surfacefun(Sm0i,dom{2},S{2});
        Sm0i2o = array_to_surfacefun(Sm0i2o,dom{1},S{1});
        Sm0o2i = array_to_surfacefun(Sm0o2i,dom{2},S{2});

        % compute n . curl Sk[m0]
        curlSm0o = taylor.dynamic.eval_curlSk(S{1},zk,m0valso.', ...
            epstaylor,targinfoo,optso);
        curlSm0i = taylor.dynamic.eval_curlSk(S{2},zk,m0valsi.', ...
            epstaylor,targinfoi,optsi);
        curlSm0i2o = taylor.dynamic.eval_curlSk(S{2},zk,m0valsi.', ...
            epstaylor,targinfoo);
        curlSm0o2i = taylor.dynamic.eval_curlSk(S{1},zk,m0valso.', ...
            epstaylor,targinfoi);
        curlSm0o = array_to_surfacefun(curlSm0o.',dom{1},S{1});
        curlSm0i = array_to_surfacefun(curlSm0i.',dom{2},S{2});
        curlSm0i2o = array_to_surfacefun(curlSm0i2o.',dom{1},S{1});
        curlSm0o2i = array_to_surfacefun(curlSm0o2i.',dom{2},S{2});

        % combine
        m0termso = 1i.*dot(vno, zk.*(Sm0o + Sm0i2o) + curlSm0o + curlSm0i2o);
        m0termsi = 1i.*dot(vni, zk.*(Sm0i + Sm0o2i) + curlSm0i + curlSm0o2i);
        Bsigma = {sigmao./2 + nogradSsigmao + nogradSsigmai2o - m0termso, ...
            sigmai./2 + nigradSsigmai + nigradSsigmao2i - m0termsi};
    else
        Bsigma = {sigmao./2 + nogradSsigmao + nogradSsigmai2o, ...
            sigmai./2 + nigradSsigmai + nigradSsigmao2i};
    end

end

end
