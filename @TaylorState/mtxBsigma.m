function Bsigma = mtxBsigma(domain,sigmavals,zk,epstaylor,epslh,varargin)
%MTXBSIGMA compute sigma-dep. terms of surface magnetic field
% 
%   Required arguments:
%     * domain: Domain object
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
S = domain.surf;
dom = domain.dom;
vn = domain.vn;
L = domain.L;

nreqarg = 5;
% torus case
if length(dom) == 1
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

    sigma = array_to_surfacefun(sigmavals,dom{1},S{1});

    if ~isa(vn,'surfacefunv')
        vn = vn{1};
    end

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

        [Sm0, gradSsigma, curlSm0] = taylor.helper.mtxBsigma_eval(S{1}, ...
            m0vals.',sigmavals,m0vals.',targinfo,[epstaylor,epslh],zk, ...
            dpars,opts,optslh);
        Sm0 = array_to_surfacefun(Sm0.',dom{1},S{1});
        gradSsigma = array_to_surfacefun(gradSsigma.',dom{1},S{1}); 
        ngradSsigma = dot(vn,gradSsigma);
        curlSm0 = array_to_surfacefun(curlSm0.',dom{1},S{1});

        % combine
        m0terms = 1i*zk.*dot(vn,Sm0) + 1i.*dot(vn,curlSm0);
        Bsigma = {sigma./2 + ngradSsigma - m0terms};
    end

% toroidal shell case    
else
    
    % targinfoo = outer surface, targinfoi = inner surface
    if nargin < nreqarg + 1
        targinfo = S;
    else
        targinfo = varargin{1};
    end

    % near quadrature corrections for +taylor routines
    if nargin < nreqarg + 2
        opts = cell(2);
        if abs(zk) < eps
            for i = 1:2
                for j = 1:2
                    Q = taylor.static.get_quadrature_correction( ...
                        S{i},epstaylor,targinfo{j});
                    opts{i,j} = [];
                    opts{i,j}.format = 'rsc';
                    opts{i,j}.precomp_quadrature = Q;
                end
            end
        else
            for i = 1:2
                for j = 1:2
                    Q = taylor.dynamic.get_quadrature_correction( ...
                        S{i},zk,epstaylor,targinfo{j});
                    opts{i,j} = [];
                    opts{i,j}.format = 'rsc';
                    opts{i,j}.precomp_quadrature = Q;
                end
            end
        end
    else
        opts = varargin{2};
    end

    % near quadrature corrections for Helmholtz layer potential
    % only needed if zk ~= 0 
    if abs(zk) >= eps
        if nargin < nreqarg + 3
            optslh = cell(2);
            for i = 1:2
                for j = 1:2
                    Qlh = helm3d.dirichlet.get_quadrature_correction( ...
                        S{i},epslh,zk,dpars,targinfo{j});
                    optslh{i,j} = [];
                    optslh{i,j}.format = 'rsc';
                    optslh{i,j}.precomp_quadrature = Qlh;
                end
            end
        else
            optslh = varargin{3};
        end
    end

    % split sigmavals into halves. top half is on outer surface, etc.
    npts = length(sigmavals)/2;
    sigvalcell = cell(1,2);
    sigvalcell{1} = sigmavals(1:npts);
    sigvalcell{2} = sigmavals(npts+1:end);

    % evaluate layer potential
    ngradSsigma = cell(2);
    if abs(zk) < eps
        for i = 1:2
            for j = 1:2
                gradSsigma = taylor.static.eval_gradS0(S{i}, ...
                    sigvalcell{i},epstaylor,targinfo{j},opts{i,j});
                gradSsigma = array_to_surfacefun(gradSsigma.', ...
                    dom{j},S{j});
                ngradSsigma{i,j} = dot(vn{j},gradSsigma);
            end
        end
    else
        for i = 1:2
            for j = 1:2
                gradSsigma = taylor.dynamic.eval_gradSk(S{i},zk, ...
                    sigvalcell{i},epstaylor,targinfo{j},opts{i,j});
                gradSsigma = array_to_surfacefun(gradSsigma.', ...
                    dom{j},S{j});
                ngradSsigma{i,j} = dot(vn{j},gradSsigma);
            end
        end
    end

    sigma = cell(1,2);
    for i=1:2
        sigma{i} = array_to_surfacefun(sigvalcell{i},dom{i},S{i});
    end

    % construct Bsigma
    Bsigma = cell(1,2);
    if abs(zk) < eps
        for i = 1:2
            Bsigma{i} = sigma{i}./2;
            for j = 1:2
                Bsigma{i} = Bsigma{i} + ngradSsigma{j,i};
            end
        end
    else
        % compute m0
        Sm0 = cell(2);
        curlSm0 = cell(2);
        for i = 1:2
            m0 = TaylorState.debyem0(sigma{i},zk,L{i},vn{i});
            m0 = surfacefun_to_array(m0,dom{i},S{i});
            for j = 1:2
                Sm0{i,j} = taylor.helper.helm_dir_vec_eval(S{i},m0.', ...
                    targinfo{j},epslh,zk,dpars,optslh{i,j});
                Sm0{i,j} = array_to_surfacefun(Sm0{i,j}.',dom{j},S{j});
                curlSm0{i,j} = taylor.dynamic.eval_curlSk(S{i},zk,m0.', ...
                    epstaylor,targinfo{j},opts{i,j});
                curlSm0{i,j} = array_to_surfacefun(curlSm0{i,j}.',dom{j},S{j});
            end
        end

        for i = 1:2
            Bsigma{i} = sigma{i}./2;
            m0terms = 0;
            for j = 1:2
                Bsigma{i} = Bsigma{i} + ngradSsigma{j,i};
                m0terms = m0terms + zk.*Sm0{j,i} + curlSm0{j,i};
            end
            Bsigma{i} = Bsigma{i} - 1i.*dot(vn{i},m0terms);
        end
    end

end

end
