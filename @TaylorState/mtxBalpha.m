function Balpha = mtxBalpha(S,dom,mH,zk,epstaylor,epslh,varargin)
%MTXBALPHA compute alpha coefficient for surface magnetic field
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * domparams: parameter describing dom and circulation [io]
%         io: [int] if 1, negate vn because inner torus
%     * mH: [surfacefunv] density for which 
%                  n . curl S0[mH]
%           is computed. Note that there is no factor of i. 
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

    if nargin < 7
        targinfo = S{1};
    else
        targinfo = varargin{1}{1};
    end
    
    if nargin < 8
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
    
    if nargin < 9
        if abs(zk) < eps
            Qlh = lap3d.dirichlet.get_quadrature_correction(S, ...
                epslh,dpars,targinfo,opts);
            optslh = [];
            optslh.format = 'rsc';
            optslh.precomp_quadrature = Qlh;
        end
    else
        optslh = varargin{3};
    end
    
    mHvals = surfacefun_to_array(mH{1},dom{1},S{1});
    mHvals = mHvals.';
    vn = normal(dom{1});
    
    if abs(zk) > eps
        % compute curl Sk[mH]
        curlSmH = taylor.dynamic.eval_curlSk(S{1},zk,mHvals,epstaylor, ...
            targinfo,opts);
        curlSmH = array_to_surfacefun(curlSmH.',dom{1},S{1});
    
        % compute Sk[mH]
        SmH = complex(zeros(size(mHvals)));
        for j=1:3
            SmH(j,:) = helm3d.dirichlet.eval(S{1},mHvals(j,:),targinfo, ...
                epslh,zk,dpars,optslh);
        end
        SmH = array_to_surfacefun(SmH.',dom{1},S{1});
    
        Balpha = {dot(vn,zk.*SmH + curlSmH)};
    else
        curlSmH = taylor.static.eval_curlS0(S{1},mHvals,epstaylor, ...
            targinfo,opts);
        curlSmH = array_to_surfacefun(curlSmH.',dom{1},S{1});
        Balpha = {dot(vn,curlSmH)};
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

    % ====

    mHvalso = surfacefun_to_array(mH{1},dom{1},S{1});
    mHvalso = mHvalso.';
    mHvalsi = surfacefun_to_array(mH{2},dom{2},S{2});
    mHvalsi = mHvalsi.';
    vno = normal(dom{1});
    vni = -normal(dom{2});
    
    if abs(zk) > eps
        % compute curl Sk[mH]
        curlSmHo = taylor.dynamic.eval_curlSk(S{1},zk,mHvalso,epstaylor, ...
            targinfoo,optso);
        curlSmHo = array_to_surfacefun(curlSmHo.',dom{1},S{1});
        curlSmHi = taylor.dynamic.eval_curlSk(S{2},zk,mHvalsi,epstaylor, ...
            targinfoi,optsi);
        curlSmHi = array_to_surfacefun(curlSmHi.',dom{2},S{2});
        curlSmHi2o = taylor.dynamic.eval_curlSk(S{2},zk,mHvalsi,epstaylor, ...
            targinfoo);
        curlSmHi2o = array_to_surfacefun(curlSmHi2o.',dom{1},S{1});
        curlSmHo2i = taylor.dynamic.eval_curlSk(S{1},zk,mHvalso,epstaylor, ...
            targinfoi);
        curlSmHo2i = array_to_surfacefun(curlSmHo2i.',dom{2},S{2});
    
        % compute Sk[mH]
        SmHo = complex(zeros(size(mHvalso)));
        SmHi = complex(zeros(size(mHvalsi)));
        SmHi2o = complex(zeros(size(mHvalso)));
        SmHo2i = complex(zeros(size(mHvalsi)));
        for j=1:3
            SmHo(j,:) = helm3d.dirichlet.eval(S{1},mHvalso(j,:),targinfoo, ...
                epslh,zk,dpars,optslho);
            SmHi(j,:) = helm3d.dirichlet.eval(S{2},mHvalsi(j,:),targinfoi, ...
                epslh,zk,dpars,optslhi);
            SmHi2o(j,:) = helm3d.dirichlet.eval(S{2},mHvalsi(j,:),targinfoo, ...
                epslh,zk,dpars);
            SmHo2i(j,:) = helm3d.dirichlet.eval(S{1},mHvalso(j,:),targinfoi, ...
                epslh,zk,dpars);
        end
        SmHo = array_to_surfacefun(SmHo.',dom{1},S{1});
        SmHi = array_to_surfacefun(SmHi.',dom{2},S{2});
        SmHi2o = array_to_surfacefun(SmHi2o.',dom{1},S{1});
        SmHo2i = array_to_surfacefun(SmHo2i.',dom{2},S{2});
    
        Balpha = {dot(vno,zk.*SmHo + curlSmHo + zk.*SmHi2o + curlSmHi2o), ...
            dot(vni,zk.*SmHi + curlSmHi + zk.*SmHo2i + curlSmHo2i)};
    else
        curlSmHo = taylor.static.eval_curlS0(S{1},mHvalso,epstaylor, ...
            targinfoo,optso);
        curlSmHo = array_to_surfacefun(curlSmHo.',dom{1},S{1});
        curlSmHi = taylor.static.eval_curlS0(S{2},mHvalsi,epstaylor, ...
            targinfoi,optsi);
        curlSmHi = array_to_surfacefun(curlSmHi.',dom{2},S{2});
        curlSmHi2o = taylor.static.eval_curlS0(S{2},mHvalsi,epstaylor, ...
            targinfoo);
        curlSmHi2o = array_to_surfacefun(curlSmHi2o.',dom{1},S{1});
        curlSmHo2i = taylor.static.eval_curlS0(S{1},mHvalso,epstaylor, ...
            targinfoi);
        curlSmHo2i = array_to_surfacefun(curlSmHo2i.',dom{2},S{2});

        Balpha = {dot(vno,curlSmHo + curlSmHi2o), ...
            dot(vni,curlSmHi + curlSmHo2i)};
    end

end

end