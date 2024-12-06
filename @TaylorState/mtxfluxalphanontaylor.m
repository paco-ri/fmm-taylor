function fluxalpha = mtxfluxalphanontaylor(S,dom,nodes,weights, ...
    mH,zk,epstaylor,epslh,varargin)
%MTXFLUXALPHANONTAYLOR compute alpha coefficient in flux condition for
%                      a non-Taylor-state magnetic field
% 
%   Required arguments:
%     * S: surfer object (see fmm3dbie/matlab README for details)
%     * dom: surfacemesh version of S (see surfacehps for details)
%     * nodes: [double(3,*)] quadrature nodes for computing flux 
%              Gauss-Legendre in r, trapezoidal in theta
%     * weights: [double(*)] corresponding quadrature weights
%     * mH: [surfacefunv] density for which 
%                  \oint curl S0[mH] . da
%              is computed
%     * zk: [double complex] wavenumber
%     * epstaylor: [double] precision requested for taylor routines
%     * epslh: [double] precision requested for lap3d/helm3d routines
% 
%   Optional arguments:
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

dpars = [1.0, 0];

nreqarg = 8;

if isa(dom,'surfacemesh')

    targinfo = [];
    targinfo.r = nodes;
    if nargin < nreqarg + 1
        targinfo = [];
        targinfo.r = nodes;
        opts = [];
        opts.format = 'rsc';
        if abs(zk) < eps
            Q = taylor.static.get_quadrature_correction(S,epstaylor, ...
                targinfo,opts);
        else
            Q = taylor.dynamic.get_quadrature_correction(S,zk,epstaylor, ...
                targinfo,opts);
        end
        opts.precomp_quadrature = Q;
    else
        opts = varargin{1};
    end
    
    if abs(zk) >= eps
        if nargin < nreqarg + 2
            optslh = [];
            optslh.format = 'rsc';
            Qlh = helm3d.dirichlet.get_quadrature_correction(S,epslh,zk, ...
                dpars,targinfo,optslh);
            optslh.precomp_quadrature = Qlh;
        else
            optslh = varargin{2};
        end
    end
    
    mHvals = surfacefun_to_array(mH,dom,S);
    mHvals = mHvals.';
    
    % Evaluate layer potential
    if abs(zk) < eps
        mHterms = taylor.static.eval_curlS0(S,mHvals,epstaylor,targinfo,opts);
    else
        curlSmH = taylor.dynamic.eval_curlSk(S,zk,mHvals,epstaylor, ...
            targinfo,opts);
        % compute Sk[mH]
        SmH = complex(zeros(size(nodes)));
        for j=1:3
            SmH(j,:) = helm3d.dirichlet.eval(S,mHvals(j,:),targinfo,epslh, ...
                zk,[1.0 0],optslh);
        end
        mHterms = curlSmH + zk.*SmH;
    end
    
    % Compute flux
    % ASSUMES y^ IS THE NORMAL TO THE CROSS-SECTION
    % fluxalpha = dot(mHterms(2,:),weights);
    fluxalpha = sum(mHterms(2,:).*weights);

else

    targtor = [];
    targtor.r = nodes{1}; 
    targpol = [];
    targpol.r = nodes{2};
    targinfo = {targtor,targpol};

    So = S{1};
    Si = S{2};

    if nargin < nreqarg + 1
        if abs(zk) < eps
            Qotor = taylor.static.get_quadrature_correction(So, ...
                epstaylor,targtor);
            Qitor = taylor.static.get_quadrature_correction(Si, ...
                epstaylor,targtor);
            Qopol = taylor.static.get_quadrature_correction(So, ...
                epstaylor,targpol);
            Qipol = taylor.static.get_quadrature_correction(Si, ...
                epstaylor,targpol);
        else
            Qotor = taylor.dynamic.get_quadrature_correction(So,zk, ...
                epstaylor,targtor);
            Qitor = taylor.dynamic.get_quadrature_correction(Si,zk, ...
                epstaylor,targtor);
            Qopol = taylor.dynamic.get_quadrature_correction(So,zk, ...
                epstaylor,targpol);
            Qipol = taylor.dynamic.get_quadrature_correction(Si,zk, ...
                epstaylor,targpol);
        end
        optsotor = [];
        optsotor.format = 'rsc';
        optsotor.precomp_quadrature = Qotor;
        optsitor = [];
        optsitor.format = 'rsc';
        optsitor.precomp_quadrature = Qitor;
        optsopol = [];
        optsopol.format = 'rsc';
        optsopol.precomp_quadrature = Qopol;
        optsipol = [];
        optsipol.format = 'rsc';
        optsipol.precomp_quadrature = Qipol;
        opts = cell(2);
        opts{1,1} = optsotor;
        opts{1,2} = optsitor;
        opts{2,1} = optsopol;
        opts{2,2} = optsipol;
    else
        opts = varargin{1};
    end
    
    if abs(zk) >= eps
        if nargin < nreqarg + 2
            Qlhotor = helm3d.dirichlet.get_quadrature_correction(So, ...
                epslh,zk,dpars,targtor);
            Qlhitor = helm3d.dirichlet.get_quadrature_correction(Si, ...
                epslh,zk,dpars,targtor);
            Qlhopol = helm3d.dirichlet.get_quadrature_correction(So, ...
                epslh,zk,dpars,targpol);
            Qlhipol = helm3d.dirichlet.get_quadrature_correction(Si, ...
                epslh,zk,dpars,targpol);
            optslhotor = [];
            optslhotor.format = 'rsc';
            optslhotor.precomp_quadrature = Qlhotor;
            optslhitor = [];
            optslhitor.format = 'rsc';
            optslhitor.precomp_quadrature = Qlhitor;
            optslhopol = [];
            optslhopol.format = 'rsc';
            optslhopol.precomp_quadrature = Qlhopol;
            optslhipol = [];
            optslhipol.format = 'rsc';
            optslhipol.precomp_quadrature = Qlhipol;
            optslh = cell(2);
            optslh{1,1} = optslhotor;
            optslh{1,2} = optslhitor;
            optslh{2,1} = optslhopol;
            optslh{2,2} = optslhipol;
        else
            optslh = varargin{2};
        end
    end

    if abs(zk) < eps
        fluxalpha = zeros(2);
        for i = 1:2
            mHvals = surfacefun_to_array(mH{i},dom{i},S{i}); % mH on surface i
            for j = 1:2
                curlSmH = taylor.static.eval_curlS0(S{i},mHvals.',...
                    epstaylor,targinfo{j},opts{j,i});
                fluxalpha(j,i) = (-1)^(j-1).*sum(curlSmH(j+1,:).*weights{j});
            end
        end
    else
        fluxalpha = zeros(2);
        for i = 1:2
            mHvals = surfacefun_to_array(mH{i},dom{i},S{i}); % mH on surface i
            for j = 1:2
                curlSmH = taylor.dynamic.eval_curlSk(S{i},zk,mHvals.',...
                    epstaylor,targinfo{j},opts{j,i});
                SmH = taylor.helper.helm_dir_vec_eval(S{i},mHvals.',...
                    targinfo{j},epslh,zk,dpars,optslh{j,i});
                integrand = curlSmH + zk.*SmH;
                fluxalpha(j,i) = (-1)^(j-1).*sum(integrand(j+1,:).*weights{j});
            end
        end
    end
end