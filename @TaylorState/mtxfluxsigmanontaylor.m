function fluxsigma = mtxfluxsigmanontaylor(domain,nodes, ...
    weights,sigma,zk,epstaylor,epslh,varargin)
%MTXFLUXSIGMANONTAYLOR compute sigma-dep. terms of flux condition for
%                      a non-Taylor-state magnetic field
% 
%   Required arguments:
%     * domain: Domain object
%     * nodes: [double(3,*)] quadrature nodes for computing flux 
%              Gauss-Legendre in r, trapezoidal in theta
%     * weights: [double(*)] corresponding quadrature weights
%     * sigma: [surfacefun] density for which 
%                  \oint (grad S_0[sigma]) . da
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

S = domain.surf;
dom = domain.dom;
vn = domain.vn;
L = domain.L;

dpars = [1.0, 0];

nreqarg = 7;

if isscalar(dom) && isa(dom{1},'surfacemesh')
    dom = dom{1};
    S = S{1};
    vn = vn{1};
    L = L{1};
end

if isa(dom,'surfacemesh')

    targinfo = [];
    targinfo.r = nodes; 

    if nargin < nreqarg + 1
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
    
    sigmavals = surfacefun_to_array(sigma,dom,S);
    sigmavals = sigmavals.';
    
    % Evaluate layer potential
    if abs(zk) < eps
        sigmaterms = taylor.static.eval_gradS0(S,sigmavals,epstaylor, ...
            targinfo,opts);
    else
        gradSsigma = taylor.dynamic.eval_gradSk(S,zk,sigmavals,epstaylor, ...
            targinfo,opts);
        m0 = TaylorState.debyem0(sigma,zk,L,vn);
        m0vals = surfacefun_to_array(m0,dom,S);
        m0vals = m0vals.';
        curlSkm0 = taylor.dynamic.eval_curlSk(S,zk,m0vals,epstaylor, ...
            targinfo,opts);
        Skm0 = complex(zeros(size(nodes)));
        for j = 1:3
            Skm0(j,:) = helm3d.dirichlet.eval(S,m0vals(j,:),targinfo,epslh, ...
                zk,dpars,optslh);
        end
        sigmaterms = gradSsigma - 1i.*(zk.*Skm0 + curlSkm0);
    end
    
    % Compute flux
    % ASSUMES y^ IS THE NORMAL TO THE CROSS-SECTION
    % fluxsigma = dot(sigmaterms(2,:),weights);
    fluxsigma = sum(sigmaterms(2,:).*weights);

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

    numfuns = size(sigma,2);
    sig = cell(2,numfuns);
    for i = 1:2
        for j = 1:numfuns
            sig{i,j} = surfacefun_to_array(sigma{i,j},dom{i},S{i});
        end
    end

    if abs(zk) < eps
        gradS0sig = cell(2,2,2);
        for i = 1:2
            for j = 1:numfuns
                for k = 1:2
                    % i source (outer/inner)
                    % j column
                    % k tor/pol
                    gradS0sig{i,j,k} = taylor.static.eval_gradS0(S{i}, ...
                        sig{i,j},epstaylor,targinfo{k},opts{k,i});
                end
            end
        end
        fluxsigma = zeros(2,numfuns);
        for i = 1:2
            for j = 1:numfuns
                for k = 1:2
                    % k + 1 is 2nd or 3rd component, depending on 
                    % toroidal vs poloidal flux
                    fluxsigma(k,j) = fluxsigma(k,j) ...
                        + (-1)^(k-1).*sum(gradS0sig{i,j,k}(k+1,:).*weights{k});
                end
            end
        end
    else
        gradSksig = cell(2,2,2);
        curlSkm0 = cell(2,2,2);
        Skm0 = cell(2,2,2);
        for i = 1:2
            for j = 1:numfuns
                m0 = TaylorState.debyem0(sigma{i,j},zk,L{i},vn{i});
                m0vals = surfacefun_to_array(m0,dom{i},S{i});
                m0vals = m0vals.';
                for k = 1:2
                    % i source (outer/inner)
                    % j column
                    % k tor/pol
                    gradSksig{i,j,k} = taylor.dynamic.eval_gradSk(S{i}, ...
                        zk,sig{i,j},epstaylor,targinfo{k},opts{k,i});
                    curlSkm0{i,j,k} = taylor.dynamic.eval_curlSk(S{i}, ...
                        zk,m0vals,epstaylor,targinfo{k},opts{k,i});
                    Skm0{i,j,k} = taylor.helper.helm_dir_vec_eval(S{i}, ...
                        m0vals,targinfo{k},epslh,zk,dpars,optslh{k,i});
                end
            end
        end
        fluxsigma = zeros(2,numfuns);
        for i = 1:2
            for j = 1:numfuns
                for k = 1:2
                    % k+1 is 2nd or 3rd component, depending on 
                    % toroidal vs poloidal flux 
                    integrand = gradSksig{i,j,k} - 1i.*(zk.*Skm0{i,j,k} ...
                        + curlSkm0{i,j,k});
                    fluxsigma(k,j) = fluxsigma(k,j) ...
                        + (-1)^(k-1).*sum(integrand(k+1,:).*weights{k});
                end
            end
        end
    end

end

end