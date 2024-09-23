function p = helm_dir_vec_eval(S,sigma,targinfo,eps,zk,rep_pars,varargin)
%
%  taylor.helper.helm_dir_vec_eval
%    Evaluates the helmholtz dirichlet layer potential at a collection 
%    of targets
%
%  Syntax
%   pot = helm3d.dirichlet.eval(S,sigma,targinfo,eps,zk,rep_pars)
%   pot = helm3d.dirichlet.eval(S,sigma,targinfo,eps,zk,rep_pars,opts)
%
%  Integral representation
%     pot = \alpha S_{k} [\sigma] + \beta D_{k} [\sigma]
%
%  S_{k}, D_{k}: helmholtz single and double layer potential
%  
%  \alpha, beta = rep_pars(1:2)
%  k = zk
%
%  Note: for targets on surface, only principal value part of the
%    layer potential is returned
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * sigma: layer potential density
%    * targinfo: target info (optional)
%       targinfo.r = (3,nt) target locations
%       targinfo.patch_id (nt,) patch id of target, = -1, if target
%          is off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%          patch if on-surface (optional)
%    * eps: precision requested
%    * zk: wave number
%    * rep_pars: kernel parameters
%        zpars(1) - single layer strength
%        zpars(2) - double layer strength
%    * opts: options struct
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%        opts.precomp_quadrature: precomputed quadrature corrections struct 
%           currently only supports quadrature corrections
%           computed in rsc format 
%    

    if(nargin < 7) 
      opts = [];
    else
      opts = varargin{1};
    end

    nonsmoothonly = false;
    if(isfield(opts,'nonsmoothonly'))
      nonsmoothonly = opts.nonsmoothonly;
    end

    isprecompq = false;
    if isfield(opts, 'precomp_quadrature')
      isprecompq = true;
      Q = opts.precomp_quadrature;
    end
    
    if(isprecompq)
      if ~(strcmpi(Q.format,'rsc'))
        fprintf('Invalid precomputed quadrature format\n');
        fprintf('Ignoring quadrature corrections\n');
        opts_qcorr = [];
        opts_qcorr.type = 'complex';
        Q = init_empty_quadrature_correction(targinfo,opts_qcorr);
      end
    end


% Extract arrays
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npat,~] = size(norders);
    npatp1 = npat+1;

    ff = 'rsc';

    [targs] = extract_targ_array(targinfo);
    [ndtarg,ntarg] = size(targs);
    ntargp1 = ntarg+1;

% Compute quadrature corrections   
    if ~isprecompq
      if ~nonsmoothonly
        opts_quad = [];
        opts_quad.format = 'rsc';
%
%  For now Q is going to be a struct with 'quad_format', 
%  'nkernels', 'pde', 'bc', 'kernel', 'ker_order',
%  and either, 'wnear', 'row_ind', 'col_ptr', or
%  with 'spmat' as a sparse matrix or a cell array of wnear/spmat
%  if nkernel is >1
%

        [Q] = helm3d.dirichlet.get_quadrature_correction(S,eps,zk,rep_pars,targinfo,opts_quad);
      else
        opts_qcorr = [];
        opts_qcorr.type = 'complex';
        Q = init_empty_quadrature_correction(targinfo,opts_qcorr);
      end
    end
    nquad = Q.iquad(end)-1;
    nnz = length(Q.col_ind);
    nnzp1 = nnz+1; 

    [novers] = get_oversampling_parameters(S,Q,eps);
    Sover = oversample(S,novers);


% Extract oversampled arrays

    [srcover,~,~,ixyzso,~,wover] = extract_arrays(Sover);
    nptso = Sover.npts; 

% Extract quadrature arrays
    row_ptr = Q.row_ptr;
    col_ind = Q.col_ind;
    iquad = Q.iquad;
    wnear = Q.wnear;

    ndim = size(sigma,1);
    p = complex(zeros(ndim,ntarg));

    zpars = complex(zeros(3,1));
    zpars(1) = zk;
    zpars(2) = rep_pars(1);
    zpars(3) = rep_pars(2);

    ndd = 0;
    dpars = [];
    ndz = 3;
    ndi = 0;
    ipars = [];
    nker = 1;
    lwork = 0;
    work = [];
    idensflag = 0;
    ipotflag = 0;
    ndim_p = 1;
% Call the layer potential evaluator
    mex_id_ = 'helm_comb_dir_eval_addsub_vec(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i double[x], i int[x], i double[x], i int[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i dcomplex[x], i int[x], i int[x], i int[x], i double[xx], i double[x], i int[x], i double[x], i int[x], i int[x], i dcomplex[xx], i int[x], i int[x], io dcomplex[xx])';
[p] = helper(mex_id_, npat, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, nker, wnear, novers, nptso, ixyzso, srcover, wover, lwork, work, idensflag, ndim, sigma, ipotflag, ndim_p, p, 1, npat, npatp1, npat, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, 1, 1, ndd, 1, ndz, 1, ndi, 1, ntargp1, nnz, nnzp1, 1, 1, nquad, npat, 1, npatp1, 12, nptso, nptso, 1, lwork, 1, 1, ndim, npts, 1, 1, ndim, ntarg);
end


