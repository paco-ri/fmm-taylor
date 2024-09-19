function [p, gradrho, curlj] = mtxBsigma_eval(S,sigma,rho,rjvec,targinfo,eps,zk,rep_pars,varargin)
%
%  taylor.helper.mtxBsigma_eval
%    Helper function for layer potential evaluations in mtxBsigma 
%
%  Syntax
%   [pot,gradrho,curlj] = helm3d.dirichlet.eval(S,sigma,rho,rjvec,targinfo,eps,zk,rep_pars)
%   [pot,gradrho,curlj] = helm3d.dirichlet.eval(S,sigma,rho,rjvec,targinfo,eps,zk,rep_pars,opts)
%
%  Integral representation
%     pot = \alpha S_{k} [\sigma] + \beta D_{k} [\sigma]
%     gradrho = grad S_{k} [\rho]
%     curlj = curl S_{k} [j]
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
%    * eps: double(2) precision requested for Sk and grad/curl Sk, resp.
%    * zk: wave number
%    * rep_pars: kernel parameters
%        zpars(1) - single layer strength
%        zpars(2) - double layer strength
%    * opts: options struct for curl and grad Sk
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%        opts.precomp_quadrature: precomputed quadrature corrections struct 
%           currently only supports quadrature corrections
%           computed in rsc format
%    * optslh: options struct for Sk
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%        opts.precomp_quadrature: precomputed quadrature corrections struct 
%           currently only supports quadrature corrections
%           computed in rsc format 
%    

    if(nargin < 9) 
      opts = [];
    else
      opts = varargin{1};
    end

    if(nargin < 10)
      optshelm = [];
    else
      optshelm = varargin{2};
    end

    nonsmoothonly = false;
    if(isfield(opts,'nonsmoothonly'))
      nonsmoothonly = opts.nonsmoothonly;
    end

    nonsmoothonlyhelm = false;
    if(isfield(optshelm,'nonsmoothonly'))
      nonsmoothonlyhelm = optshelm.nonsmoothonly;
    end

    isprecompq = false;
    if isfield(opts, 'precomp_quadrature')
      isprecompq = true;
      Q = opts.precomp_quadrature;
    end

    isprecompqhelm = false;
    if isfield(optshelm, 'precomp_quadrature')
      isprecompqhelm = true;
      Qhelm = optshelm.precomp_quadrature; 
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

    if(isprecompqhelm)
      if ~(strcmpi(Qhelm.format,'rsc'))
        fprintf('Invalid precomputed quadrature format\n');
        fprintf('Ignoring quadrature corrections\n');
        opts_qcorr = [];
        opts_qcorr.type = 'complex';
        Qhelm = init_empty_quadrature_correction(targinfo,opts_qcorr);
      end
    end


% Extract arrays
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;

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

        [Q] = taylor.dynamic.get_quadrature_correction(S,zk,eps(1),targinfo,opts_quad);
      else
        opts_qcorr = [];
        opts_qcorr.type = 'complex';
        Q = init_empty_quadrature_correction(targinfo,opts_qcorr);
      end
    end

    if ~isprecompqhelm
      if ~nonsmoothonlyhelm
        opts_quad = [];
        opts_quad.format = 'rsc';
        [Qhelm] = helm3d.dirichlet.get_quadrature_correction(S,eps(2),zk,rep_pars,targinfo,opts_quad);
      else
        opts_qcorr = [];
        opts_qcorr.type = 'complex';
        Qhelm = init_empty_quadrature_correction(targinfo,opts_qcorr);
      end
    end

    nquad = Q.iquad(end)-1;
    nquadhelm = Qhelm.iquad(end)-1;
    nnz = length(Q.col_ind);
    nnzhelm = length(Qhelm.col_ind);
    nnzp1 = nnz+1;
    nnzhelmp1 = nnzhelm+1; 

    [novers] = get_oversampling_parameters(S,Q,eps(1));
    [novershelm] = get_oversampling_parameters(S,Qhelm,eps(2));
    Sover = oversample(S,novers);
    Soverhelm = oversample(S,novershelm);


% Extract oversampled arrays

    [srcover,~,~,ixyzso,~,wover] = extract_arrays(Sover);
    [srcoverhelm,~,~,ixyzsohelm,~,woverhelm] = extract_arrays(Soverhelm);
    nptso = Sover.npts;
    nptsohelm = Soverhelm.npts; 

% Extract quadrature arrays
    row_ptr = Q.row_ptr;
    col_ind = Q.col_ind;
    iquad = Q.iquad;
    wnear = Q.wnear;

    row_ptr_helm = Qhelm.row_ptr;
    col_ind_helm = Qhelm.col_ind;
    iquadhelm = Qhelm.iquad;
    wnearhelm = Qhelm.wnear;

    ndim = size(sigma,1);
    p = complex(zeros(ndim,ntarg));
    gradrho = complex(zeros(3,ntarg));
    curlj = complex(zeros(3,ntarg));

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
    eps_helm = eps(2);
    eps_gc = eps(1);
% Call the layer potential evaluator
    mex_id_ = 'lpcomp_gradcurlhelm_addsub(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i double[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i dcomplex[xx], i dcomplex[xx], i dcomplex[x], i int[x], i int[x], i int[x], i double[xx], i double[x], io dcomplex[xx], io dcomplex[xx])';
[curlj, gradrho] = helper(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps_gc, zk, nnz, row_ptr, col_ind, iquad, nquad, wnear, rjvec, rho, novers, nptso, ixyzso, srcover, wover, curlj, gradrho, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, 1, 1, 1, ntargp1, nnz, nnzp1, 1, nquad, 3, 3, npts, npts, npatches, 1, npatp1, 12, nptso, nptso, 3, ntarg, 3, ntarg);
    mex_id_ = 'helm_comb_dir_eval_addsub_vec(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i double[x], i int[x], i double[x], i int[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i dcomplex[x], i int[x], i int[x], i int[x], i double[xx], i double[x], i int[x], i double[x], i int[x], i int[x], i dcomplex[xx], i int[x], i int[x], io dcomplex[xx])';
[p] = helper(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps_helm, ndd, dpars, ndz, zpars, ndi, ipars, nnzhelm, row_ptr_helm, col_ind_helm, iquadhelm, nquadhelm, nker, wnearhelm, novershelm, nptsohelm, ixyzsohelm, srcoverhelm, woverhelm, lwork, work, idensflag, ndim, sigma, ipotflag, ndim_p, p, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, 1, 1, ndd, 1, ndz, 1, ndi, 1, ntargp1, nnzhelm, nnzhelmp1, 1, 1, nquadhelm, npatches, 1, npatp1, 12, nptsohelm, nptsohelm, 1, lwork, 1, 1, ndim, npts, 1, 1, ndim, ntarg);
end

