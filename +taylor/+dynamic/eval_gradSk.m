function gradrho = eval_gradSk(S,zk,rho,eps,varargin)
%EVAL_GRADSK compute grad S0[rho]
%
%  Input arguments:
%    * S: surfer object, see README.md in matlab for details
%    * zk: wavenumber k for layer potential Sk
%    * rho: layer potential density for which grad Sk[rho] is computed
%    * eps: precision requested
%    * targinfo: target info (optional)
%       targinfo.r = (3,nt) target locations
%       targinfo.du = u tangential derivative info
%       targinfo.dv = v tangential derivative info
%       targinfo.n = normal info
%       targinfo.patch_id (nt,) patch id of target, = -1, if target
%          is off-surface (optional)
%       targinfo.uvs_targ (2,nt) local uv ccordinates of target on
%          patch if on-surface (optional)
%    * opts: options struct
%        opts.nonsmoothonly - use smooth quadrature rule for evaluating
%           layer potential (false)
%        opts.precomp_quadrature - computed quadrature corrections struct 
%           currently only supports quadrature corrections
%           computed in rsc format 
% 

    if(nargin < 6) 
      opts = [];
    else
      opts = varargin{2};
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
        opts_qcorr.type = 'double';
        Q = init_empty_quadrature_correction(targinfo,opts_qcorr);
      end
    end

% Extract arrays
    [srcvals,srccoefs,norders,ixyzs,iptype,wts] = extract_arrays(S);
    [n12,npts] = size(srcvals);
    [n9,~] = size(srccoefs);
    [npatches,~] = size(norders);
    npatp1 = npatches+1;

    if(nargin < 5)
      targinfo = [];
      targinfo.r = S.r;
      targinfo.du = S.du;
      targinfo.dv = S.dv;
      targinfo.n = S.n;
      patch_id  = zeros(npts,1);
      uvs_targ = zeros(2,npts);
      mex_id_ = 'get_patch_id_uvs(i int[x], i int[x], i int[x], i int[x], i int[x], io int[x], io double[xx])';
[patch_id, uvs_targ] = gradcurlSk(mex_id_, npatches, norders, ixyzs, iptype, npts, patch_id, uvs_targ, 1, npatches, npatp1, npatches, 1, npts, 2, npts);
      targinfo.patch_id = patch_id;
      targinfo.uvs_targ = uvs_targ;
      opts = [];
    else
      targinfo = varargin{1};
    end

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

        [Q] = taylor.dynamic.get_quadrature_correction(S,zk,eps,targinfo,opts_quad);
      else
        opts_qcorr = [];
        opts_qcorr.type = 'double';
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

    gradrho = complex(zeros(3,ntarg));

% Call layer potential evaluator
    mex_id_ = 'lpcomp_gradhelm_addsub(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i double[x], i dcomplex[x], i int[x], i int[x], i int[x], i int[x], i int[x], i dcomplex[xx], i dcomplex[x], i int[x], i int[x], i int[x], i double[xx], i double[x], io dcomplex[xx])';
[gradrho] = gradcurlSk(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, zk, nnz, row_ptr, col_ind, iquad, nquad, wnear, rho, novers, nptso, ixyzso, srcover, wover, gradrho, 1, npatches, npatp1, npatches, 1, n9, npts, n12, npts, 1, 1, ndtarg, ntarg, 1, 1, 1, ntargp1, nnz, nnzp1, 1, nquad, 3, npts, npatches, 1, npatp1, 12, nptso, nptso, 3, ntarg);
    
end


