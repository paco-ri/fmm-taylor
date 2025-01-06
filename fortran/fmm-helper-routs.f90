subroutine lap_comb_dir_eval_addsub_vec(npatches, norders, ixyzs, &
     iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, &
     dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
     nquad, nker, wnear, novers, nptso, ixyzso, srcover, whtsover, &
     lwork, work, idensflag, ndim, sigma, ipotflag, ndim_p, pot)

  ! Vectorized and complex sigma version of lpcomp_helm_comb_dir_addsub
  !
  ! Changed input arguments:
  !   ndim: integer
  !     number of densities for which to evaluate layer potential
  !   sigma: complex *16 (ndim,npts)
  !
  ! Changed output argument:
  !   pot: complex *16 (ndim,npts)

  implicit none
  integer npatches,npts,ndd,ndz,ndi,nnz,ndtarg
  integer norders(npatches),ixyzs(npatches+1),iptype(npatches),ipars(ndi)
  real *8 srccoefs(9,npts),srcvals(12,npts),eps,dpars(ndd)
  complex *16 zpars(ndz)
  integer ntarg,nquad,nker,nptso,lwork,ndim
  integer row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1),novers(npatches+1)
  complex *16 wnear(nquad),sigma(ndim,npts)
  integer ixyzso(npatches+1)
  real *8 targs(ndtarg,ntarg),srcover(12,nptso),whtsover(nptso),work(lwork)
  integer idensflag, ipotflag, ndim_p
  complex *16 pot(ndim,ntarg)

  real *8 potr(ndim,ntarg),poti(ndim,ntarg)

  integer i
  
  do i = 1,ndim
     call lap_comb_dir_eval_addsub(npatches, norders, ixyzs, &
     iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, &
     dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
     nquad, nker, wnear, novers, nptso, ixyzso, srcover, whtsover, &
     lwork, work, idensflag, 1, real(sigma(i,:)), ipotflag, ndim_p, potr(i,:))
     call lap_comb_dir_eval_addsub(npatches, norders, ixyzs, &
     iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, &
     dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
     nquad, nker, wnear, novers, nptso, ixyzso, srcover, whtsover, &
     lwork, work, idensflag, 1, aimag(sigma(i,:)), ipotflag, ndim_p, poti(i,:))
  enddo

  pot = dcmplx(potr,poti)

  return
end subroutine lap_comb_dir_eval_addsub_vec


subroutine helm_comb_dir_eval_addsub_vec(npatches, norders, ixyzs, &
     iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, &
     dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
     nquad, nker, wnear, novers, nptso, ixyzso, srcover, whtsover, &
     lwork, work, idensflag, ndim, sigma, ipotflag, ndim_p, pot)

  ! Vectorized version of lpcomp_helm_comb_dir_addsub
  !
  ! Changed input arguments:
  !   ndim: integer
  !     number of densities for which to evaluate layer potential
  !   sigma: complex *16 (ndim,npts)
  !
  ! Changed output argument:
  !   pot: complex *16 (ndim,npts)

  implicit none
  integer npatches,npts,ndd,ndz,ndi,nnz,ndtarg
  integer norders(npatches),ixyzs(npatches+1),iptype(npatches),ipars(ndi)
  real *8 srccoefs(9,npts),srcvals(12,npts),eps,dpars(ndd)
  complex *16 zpars(ndz)
  integer ntarg,nquad,nker,nptso,lwork,ndim
  integer row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1),novers(npatches+1)
  complex *16 wnear(nquad),sigma(ndim,npts)
  integer ixyzso(npatches+1)
  real *8 targs(ndtarg,ntarg),srcover(12,nptso),whtsover(nptso),work(lwork)
  integer idensflag, ipotflag, ndim_p
  complex *16 pot(ndim,ntarg)

  integer i
  
  do i = 1,ndim
     call helm_comb_dir_eval_addsub(npatches, norders, ixyzs, &
     iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, eps, ndd, &
     dpars, ndz, zpars, ndi, ipars, nnz, row_ptr, col_ind, iquad, &
     nquad, nker, wnear, novers, nptso, ixyzso, srcover, whtsover, &
     lwork, work, idensflag, 1, sigma(i,:), ipotflag, ndim_p, pot(i,:))
  enddo

  return
end subroutine helm_comb_dir_eval_addsub_vec

subroutine mtxbsigma_eval_addsub(npatches,norders,ixyzs, &
     iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs, &
     eps,zpars,nnz_gc,row_ptr_gc, &
     col_ind_gc,iquad_gc,nquad_gc,wnear_gc, &
     novers_gc,nptso_gc,ixyzso_gc,srcover_gc,whtsover_gc, &
     nnz_helm,row_ptr_helm,col_ind_helm,iquad_helm, &
     nquad_helm,wnear_helm,novers_helm,nptso_helm, &
     ixyzso_helm,srcover_helm,whtsover_helm, &
     ndim,sigma,rjvec,rho,pot,curlgrad)
  
  implicit none
  integer npatches,npts,nnz_gc,ndtarg
  integer norders(npatches),ixyzs(npatches+1)
  integer iptype(npatches)
  integer ipars
  real *8 srccoefs(9,npts),srcvals(12,npts),eps(2)
  real *8 dpars
  complex *16 zpars(3)
  integer ntarg,nquad_gc,nptso_gc
  integer ndim
  
  integer row_ptr_gc(ntarg+1),col_ind_gc(nnz_gc)
  integer iquad_gc(nnz_gc+1),novers_gc(npatches)
  complex *16 wnear_gc(nquad_gc,3)
  integer ixyzso_gc(npatches+1)
  real *8 srcover_gc(12,nptso_gc)
  real *8 whtsover_gc(nptso_gc)
  
  complex *16 sigma(ndim,npts)
  
  real *8 targs(ndtarg,ntarg)
  real *8 work
  complex *16 pot(ndim,ntarg)

  integer ixyzso_helm(npatches+1)
  integer nnz_helm,row_ptr_helm(ntarg+1)
  integer col_ind_helm(nnz_helm),nquad_helm
  integer iquad_helm(nnz_helm+1)
  complex *16 rjvec(3,npts),rho(npts)
  complex *16 wnear_helm(nquad_helm)

  integer novers_helm(npatches),nptso_helm
  real *8 srcover_helm(12,nptso_helm),whtsover_helm(nptso_helm)
  complex *16 curlj(3,ntarg),gradrho(3,ntarg),curlgrad(6,ntarg)

  integer i
  
  call lpcomp_gradcurlhelm_addsub(npatches,norders,ixyzs,&
     iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs, &
     eps(1),zpars(1),nnz_gc,row_ptr_gc,col_ind_gc,iquad_gc,&
     nquad_gc,wnear_gc,rjvec,rho,novers_gc,nptso_gc,ixyzso_gc,&
     srcover_gc,whtsover_gc,curlj,gradrho)
  !$OMP PARALLEL DO DEFAULT(SHARED)
  do i = 1,npts
     curlgrad(1:3,i) = curlj(1:3,i)
     curlgrad(4:6,i) = gradrho(1:3,i)
  enddo
  !$OMP END PARALLEL DO
  call helm_comb_dir_eval_addsub_vec(npatches,norders,ixyzs, &
     iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,eps(2),0, &
     dpars,3,zpars,0,ipars,nnz_helm,row_ptr_helm,col_ind_helm, &
     iquad_helm,nquad_helm,1,wnear_helm,novers_helm,nptso_helm, &
     ixyzso_helm,srcover_helm,whtsover_helm,0,work,0,ndim,sigma, &
     0,1,pot)

  return
end subroutine mtxbsigma_eval_addsub










