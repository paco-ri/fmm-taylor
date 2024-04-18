function [curlj,gradrho] = gradcurlS0(norders,ixyzs,iptype,srccoefs,srcvals,targs,ipatch_id,uvs_targ,eps,rjvec,rho)
%GRADCURLS0 compute grad S0[rho] and curl S0[J]

npatches = length(norders);
npts = size(srccoefs,2);
[ndtarg,ntarg] = size(targs);
npatchesplusone = npatches+1;

mex_id_ = 'lpcomp_virtualcasing(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i int[x], i double[xx], i double[x], i double[xx], i double[x], io double[xx], io double[xx])';
[curlj, gradrho] = gradcurlS0(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, ipatch_id, uvs_targ, eps, rjvec, rho, curlj, gradrho, 1, npatches, npatchesplusone, npatches, 1, 9, npts, 12, npts, 1, 1, ndtarg, ntarg, ntarg, 2, ntarg, 1, 3, npts, npts, 3, npts, 3, npts);

end