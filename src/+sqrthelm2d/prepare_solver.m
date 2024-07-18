function [S] = prepare_solver(lmax, npts, dx, ckb, cfs)
%
%  Prepare solver either iterative or direvt
%  lmax: max excursion in l
%  n: number of points in each direction
%  dx: spacing in each direction
%  ckb: wave number
%  cfs: quadrature correction
%  opts: options struct (optional)
%    opts.ifflam = whether to use flam or not (false)

    S = [];
    S.lmax = lmax;
    S.npts = npts;
    S.dx = dx;
    S.ckb = ckb;
    S.cfs = cfs;
    n = npts*npts;
    S.n = n;
    
    
    S.spmat = sqrthelm2d.get_sparse_correction(npts, dx, ckb, cfs);    
    S.gmat = sqrthelm2d.build_gmat(lmax, dx, ckb, cfs);    


end