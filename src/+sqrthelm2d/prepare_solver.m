function [S] = prepare_solver(lmax, npts, ckb, cfs)
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
    
    S.ckb = ckb;
    S.cfs = cfs;
    n = npts*npts;
    S.n = n;

    xs = linspace(-lmax, lmax, npts);
    ys = xs;
    [X,Y] = ndgrid(xs,ys);
    xpts = [X(:).';Y(:).'];
    S.xpts = xpts;
    S.dx = xs(2)-xs(1);
    
    S.spmat = sqrthelm2d.get_sparse_correction(npts, S.dx, ckb, cfs);    
    S.gmat = sqrthelm2d.build_gmat(lmax, S.dx, ckb, cfs);    


end