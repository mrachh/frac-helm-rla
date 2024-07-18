function Sout = compute_factorization(S, V, opts)
    lmax = S.lmax;
    npts = S.npts;
    xs = linspace(-lmax, lmax, npts);
    ys = xs;
    [X,Y] = ndgrid(xs,ys);
    dx = S.dx;
    ckb = S.ckb;

    if nargin < 3
        opts = [];
    end
    p = 64;
    if isfield(opts, 'p')
        p = opts.p;
    end

    nr = 10;
    if isfield(opts, 'nr')
        nr = opts.nr;
    end
    
    tol = 1e-7;
    if isfield(opts, 'tol')
        tol = opts.tol;
    end

    occ = 2000;
    if isfield(opts, 'occ')
        occ = 2000;
    end

    verb = true;
    if isfield(opts, 'verb')
        verb = opts.verb;
    end

    [proxy, ~, pw] = proxy_circ_pts(p, nr);

    xpts = [X(:).';Y(:).'];

    pfun = @(xfun,slf,nbr,l,ctr) pxyfun_frac(xfun, slf, nbr, proxy, pw,...
        l, ctr, V(:), ckb, dx);

    rank_or_tol = tol;
    
    opts_use = [];
    opts_use.verb = verb;

    matind = @(i,j) frac_mat_fun(i, j, X, Y, V, ckb, dx, S.spmat_with_v);
    F = rskelf(matind, xpts, occ, rank_or_tol, pfun, opts_use); 
    Sout = S;
    Sout.F = F;

end



% proxy function
function [Kpxy,nbr] = pxyfun_frac(x,slf,nbr,proxy,pw,l,ctr,vin,ckb,dx)
% PXYFUN(X,SLF,NBR,L,CTR) computes interactions between the points
% X(:,SLF) and the set of proxy points by scaling the proxy sphere to 
% appropriately contain a box at level L centered at CTR and then
% calling helm_dirichlet_kernel

    pxy = bsxfun(@plus,proxy.*l,ctr);
    
    [xt,xs] = ndgrid(pxy(1,:),x(1,slf));
    [yt,ys] = ndgrid(pxy(2,:),x(2,slf));
    [wt,~ ] = ndgrid(pw,x(1,slf));
    dx1 = xt-xs;
    dx2 = yt-ys;
    Kpxy1 = sqrthelm2d.green(dx1,dx2,ckb)*dx^2;
    Kpxy2 = sqrthelm2d.green(dx1,dx2,ckb)*diag(ckb*vin(slf));
    Kpxy2 = Kpxy2.*wt;
    Kpxy = [Kpxy1;Kpxy2];
    dx = x(1,nbr) - ctr(1);
    dy = x(2,nbr) - ctr(2);
    dist = sqrt(dx.^2 + dy.^2);
    nbr = nbr(dist/norm(l) < 1.5);

end


%%
function [amat] = frac_mat_fun(i,j,x,y,v,ckb,dx,scorr)
    [XS,XT] = ndgrid(x(i),x(j));
    [YS,YT] = ndgrid(y(i),y(j));
    dx1 = XS-XT;
    dx2 = YS-YT;
    amat = sqrthelm2d.green(dx1,dx2,ckb);
    amat(isnan(amat)) = 0;
    amat = diag(-ckb*dx^2*v(i).')*amat;
    amat = amat + scorr(i,j);
end

%%


function [proxy, pnorm, pw] = proxy_circ_pts(p, nr)
    if nargin < 1
        p = 64;
    end

    if nargin < 2
        nr = 5;
    end
    
    proxy = [];
    pnorm = [];
    pw    = [];
    
    for ii=1:nr
        theta = (0:(p-1))*2*pi/p;
        rad = (1.5+(ii-1)*0.5/(nr-1));
        proxyt = rad*[cos(theta);sin(theta)];
        pnormt = [cos(theta); sin(theta)];
        pwt = ones(p,1)*(2.0*pi*rad)/p;
        
        proxy = [proxy,proxyt];
        pnormt = [pnorm,pnormt];
        pw   = [pw,pwt];
    end
end
