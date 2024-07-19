function [xout,params] = solve_it(ckb,npts,lmax)
xs = linspace(-lmax,lmax,npts);
ys = xs;
[X,Y] = ndgrid(xs,ys);
dx = xs(2)-xs(1);


cfs = get_corrs(ckb,lmax,npts);
cfs_hd = get_corrs_hd(ckb,lmax,npts);
cfs = cfs_hd;

nbig = 2^ceil(log(ceil((4*lmax)/dx))/log(2))
xlrg = (((-nbig)/2+1):nbig/2)*dx;
cent = find(xlrg == 0);
xlrg = circshift(xlrg,cent+1);

[XL,YL]=ndgrid(xlrg,xlrg);
gmat = gfunc(XL,YL,ckb);
gmat(isnan(gmat)) = 0;

gmat(1,1)
gmat(1,1) = cfs(1)/(dx^2);
[XL(2,1),YL(2,1)]
gmat(2,1) = gmat(2,1)+cfs(2)/(2*dx^2);
[XL(end,1),YL(end,1)]
gmat(end,1)= gmat(end,1)+cfs(2)/(2*dx^2);
[XL(1,2),YL(1,2)]
gmat(3,1) = gmat(3,1)+cfs(3)/(2*dx^2);
[XL(1,end),YL(1,end)]
gmat(end-1,1) = gmat(end-1,1)+cfs(3)/(2*dx^2);

gmat(1,2) = gmat(1,2)+cfs(4)/(2*dx^2);
[XL(end,1),YL(end,1)]
gmat(1,end)= gmat(1,end)+cfs(4)/(2*dx^2);
[XL(1,2),YL(1,2)]
gmat(1,3) = gmat(1,3)+cfs(5)/(2*dx^2);
[XL(1,end),YL(1,end)]
gmat(1,end-1) = gmat(1,end-1)+cfs(5)/(2*dx^2);

gmat(2,2)     = gmat(2,2)    +cfs(6)/(4*dx^2);
gmat(2,end)   = gmat(2,end)  +cfs(6)/(4*dx^2);
gmat(end,2)   = gmat(end,2)  +cfs(6)/(4*dx^2);
gmat(end,end) = gmat(end,end)+cfs(6)/(4*dx^2);

gmat = fft2(gmat);


V = get_V(X,Y);
%[vecout] = apply_op_fast(xvec,gmat,npts,dx,ckb,V);


[y0] = get_y0(X,Y,V,ckb);

tic
fun_mat = @(x) apply_op_fast(x,gmat,npts,dx,ckb,V);
[xout,FLAG,RELRES,ITER,RESVEC] = gmres(fun_mat,y0(:),200,10^(-9),20);
xout = reshape(xout,[npts,npts]);
xout = xout(:);

params = [];
params.npts = npts;
params.ckb  = ckb;
params.lmax = lmax;
params.dx   = dx;
params.X    = X(:);
params.Y    = Y(:);
params.V    = V(:);
params.Y    = Y(:);

end