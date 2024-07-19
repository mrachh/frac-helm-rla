
ckb = pi/2;
lmax= 10;
npts= 51;

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


theta_in = pi/3;
y0 = ckb*exp(1i*X*cos(theta_in)*ckb+1i*Y*sin(theta_in)*ckb).*V;
tic
fun_mat = @(x) apply_op_fast(x,gmat,npts,dx,ckb,V);
[xout,FLAG,RELRES,ITER,RESVEC] = gmres(fun_mat,y0(:),200,10^(-9),20);
xout = reshape(xout,[npts,npts]);
toc

%%
sz = size(X);
X = X(:);
Y = Y(:);

%%
V2 = V.';
V2 = V;

n = numel(xs);
e = ones(n,1);
D = spdiags([e], [0], n, n);

x0 = kron(D,D);

xp1 = kron(spdiags(e,1,n,n),D);
xp2 = kron(spdiags(e,2,n,n),D);
xm1 = kron(spdiags(e,-1,n,n),D);
xm2 = kron(spdiags(e,-2,n,n),D);
yp1 = kron(D,spdiags(e,1,n,n));
yp2 = kron(D,spdiags(e,2,n,n));
ym1 = kron(D,spdiags(e,-1,n,n));
ym2 = kron(D,spdiags(e,-2,n,n));

xp1yp1 = kron(spdiags(e,1,n,n),spdiags(e,1,n,n));
xm1yp1 = kron(spdiags(e,-1,n,n),spdiags(e,1,n,n));
xp1ym1 = kron(spdiags(e,1,n,n),spdiags(e,-1,n,n));
xm1ym1 = kron(spdiags(e,-1,n,n),spdiags(e,-1,n,n));

A2 = xp1*cfs(2)/(2*dx^2)+xm1*cfs(2)/(2*dx^2) + ...
    xp2*cfs(3)/(2*dx^2)+xm2*cfs(3)/(2*dx^2) + ...
    yp1*cfs(4)/(2*dx^2)+ym1*cfs(4)/(2*dx^2) + ...
    yp2*cfs(5)/(2*dx^2)+ym2*cfs(5)/(2*dx^2) + ...
    xp1yp1*cfs(6)/(4*dx^2)+xm1yp1*cfs(6)/(4*dx^2) + ...
    xp1ym1*cfs(6)/(4*dx^2)+xm1ym1*cfs(6)/(4*dx^2) + ...
    x0*cfs(1)/(dx^2);

nl = size(A2,1);
A2 = -ckb*spdiags(sqrt(V2(:)),0,nl,nl)*A2*dx^2*spdiags(sqrt(V2(:)),0,nl,nl) + ...
    spdiags(ones(nl,1),0,nl,nl);
A2 = A2.';

%%
p = 128;
[proxy,pnorm,pw] = proxy_circ_pts(p);

xpts = [X.';Y.'];

pfun = @(xfun,slf,nbr,l,ctr) pxyfun_frac(xfun,slf,nbr,proxy,pw,...
    l,ctr,sqrt(V(:)),ckb,dx);

occ = 4000;
rank_or_tol = 1E-10;
pxyfun = [];
opts = [];
opts.verb = 1;

matind = @(i,j) frac_mat_fun(i,j,X,Y,sqrt(V(:)),ckb,dx,A2);
  F = rskelf(matind,xpts,occ,rank_or_tol,pfun,opts); 

  %%
 y02 = ckb*exp(1i*X*cos(theta_in)*ckb+1i*Y*sin(theta_in)*ckb).*sqrt(V(:));
 sig = rskelf_sv(F,y02(:));
 sig = sig.*sqrt(V(:));

%% 

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
gf = gfunc(dx1,dx2,ckb);
%Kpxy1 = gf*dx;
Kpxy2 = gf*diag(ckb*(vin(slf)));
Kpxy2 = dx*Kpxy2.*sqrt(wt);
Kpxy = [Kpxy2];
dx = x(1,nbr) - ctr(1);
dy = x(2,nbr) - ctr(2);
dist = sqrt(dx.^2 + dy.^2);
nbr = nbr(dist/norm(l) < 1.5);
end


%%
function [amat] = frac_mat_fun(i,j,x,y,v,ckb,dx,scorr)
    ni = numel(i);
    nj = numel(j);
    if (ni*nj == 0)
        amat = zeros(ni,nj);
        return
    end
    [XS,XT] = ndgrid(x(i),x(j));
    [YS,YT] = ndgrid(y(i),y(j));
    dx1 = XS-XT;
    dx2 = YS-YT;
    amat = gfunc(dx1,dx2,ckb);
    amat(isnan(amat)) = 0;
    vmat = (-ckb*dx^2*v(i))*v(j).';
    amat = vmat.*amat;
    amat = amat + scorr(i,j);

end

%%


function [proxy,pnorm,pw] = proxy_circ_pts(p)

if nargin < 1
    p = 64;
end

proxy = [];
pnorm = [];
pw    = [];

for ii=1:5
theta = (0:(p-1))*2*pi/p;
rad = (1.5+(ii-1)*0.5/4);
proxyt = rad*[cos(theta);sin(theta)];
pnormt = [cos(theta); sin(theta)];
pwt = ones(p,1)*(2.0*pi*rad)/p;

proxy = [proxy,proxyt];
pnormt = [pnorm,pnormt];
pw   = [pw,pwt];
end 

end
