function [cfs] = get_corrs(ckb,lmax,npts)
xs = linspace(-lmax,lmax,npts);
ys = xs;
[X,Y] = ndgrid(xs,ys);

xvec = exp(-X.^2-Y.^2);
gval = gfunc(X,Y,ckb);
gval(isnan(gval)) = 0;



dx = xs(2)-xs(1);
qquad = sum(gval(:).*xvec(:))*dx^2;
qadap = integral2(@(x,y) gfunc(x,y,ckb).*exp(-x.^2-y.^2),-10,10,-10,10);
dcorr = (qadap-qquad)/dx^2;

qquadx = sum(gval(:).*xvec(:).*(X(:).^2))*dx^2;
qadapx = integral2(@(x,y) gfunc(x,y,ckb).*exp(-x.^2-y.^2).*(x.^2),-10,10,-10,10);
qquady = sum(gval(:).*xvec(:).*(Y(:).^2))*dx^2;
qadapy = integral2(@(x,y) gfunc(x,y,ckb).*exp(-x.^2-y.^2).*(y.^2),-10,10,-10,10);

xts = [0,dx,0 ];
yts = [0,0,dx];

f0 = exp(-xts.^2-yts.^2);
f1 = f0.*(xts.^2);
f2 = f0.*(yts.^2);

cfs = [f0.',f1.',f2.'].'\[qadap-qquad;qadapx-qquadx;qadapy-qquady];

gval((npts+1)/2,(npts+1)/2) = cfs(1)/(dx^2);
gval((npts+1)/2+1,(npts+1)/2) = gval((npts+1)/2+1,(npts+1)/2)+cfs(2)/(2*dx^2);
gval((npts+1)/2-1,(npts+1)/2) = gval((npts+1)/2-1,(npts+1)/2)+cfs(2)/(2*dx^2);
gval((npts+1)/2,(npts+1)/2+1) = gval((npts+1)/2,(npts+1)/2+1)+cfs(3)/(2*dx^2);
gval((npts+1)/2,(npts+1)/2-1) = gval((npts+1)/2,(npts+1)/2-1)+cfs(3)/(2*dx^2);

qquad = sum(gval(:).*xvec(:))*dx^2;
qadap = integral2(@(x,y) gfunc(x,y,ckb).*exp(-x.^2-y.^2),-10,10,-10,10);
qquadx = sum(gval(:).*xvec(:).*(X(:).^2))*dx^2;
qadapx = integral2(@(x,y) gfunc(x,y,ckb).*exp(-x.^2-y.^2).*(x.^2),-10,10,-10,10);
qquady = sum(gval(:).*xvec(:).*(Y(:).^2))*dx^2;
qadapy = integral2(@(x,y) gfunc(x,y,ckb).*exp(-x.^2-y.^2).*(y.^2),-10,10,-10,10);
[qadap-qquad;qadapx-qquadx;qadapy-qquady]

qquadt = sum(gval(:).*xvec(:).*(sin(X(:)*2.1)+cos(Y(:)*3.4)+1.31))*dx^2;
qadapt = integral2(@(x,y) gfunc(x,y,ckb).*exp(-x.^2-y.^2).*(sin(x*2.1)+cos(y*3.4)+1.31),-10,10,-10,10);
qquadt-qadapt
end
