
ckb = 2*pi;
lmax= 10;
npts= 101;

xs = linspace(-lmax, lmax, npts);
ys = xs;
[X,Y] = ndgrid(xs,ys);
dx = xs(2)-xs(1);

cfs = sqrthelm2d.get_diag_correction(ckb, lmax, npts);
%%

Sprecomp = sqrthelm2d.prepare_solver(lmax, npts, dx, ckb, cfs);

V = qfuns.gaussian(X,Y);
theta_in = pi/3;
y0 = ckb*exp(1i*X*cos(theta_in)*ckb+1i*Y*sin(theta_in)*ckb).*V;

opts = [];
opts.ifflam = false;
tic, xout = sqrthelm2d.solve(Sprecomp, V, y0, opts); toc;

opts.ifflam = true;
tic, [xout2, Sprecomp] = sqrthelm2d.solve(Sprecomp, V, y0, opts); toc;

%%
opts.ifcompute_f = false;
tic, [xout3] = sqrthelm2d.solve(Sprecomp, V, y0, opts); toc;