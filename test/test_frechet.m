ckb = 10;
npts = 151;
lmax = pi/2;


% setup targinfo
nt = 100;
thet = 0:2*pi/nt:2*pi-2*pi/nt;
rtarg = 10;
targinfo = rtarg*[cos(thet); sin(thet)];


cfs = sqrthelm2d.get_diag_correction(ckb, lmax, npts);
S = sqrthelm2d.prepare_solver(lmax, npts, ckb, cfs);

V = qfuns.gaussian(S.xpts(1,:), S.xpts(2,:));
V = V(:);
theta_in = pi/3;
y0 = exp(1i*ckb*(S.xpts(1,:)*cos(theta_in) + S.xpts(2,:)*sin(theta_in)));
y0 = y0(:);
y = y0.*V*ckb;

opts = [];
opts.ifflam = false;
tic, xout = sqrthelm2d.solve(S, V, y, opts); toc;


pmat = sqrthelm2d.prepare_postprocess(S, targinfo);
pot = sqrthelm2d.postprocess(S, xout, targinfo, pmat);

dVfun = @(x,y) exp(-((x+0.01).^2 + (y-0.03).^2)/0.25^2);
dV = dVfun(S.xpts(1,:), S.xpts(2,:));
dV = dV(:);
h = 0.01;

VpdV = V + dV*h;
ypdV = y0.*VpdV*ckb;
xoutp = sqrthelm2d.solve(S, VpdV, ypdV, opts); toc;
pot_p = sqrthelm2d.postprocess(S, xoutp, targinfo, pmat);

VmdV = V - dV*h;
ymdV = y0.*VmdV*ckb;
xoutm = sqrthelm2d.solve(S, VmdV, ymdV, opts); toc;
pot_m = sqrthelm2d.postprocess(S, xoutm, targinfo, pmat);

df = (pot_p - pot_m)/2/h;


usc = sqrthelm2d.compute_scattered_field(S, V, xout, opts);
ut = usc + y0;

y1 = ckb*dV(:).*ut;
xout_df = sqrthelm2d.solve(S, V, y1, opts); 
pot_df = sqrthelm2d.postprocess(S, xout_df, targinfo, pmat);
