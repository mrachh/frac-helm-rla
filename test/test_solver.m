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
y0 = y0.*V*ckb;

opts = [];
opts.ifflam = false;
tic, xout = sqrthelm2d.solve(S, V, y0, opts); toc;


pmat = sqrthelm2d.prepare_postprocess(S, targinfo);
% evaluate potential at all targets
pot = sqrthelm2d.postprocess(S, xout, targinfo, pmat);






%%
opts.ifflam = true;
tic, [xout2, S] = sqrthelm2d.solve(S, V, y0, opts); toc;

%%
opts.ifcompute_f = false;
tic, [xout3] = sqrthelm2d.solve(S, V, y0, opts); toc;