ckb = pi/2;
lmax= 10;
npts= 51;

[xout,  ~] = solve_it(ckb,npts,lmax);
[xout2, ~] = solve_flam(ckb,npts,lmax);
norm(xout-xout2)