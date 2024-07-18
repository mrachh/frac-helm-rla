ckb = 2*pi;
lmax= 10;
rs = [1,50;50,1];
npts = [101,201,401];

vouts = zeros(size(rs,1),numel(npts));

for ii=1:numel(npts)
    ii
    npt = npts(ii);
    [xout,params] = solve_it(ckb,npt,lmax);
    vouts(:,ii) = post_proc(rs,xout,params);
end



function [vout] = post_proc(rs,sig,params)

[x1,x2] = ndgrid(rs(1,:),params.X);
[y1,y2] = ndgrid(rs(2,:),params.Y);
xdiff = x1-x2;
ydiff = y1-y2;
gf = gfunc(xdiff,ydiff,params.ckb);

vout = gf*sig*(params.dx^2);


end