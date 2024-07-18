function [cvals] = gfunc(xs, ys, ckb)

    rs_k = ckb*sqrt(xs.^2 + ys.^2);
    [h0] = struve101(rs_k);
    cvals = 1./(2*pi*rs_k) + 1i/4*(h0+besselh(0,rs_k));
    cvals = cvals*ckb;
    
end

