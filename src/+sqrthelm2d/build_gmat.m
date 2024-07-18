function [gmat] = build_gmat(lmax, dx, ckb, cfs)

    nbig = 2^ceil(log(ceil((4*lmax)/dx))/log(2));
    xlrg = (((-nbig)/2+1):nbig/2)*dx;
    cent = find(xlrg == 0);
    xlrg = circshift(xlrg, cent+1);
    
    [XL,YL] = ndgrid(xlrg,xlrg);
    gmat = sqrthelm2d.green(XL, YL, ckb);
    gmat(isnan(gmat)) = 0;
    
    gmat(1,1) = cfs(1)/(dx^2);
    gmat(2,1) = gmat(2,1) + cfs(2)/(2*dx^2);
    gmat(end,1)= gmat(end,1) + cfs(2)/(2*dx^2);
    gmat(1,2) = gmat(1,2) + cfs(3)/(2*dx^2);
    gmat(1,end) = gmat(1,end) + cfs(3)/(2*dx^2);
    
    gmat = fft2(gmat);

end