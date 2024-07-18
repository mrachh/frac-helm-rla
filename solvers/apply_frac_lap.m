function [vecout] = apply_frac_lap(vecin,gmat,npts,dx)
    svin = size(vecin);
    vecin = reshape(vecin,[npts,npts]);
    vmat = zeros(size(gmat));
    vmat(1:npts,1:npts) = vecin;
    wmat = ifft2(gmat.*fft2(vmat));
    vecout = wmat(1:npts,1:npts)*dx^2;
    vecout = reshape(vecout,svin);
end

