function [vecout] = apply_op_fast(vecin, gmat, npts, dx, ckb, V)
    vuse = reshape(V, [npts, npts]);
    svin = size(vecin);
    vecin = reshape(vecin, [npts,npts]);
    vmat = zeros(size(gmat));
    vmat(1:npts,1:npts) = vecin;
    wmat = ifft2(gmat.*fft2(vmat));
    vecout = -ckb*vuse.*wmat(1:npts,1:npts)*dx^2 + vecin;
    vecout = reshape(vecout,svin);
    
end

