function [x, Sout] = compute_scattered_field(S, V, x0, opts)
    if nargin < 4
        opts = [];
    end

    ifflam = false;
    tol = 1e-9;
    if isfield(opts, 'ifflam')
        ifflam = opts.ifflam;
    end

    if isfield(opts, 'tol')
        tol = opts.tol;
    end

    ifupdate_spmat = false;
    if isfield(opts, 'ifupdate_spmat')
        ifupdate_spmat = opts.ifupdate_spmat;
    end

    ifcompute_f = false;
    if isfield(opts, 'ifcompute_f')
        ifcompute_f = opts.ifcompute_f;
    end
    Sout = S;
    if ifflam
        if ifupdate_spmat
            Sout = sqrthelm2d.update_spmat_with_v(S, V);
        end

        if ifcompute_f
            Sout = sqrthelm2d.compute_factorization(Sout, V);
        end
        x02 = x0(:).*sqrt(V(:));
        x = rskelf_mv(Sout.F, x02(:));
        x = x./sqrt(V(:));
        x = -(x - x0)/S.ckb./V(:);
    else 
        x = sqrthelm2d.apply_op_fast(x0, S.gmat, S.npts, S.dx, S.ckb, V);
        x = -(x - x0)/S.ckb./V(:);
    end
end