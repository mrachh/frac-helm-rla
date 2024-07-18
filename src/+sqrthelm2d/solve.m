function [x, Sout] = solve(S, V, y0, opts)
    if nargin < 4
        opts = [];
    end

    ifflam = false;
    tol = 1e-7;
    if isfield(opts, 'ifflam')
        ifflam = opts.ifflam;
    end

    if isfield(opts, 'tol')
        tol = opts.tol;
    end

    ifupdate_spmat = true;
    if isfield(opts, 'ifupdate_spmat')
        ifupdate_spmat = opts.ifupdate_spmat;
    end

    ifcompute_f = true;
    if isfield(opts, 'ifcompute_f')
        ifcompute_f = opts.ifcompute_f;
    end
    if ifflam
        if ifupdate_spmat
            Sout = sqrthelm2d.update_spmat_with_v(S, V);
        end

        if ifcompute_f
            Sout = sqrthelm2d.compute_factorization(Sout, V);
        end

        x = rskelf_sv(Sout.F, y0(:));

    else
        Sout = S;
        fun_mat = @(x) sqrthelm2d.apply_op_fast(x, S.gmat, S.npts, S.dx, S.ckb, V);
        x = gmres(fun_mat, y0(:), [], tol, S.n);
    end
end