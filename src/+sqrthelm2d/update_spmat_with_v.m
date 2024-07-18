function [Sout] = update_spmat_with_v(Sin, V)
    Sout = Sin;
    spmat = Sin.spmat;
    ckb = Sin.ckb;
    dx = Sin.dx;
    n = Sin.n;

    spmat = spdiags(-ckb*sqrt(V(:)),0,n,n)*spmat*dx^2*spdiags(sqrt(V(:)),0,n,n)...
              + spdiags(ones(n,1),0,n,n);

    Sout.spmat_with_v = spmat;
end