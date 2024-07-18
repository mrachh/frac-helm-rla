function [spmat] = get_sparse_correction(n, dx, ckb, cfs)
    e = ones(n,1);
    A = spdiags([e e], [-1,1], n, n);
    D = spdiags(e, 0, n, n);
    spmat = kron(A,D) + kron(D,A);
    m = size(spmat,1);
    spmat = spmat*(cfs(2)/(2*dx^2)) + spdiags(ones(m,1),0,m,m)*(cfs(1)/(dx^2));
end