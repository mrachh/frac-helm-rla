function pot = postprocess(S, x, targs, pmat)
    pot = (pmat*x)*S.dx*S.dx;
end