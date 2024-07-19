function pmat = prepare_postprocess(S, targs)

    [~, ns] = size(S.xpts);
    [~, nt] = size(targs);

    xs = repmat(S.xpts(1,:), nt, 1);
    ys = repmat(S.xpts(2,:), nt, 1);

    xt = repmat(targs(1,:).', 1, ns);
    yt = repmat(targs(2,:).', 1, ns);
    
    rx = xt - xs;
    ry = yt - ys;
    pmat = sqrthelm2d.green(rx, ry, S.ckb);
   
end