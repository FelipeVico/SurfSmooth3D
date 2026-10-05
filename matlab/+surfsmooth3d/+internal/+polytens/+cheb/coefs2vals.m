function bmat = coefs2vals(norder,uvs)
    bmat = surfsmooth3d.internal.polytens.cheb.pols(norder,uvs);
    bmat = bmat.';
end