function bmat = coefs2vals(norder,uvs)
    bmat = surfsmooth3d.internal.polytens.lege.pols(norder,uvs);
    bmat = bmat.';
end