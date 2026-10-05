function amat = vals2coefs(norder,uvs)
    bmat = surfsmooth3d.internal.polytens.lege.coefs2vals(norder,uvs);
    amat = inv(bmat);
end