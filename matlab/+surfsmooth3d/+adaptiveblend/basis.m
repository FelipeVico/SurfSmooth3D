function b = basis(p)
b.rv = surfsmooth3d.internal.koorn.rv_nodes(p);
b.transform = surfsmooth3d.internal.koorn.vals2coefs(p,b.rv);
[~,b.du,b.dv] = surfsmooth3d.internal.koorn.ders(p,b.rv);
children = surfsmooth3d.edgepreserve.triangle_maps(1,'smoother');
b.check = zeros(2,0);
for k = 1:4
    b.check = [b.check children(:,1,k)+children(:,2:3,k)*b.rv];
end
t = (surfsmooth3d.internal.polytens.lege.pts(p+2)+1)/2; t = reshape(t,1,[]);
b.check = [b.check [t 1-t zeros(size(t)); zeros(size(t)) t 1-t] [0 1 0;0 0 1]];
[b.pol,b.check_du,b.check_dv] = surfsmooth3d.internal.koorn.ders(p,b.check);
b.np = size(b.rv,2);
b.nc = size(b.check,2);
b.position_first = max(2,p-1)*(max(2,p-1)+1)/2+1;
b.normal_first = max(1,p-1)*(max(1,p-1)+1)/2+1;
end
