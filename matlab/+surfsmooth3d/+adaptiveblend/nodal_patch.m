function rec = nodal_patch(r,cad,beta,sigma,b,diameter,radius)
%NODAL_PATCH Keep invalid patches assessable rather than aborting refinement.
c = (r-r(:,1))*b.transform.';
du = c*b.du; dv = c*b.dv;
crossn = cross(du/radius,dv/radius,1);
jac = vecnorm(crossn); normal = crossn./jac;
src = [r;du;dv;normal];
eta = [inf inf -1 -1];
if all(isfinite(src(:))) && all(jac>1e-28) && diameter>0
    nc = normal*b.transform.';
    eta(1) = sqrt(2*sum(abs(c(:,b.position_first:end)).^2,'all'))/diameter;
    eta(2) = sqrt(2*sum(abs(nc(:,b.normal_first:end)).^2,'all'));
end
rec = struct('src',src,'cad',cad,'beta',beta,'sigma',sigma, ...
    'indicators',eta,'checked',false);
end
