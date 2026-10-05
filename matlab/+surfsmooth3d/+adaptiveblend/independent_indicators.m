function eta = independent_indicators(rec,b,reference,refu,refv,diameter,radius)
%INDEPENDENT_INDICATORS Compare the partial polynomial with an independent map.
rn = cross(refu/radius,refv/radius,1); rjac = vecnorm(rn); rn = rn./rjac;
coeff = (rec.src(1:3,:)-rec.src(1:3,1))*b.transform.';
rpoly = coeff*b.pol+rec.src(1:3,1);
pn = cross(coeff*b.check_du/radius,coeff*b.check_dv/radius,1);
pjac = vecnorm(pn); pn = pn./pjac;
nn = rec.src(10:12,:)*b.transform.'*b.pol;
nlen = vecnorm(nn); nn = nn./nlen;
eta = [max(vecnorm(rpoly-reference))/diameter, ...
    max([vecnorm(pn-rn),vecnorm(nn-rn)])];
invalid = any(~isfinite([reference(:);rpoly(:);rn(:);pn(:);nn(:)])) || ...
    any(rjac<=1e-28 | pjac<=1e-28 | nlen<=1e-14) || ...
    any(sum(pn.*rn,1)<=0 | sum(nn.*rn,1)<=0);
if invalid || any(~isfinite(eta)), eta(:) = inf; end
end
