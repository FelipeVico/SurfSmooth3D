function [beta,gradient,dsoft] = beta_values(r,edges,s,sigma,gradSigma,diameter)
%BETA_VALUES Existing edgepreserve formula, with its spatial derivative.
% Clipping is unchanged. NaN derivatives identify nondifferentiable joins.
n = size(r,2);
beta = zeros(1,n); gradient = zeros(3,n); dsoft = inf(1,n);
if isempty(s.selectedEdgeIds), return; end
xyz = zeros(3,0); weights = zeros(1,0);
for id = unique(s.selectedEdgeIds)
    edge = edges(find([edges.edge_id]==id,1));
    for j = 1:numel(edge.panels)
        xyz = [xyz edge.panels(j).xyz]; %#ok<AGROW>
        weights = [weights reshape(edge.panels(j).w,1,[])]; %#ok<AGROW>
    end
end
if isempty(weights), error('adaptiveblend:edges','Empty selected edge quadrature.'); end
sigma = reshape(sigma,1,[]);
floorSigma = max(s.sigmaFloorRel*diameter,realmin);
gradSigma(:,sigma < floorSigma) = 0;
sigma = max(sigma,floorSigma);
if any(~isfinite(sigma) | sigma<=0), error('adaptiveblend:sigma','Invalid sigma.'); end
e = s.c_eps*sigma; ge = s.c_eps*gradSigma;
A = zeros(1,n); gA = zeros(3,n);
chunk = max(1,min(2000,floor(2000000/numel(weights))));
for first = 1:chunk:n
    ix = first:min(n,first+chunk-1);
    dx = r(1,ix).'-xyz(1,:); dy = r(2,ix).'-xyz(2,:); dz = r(3,ix).'-xyz(3,:);
    q = dx.^2+dy.^2+dz.^2;
    el = e(ix).';
    weighted = exp(-q./(2*el.^2)).*weights;
    normer = 1./(sqrt(2*pi)*el);
    A(ix) = (normer.*sum(weighted,2)).';
    spatial = -[sum(weighted.*dx,2) sum(weighted.*dy,2) sum(weighted.*dz,2)]./(el.^2);
    scale = sum(weighted.*(q./el.^3-1./el),2);
    gA(:,ix) = (normer.*spatial).' + (normer.*scale).'.*ge(:,ix);
end
floorA = max(s.distanceFloor,realmin);
L = max(0,-2*log(max(A,floorA))); root = sqrt(L);
dsoft = e.*root;
gd = ge.*root;
moving = A > floorA & A < 1;
gd(:,moving) = gd(:,moving)-gA(:,moving).*(e(moving)./(root(moving).*A(moving)));
% The square-root distance is not differentiable at A=1 in general.
join = abs(A-1)<=32*eps;
gd(:,join) = NaN;
radius = max(s.c_rad*sigma,eps); tau = max(s.c_tau*sigma,eps);
gr = s.c_rad*gradSigma; gt = s.c_tau*gradSigma;
gr(:,s.c_rad*sigma < eps) = 0; gt(:,s.c_tau*sigma < eps) = 0;
z = (dsoft-radius)./(sqrt(2)*tau);
beta = min(max(.5*erfc(z),0),1);
gradient = -exp(-z.^2)/sqrt(pi).*((gd-gr)./(sqrt(2)*tau) - ...
    gt.*((dsoft-radius)./(sqrt(2)*tau.^2)));
end
