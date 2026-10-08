function validate_settings(s,pkg)
if isfield(s,'newton_recovery')
    surfsmooth3d.validate_newton_recovery(s.newton_recovery);
end
validateattributes(s.order,{'double'},{'scalar','integer','>=',1,'<=',20});
validateattributes(s.rlam,{'double'},{'scalar','finite','positive'});
validateattributes(s.adapt_sigma,{'double'},{'scalar','integer','>=',0,'<=',3});
validateattributes(s.eps_adapt,{'double'},{'scalar','finite','>',0,'<',1});
validateattributes(s.max_refine,{'double'},{'scalar','integer','>=',0,'<=',20});
validateattributes(s.max_points,{'double'},{'scalar','integer','positive','<=',flintmax});
for name = {'c_eps','c_rad','c_tau','sigmaFloorRel','distanceFloor'}
    validateattributes(s.(name{1}),{'double'},{'scalar','finite','positive'});
end
if s.distanceFloor >= 1, error('adaptiveblend:floor','distanceFloor must be below one.'); end
if ~isempty(s.selectedEdgeIds)
    validateattributes(s.selectedEdgeIds,{'double'},{'vector','integer','positive'});
end
if any(~ismember(s.selectedEdgeIds,[pkg.edges.edge_id]))
    error('adaptiveblend:edgeID','Unknown selected CAD edge ID.');
end
end
