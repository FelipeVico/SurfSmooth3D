function [S_cad, mapping] = resample_cad(pkg, S_output, outputLevel)
%RESAMPLE_CAD Evaluate one fixed piecewise CAD polynomial on an output layout.
order = double(S_output.norders(1));
if any(S_output.norders(:) ~= order) || any(S_output.iptype(:) ~= 1) || ...
        S_output.npatches ~= pkg.coarse_triangle_count*4^outputLevel
    error('edgepreserve:outputLayout', 'Output is not the requested uniform RV layout.');
end
mapping = surfsmooth3d.edgepreserve.reference_map(pkg.source_level,outputLevel,order);
source = pkg.surface;
nsource = 4^pkg.source_level;
nperparent = numel(mapping.source_patch);
r = zeros(3,S_output.npts);
evaluation = cell(1,nsource);
indices = cell(1,nsource);
for child = 1:nsource
    indices{child} = find(mapping.source_patch == child);
    if ~isempty(indices{child})
        evaluation{child} = surfsmooth3d.internal.koorn.pols(pkg.source_order, ...
            mapping.source_uv(:,indices{child}));
    end
end
for parent = 1:pkg.coarse_triangle_count
    for child = 1:nsource
        if isempty(indices{child}), continue; end
        patch = (parent-1)*nsource+child;
        destination = (parent-1)*nperparent+indices{child};
        r(:,destination) = source.srccoefs{patch}(1:3,:)*evaluation{child};
    end
end
srcvals = surfsmooth3d.edgepreserve.srcvals_from_positions(S_output,r);
S_cad = surfsmooth3d.surfer(S_output.npatches,S_output.norders,srcvals,S_output.iptype);
end
