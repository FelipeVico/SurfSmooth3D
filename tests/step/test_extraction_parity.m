function report = test_extraction_parity(baselineFile, outputDirectory, packageRoot)
%TEST_EXTRACTION_PARITY Replay frozen V3 MATLAB results with SurfSmooth3D.
% The baseline is captured in a separate process from the original package.
% Numerical arrays, diagnostics, counters and exported GO3 data are compared.
% Input/helper paths are the only excluded metadata fields.
% packageRoot optionally selects an isolated build containing matlab/; STEP
% fixtures still come from this repository. Use separate MATLAB processes
% when comparing different native builds.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if nargin < 3 || isempty(packageRoot), packageRoot = root; end
packageRoot = canonical_path(packageRoot);
matlabRoot = fullfile(packageRoot,'matlab');
addpath(matlabRoot,'-begin');
resolved = verify_package_route(matlabRoot);
reference = load(baselineFile,'cases');
if nargin < 2 || isempty(outputDirectory)
    outputDirectory = fullfile(root,'build','validation','step');
end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
report = struct('passed',true,'cases',[], 'package_root',packageRoot, ...
    'fixture_root',root,'baseline_file',canonical_path(baselineFile), ...
    'resolved_functions',resolved,'loaded_mex','', 'route_verified',false);
for k=1:numel(reference.cases)
    c=reference.cases{k};
    mesh=surfsmooth3d.stepmesher.mesh_step(fullfile(root,'step_mesher','examples', ...
        'step_files',[c.name,'.step']),c.options);
    report.loaded_mex = verify_loaded_mex(resolved.expected_mex);
    report.route_verified = true;
    [srcvals,norders,iptype]=surfsmooth3d.stepmesher.to_srcvals(mesh);
    stats=surfsmooth3d.stepmesher.validate_area_flux(mesh);
    [count,exact,maxScaled]=compare_struct(c.mesh,mesh,'mesh');
    [n,e,d]=compare_struct(c.stats,stats,'stats');
    count=count+n; exact=exact+e; maxScaled=max(maxScaled,d);
    [e,d]=compare_array(c.srcvals,srcvals,'srcvals');
    count=count+1; exact=exact+e; maxScaled=max(maxScaled,d);
    assert(isequal(c.norders,norders) && isequal(c.iptype,iptype));
    packageDir=fullfile(outputDirectory,c.name);
    package=surfsmooth3d.stepmesher.export_geom_package(mesh,packageDir,c.name);
    loaded=surfsmooth3d.surfer.load_from_file(package.surface_file);
    [values,~,~,~,~,~]=extract_arrays(loaded);
    [e,d]=compare_array(srcvals,values,'GO3 round trip');
    count=count+1; exact=exact+e; maxScaled=max(maxScaled,d);
    entry=struct('name',c.name,'patches',mesh.npatches,'area',stats.area, ...
        'numeric_comparisons',count,'exact_comparisons',exact, ...
        'maximum_scaled_difference',maxScaled,'output',packageDir);
    report.cases=[report.cases,entry]; %#ok<AGROW>
    fprintf('STEP_PARITY %s: %d/%d exact; max scaled difference %.3g\n', ...
        c.name,exact,count,maxScaled);
end
fid=fopen(fullfile(outputDirectory,'parity.json'),'w');
assert(fid>=0); cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('STEP_EXTRACTION_PARITY_PASSED\n');
end

function resolved=verify_package_route(matlabRoot)
names={'surfsmooth3d.stepmesher.mesh_step', ...
    'surfsmooth3d.stepmesher.to_srcvals', ...
    'surfsmooth3d.stepmesher.validate_area_flux', ...
    'surfsmooth3d.stepmesher.export_geom_package','surfsmooth3d.surfer'};
resolved=struct('names',{names},'paths',{cell(size(names))}, ...
    'expected_mex',canonical_path(fullfile(matlabRoot,'+surfsmooth3d', ...
        '+stepmesher','private',['step_mesher_mex.',mexext])));
for k=1:numel(names)
    found=which(names{k});
    assert(~isempty(found),'stepmesher:parityRoute','Missing %s.',names{k});
    found=canonical_path(found);
    assert(startsWith(found,[matlabRoot,filesep]),'stepmesher:parityRoute', ...
        '%s resolves outside requested build: %s',names{k},found);
    resolved.paths{k}=found;
end
assert(isfile(resolved.expected_mex),'stepmesher:parityMex', ...
    'Requested build has no STEP MEX: %s',resolved.expected_mex);
end

function loaded=verify_loaded_mex(expected)
[~,files]=inmem('-completenames');
loaded='';
for k=1:numel(files)
    [~,stem]=fileparts(files{k});
    if strcmp(stem,'step_mesher_mex')
        found=canonical_path(files{k});
        assert(strcmp(found,expected),'stepmesher:parityMex', ...
            'STEP MEX loaded from the wrong build: %s',found);
        loaded=found;
    end
end
assert(~isempty(loaded),'stepmesher:parityMex', ...
    'Requested STEP MEX was not found among loaded native modules.');
end

function value=canonical_path(value)
value=char(java.io.File(char(value)).getCanonicalPath());
end

function [count,exact,maxScaled]=compare_struct(expected,actual,label)
count=0; exact=0; maxScaled=0;
fields=fieldnames(expected);
assert(isequal(sort(fields),sort(fieldnames(actual))), ...
    'stepmesher:parityFields','%s fields changed',label);
for j=1:numel(fields)
    name=fields{j};
    if any(strcmp(name,{'step_file','occ_info_exe'})), continue; end
    a=expected.(name); b=actual.(name); tag=[label,'.',name];
    if isstruct(a)
        [n,e,d]=compare_struct(a,b,tag);
        count=count+n; exact=exact+e; maxScaled=max(maxScaled,d);
    elseif isnumeric(a) || islogical(a)
        if islogical(a) || isinteger(a) || is_discrete_field(name)
            assert(isequal(a,b),'stepmesher:parityDiscrete', ...
                '%s discrete topology/order/counter data changed',tag);
        end
        [e,d]=compare_array(a,b,tag);
        count=count+1; exact=exact+e; maxScaled=max(maxScaled,d);
    else
        assert(isequal(a,b),'stepmesher:parityMetadata','%s changed',tag);
    end
end
end

function yes=is_discrete_field(name)
yes=any(strcmp(name,{'nodes_per_patch','npatches','norder','iptype', ...
    'parent_coarse_tri','scaffold_triangles','edge_order','edge_nodes_per_panel', ...
    'ncad_edges','npanels_total','edge_ids','edge_adjacent_faces','edge_panel_start', ...
    'native_identity_faces','native_identity_occurrences','refinement_level', ...
    'edge_gll_order','curvature_points','worst_quality_triangle'})) || ...
    endsWith(name,{'_tags','_count','_success','_failures','_fallback_points'});
end

function [exact,scaled]=compare_array(a,b,label)
assert(isequal(size(a),size(b)),'stepmesher:parityShape','%s shape changed',label);
exact=double(isequaln(a,b));
scaled=0;
if exact, return; end
assert(isequal(isfinite(a),isfinite(b)) && isequal(isnan(a),isnan(b)), ...
    'stepmesher:parityFinite','%s finite-value pattern changed',label);
use=isfinite(a);
if any(use(:))
    delta=max(abs(double(a(use))-double(b(use))));
    scale=max(1,max(abs(double(a(use)))));
    scaled=delta/scale;
end
% This is a roundoff-level packaging comparison, not a geometric error bound.
assert(scaled<=1e-12,'stepmesher:parityNumeric', ...
    '%s changed beyond roundoff: scaled difference %.17g',label,scaled);
end
