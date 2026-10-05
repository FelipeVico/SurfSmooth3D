function result = solve(pkg, settings)
%SOLVE Stage immutable inputs and compute one complete adaptive result.
[work,cleanup] = surfsmooth3d.edgepreserve.stage_source(pkg); %#ok<ASGLU>
opts = settings;
opts.fcad = work.cad_file;
opts.nquad = pkg.source_order;
opts.filetype = 3;
[S,info] = surfsmooth3d.multiscale_mesher_adaptive(work.scaffold_file,settings.order,opts);
% Do not expose a deleted temporary source path as persistent provenance.
info.settings = rmfield(info.settings,'fcad');
reference = surfsmooth3d.adaptivesmoother.resample_cad(pkg,S,info);
report = surfsmooth3d.edgepreserve.validity_report(S,'adaptive fully smooth',reference);
result = struct('surface',S,'info',info,'reference',reference,'report',report, ...
    'settings',settings);
end
