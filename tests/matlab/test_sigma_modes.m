function test_sigma_modes(packageDir, legacyFile)
%TEST_SIGMA_MODES Check the public evaluator, fixed-source modes and provenance.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'matlab')); setup_surfsmooth3d;
pkg = surfsmooth3d.edgepreserve.load_package(packageDir);
[work,cleanup] = surfsmooth3d.edgepreserve.stage_source(pkg); %#ok<ASGLU>
settings = surfsmooth3d.edgepreserve.default_settings(); settings.rlam = 4;
targets = pkg.surface.r(:,1:117:end);
opts = surfsmooth3d.edgepreserve.mesher_options(pkg,work,settings,0);
v = pkg.scaffold.nodes; t = pkg.scaffold.tris;
a = v(:,t(1,:)); b = v(:,t(2,:)); c = v(:,t(3,:));
lengths = [vecnorm(b-a);vecnorm(c-b);vecnorm(a-c)];
centers = (a+b+c)/3;
legacy = [];
if nargin > 1 && ~isempty(legacyFile), legacy = load(legacyFile); end
for flag = [0 1 2 3]
    opts.adapt_sigma = flag;
    [sigma,gradient,meta] = surfsmooth3d.multiscale_mesher_sigma_eval(work.scaffold_file,6,opts,targets);
    assert(meta.adapt_sigma == flag && all(isfinite(sigma) & sigma>0));
    assert(all(isfinite(gradient),'all'));
    if flag == 0
        expected = mean(max(lengths,[],1))/settings.rlam;
        assert(max(abs(sigma-expected))<1e-13 && all(gradient==0,'all'));
    elseif flag == 1 || flag == 3
        seeds = max(lengths,[],1)/settings.rlam;
        if flag == 3, seeds = min(lengths,[],1)/settings.rlam; end
        width = sqrt(5/2)*max(seeds);
        delta = reshape(targets,3,[],1)-reshape(centers,3,1,[]);
        weights = reshape(exp(-sum(delta.^2,1)/(2*width^2)),size(targets,2),[]);
        expected = weights*seeds(:)./sum(weights,2);
        assert(max(abs(sigma-expected))<1e-12);
        for axis = 1:3
            dw = -reshape(delta(axis,:,:),size(weights))/width^2.*weights;
            expectedGradient = sum(dw.*(seeds-expected),2)./sum(weights,2);
            assert(max(abs(gradient(axis,:).'-expectedGradient))<1e-12);
        end
    end
    if ~isempty(legacy) && flag<=2
        assert(isequal(targets,legacy.targets));
        assert(isequal(sigma,legacy.sigma(:,flag+1)));
        assert(isequal(gradient,legacy.gradient(:,:,flag+1)));
    end
    for order = [4 8]
        opts.nrefine = 2;
        [sameSigma,sameGradient] = surfsmooth3d.multiscale_mesher_sigma_eval(work.scaffold_file,order,opts,targets);
        assert(isequal(sigma,sameSigma) && isequal(gradient,sameGradient));
    end
    settings.adapt_sigma = flag;
    directory = surfsmooth3d.edgepreserve.new_run_directory(work.directory,pkg,settings,'sigma_test');
    stored = load(fullfile(directory,'parameters.mat'));
    assert(stored.parameters.settings.adapt_sigma == flag);
    assert(strcmp(stored.parameters.sigma_mode,surfsmooth3d.edgepreserve.sigma_mode_label(flag)));
    filename = surfsmooth3d.edgepreserve.surface_filename(pkg.geometry_name,settings,6,0,'fullysmooth');
    assert(contains(filename,sprintf('_adapt%d',flag)));
end
for bad = {-1,4,1.5,nan,inf,[0 1],'1',1+1i}
    opts.adapt_sigma = bad{1};
    rejected(@() surfsmooth3d.multiscale_mesher(work.scaffold_file,4,opts),'MULTISCALE_MESHER:InvalidSigmaMode');
    rejected(@() surfsmooth3d.multiscale_mesher_sigma_eval(work.scaffold_file,4,opts,targets), ...
        'MULTISCALE_MESHER_SIGMA_EVAL:InvalidSigmaMode');
end
% Normalize logical flags before passing them to the double-only MEX gateway.
opts.adapt_sigma = true;
sigmaLogical = surfsmooth3d.multiscale_mesher_sigma_eval(work.scaffold_file,4,opts,targets);
opts.adapt_sigma = 1;
assert(isequal(sigmaLogical,surfsmooth3d.multiscale_mesher_sigma_eval(work.scaffold_file,4,opts,targets)));
for k = 1:numel(pkg.source_paths)
    assert(strcmp(surfsmooth3d.edgepreserve.file_sha256(pkg.source_paths{k}),pkg.source_hashes{k}));
end
fprintf('PASS: public sigma modes, gradients, legacy values, output-resolution independence and provenance.\n');
end

function rejected(call,identifier)
try
    call();
catch exception
    assert(strcmp(exception.identifier,identifier),exception.message);
    return
end
error('sigma_modes:validation','Invalid sigma mode was accepted.');
end
