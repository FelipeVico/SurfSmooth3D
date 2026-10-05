function report = test_bie_consumer(bieRoot, filename, geometryScaled)
%TEST_BIE_CONSUMER Solve using only unchanged fmm3dbie and a saved GO3 file.
restoredefaultpath;
addpath(fullfile(bieRoot,'matlab'),fullfile(bieRoot,'matlab','src'));
assert(isempty(which('surfsmooth3d.multiscale_mesher')));
assert(isempty(which('surfsmooth3d_routs')));
S=surfer.load_from_file(filename);
source=struct('r',[.1;.2;-.1]);
targets=struct('r',[2.5 -2 0;.2 .5 -2.5;.1 .3 .4]);
if nargin>=3 && geometryScaled
    % The STEP regression is a convex intersection of two balls. Its center
    % supplies an interior source; radius-scaled targets are well outside it.
    center=mean(S.r,2); radius=max(vecnorm(S.r-center));
    source.r=center+radius*[.01;.02;-.01];
    targets.r=center+radius*targets.r;
end
eps=1e-6; pars=[1 1];
opts=struct('ifout',1,'eps_gmres',1e-7,'maxit',80);
rhs=lap3d.kern(source,S,'s');
[sigma,history,residual]=lap3d.dirichlet.solver(S,rhs,eps,pars,opts);
computed=lap3d.dirichlet.eval(S,sigma,targets,eps,pars);
exact=1./(4*pi*vecnorm(targets.r-source.r,2,1).');
error=norm(computed(:)-exact)/norm(exact);
assert(isfinite(residual) && residual<1e-5);
assert(error<5e-4,'BIE relative error %.4g.',error);
report=struct('nodes',S.npts,'patches',S.npatches,'relative_error',error, ...
    'residual',residual,'iterations',numel(history),'surfsmooth3d_on_path',false, ...
    'source',source.r,'targets',targets.r);
fid=fopen(fullfile(fileparts(filename),'bie-solve-report.json'),'w');
c=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report,PrettyPrint=true));
fprintf('BIE_CONSUMER_PASSED: residual %.4g, relative error %.4g.\n',residual,error);
clear mex;
end
