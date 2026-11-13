function coverage_scope(suites, outName)
% coverage_scope  Line coverage of the production .m files under the END-TO-END exact-engine suites.
% The covered lines are the Python port's true scope (python/PORT_PLAN.md, Phase 0).
%   coverage_scope({'conjQTest','biconjQTest'}, 'exact_scope.xml')
    import matlab.unittest.TestSuite
    import matlab.unittest.TestRunner
    import matlab.unittest.plugins.CodeCoveragePlugin
    import matlab.unittest.plugins.codecoverage.CoberturaFormat
    p = strsplit(path, pathsep); m = p(contains(p, [filesep 'toolbox' filesep 'maple']));
    if ~isempty(m), rmpath(m{:}); addpath(fullfile(matlabroot,'toolbox','symbolic','symbolic')); end
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % CCA2
    src = fullfile(root, 'MATLAB');                                    % the production .m files
    addpath(src);
    s = [];
    for k = 1:numel(suites), s = [s, TestSuite.fromClass(meta.class.fromName(suites{k}))]; end %#ok<AGROW>
    runner = TestRunner.withTextOutput;
    runner.addPlugin(CodeCoveragePlugin.forFolder(src, 'IncludingSubfolders', false, ...
        'Producing', CoberturaFormat(fullfile(root, '.claude', 'coverage', outName))));
    r = runner.run(s);
    fprintf('COUNTS %d %d %d\n', sum([r.Passed]), sum([r.Failed]), sum([r.Incomplete]));
    for k = 1:numel(r), if r(k).Failed, fprintf('FAILED %s\n', r(k).Name); end, end
end
