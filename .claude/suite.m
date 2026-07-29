% Run the whole suite in $CCA2DIR and print a per-suite pass/fail table plus the totals.
DIR = getenv('CCA2DIR');
cd(DIR);
fprintf('=== suite: %s ===\n', DIR);
d = [dir('*Test.m'); dir('test*.m')];
names = unique(cellfun(@(n) n(1:end-2), {d.name}, 'UniformOutput', false));
tp = 0; tf = 0; ti = 0;
for i = 1:numel(names)
    try
        r = runtests(names{i});
        np = sum([r.Passed]); nf = sum([r.Failed]); ni = sum([r.Incomplete]);
    catch ME
        fprintf('SUITE %-28s ERRORED (%s)\n', names{i}, ME.identifier);
        continue
    end
    tp = tp + np; tf = tf + nf; ti = ti + ni;
    fprintf('SUITE %-28s pass=%3d fail=%3d incomplete=%3d\n', names{i}, np, nf, ni);
    for k = 1:numel(r)
        if r(k).Failed
            fprintf('  FAIL %s\n', r(k).Name);
        end
    end
end
fprintf('TOTALS pass=%d fail=%d incomplete=%d over %d suites\n', tp, tf, ti, numel(names));
