% 5b: does the nearest-root probe ever pick a DIFFERENT branch than "first real root"?
% Run over every conic constraint x vertex the three testRegion fixtures and the captured merge
% operands present -- the same shapes removeTangent actually sees.
cd(getenv('CCA2DIR')); addpath(getenv('SPDIR'));
warning('off','symbolic:sym:isAlways:TruthUnknown');
D = dir(fullfile(getenv('CCA2_DUMP_MERGE'), 'mg_*.mat'));
regs = {};
for i = 1:numel(D)
    S = load(fullfile(D(i).folder, D(i).name));
    regs{end+1} = S.A; regs{end+1} = S.B; %#ok<AGROW>
end
x = sym('x'); y = sym('y');
regs{end+1} = region([16*x - 4*x*y - x^2 - 4*y^2, -x-2*y, (2*y)/3 - x - sym(1)/3, x + 2*y - 2], [x y]);
regs{end+1} = region([(x+y)^2 - 4*x, -x, -y], [x y]);
regs{end+1} = region([x^2 - y, y - 1, -x - 1], [x y]);

nTot = 0; nDiff = 0; nOneRoot = 0;
for r = 1:numel(regs)
    R = regs{r}; v = R.vars;
    for i = 1:size(R.ineqs,2)
        if ~R.ineqs(i).isQuad, continue, end
        for j = 1:R.nv
            if ~R.ineqs(i).subsF(v,[R.vx(j),R.vy(j)]).isZero, continue, end
            sx = R.vx(j) + 0.1;
            p = R.ineqs(i).subsF(v(1), sx);
            rr = region.rootsIn(p.f, v(2));
            real0 = [];
            for k = 1:numel(rr)
                try, d = double(rr(k)); catch, continue, end
                if isreal(d) && isfinite(d), real0(end+1) = d; end %#ok<AGROW>
            end
            if isempty(real0), continue, end
            nTot = nTot + 1;
            if numel(real0) == 1, nOneRoot = nOneRoot + 1; continue, end
            first = real0(1);
            [~, w] = min(abs(real0 - double(R.vy(j))));
            near = real0(w);
            if abs(first - near) > 1e-12
                nDiff = nDiff + 1;
                fprintf('  reg %d ineq %d v%d: py=%.4f  firstReal=%.4f  nearest=%.4f  gap=%.4f\n', ...
                        r, i, j, double(R.vy(j)), first, near, abs(first-near));
            end
        end
    end
end
fprintf('\nvertex/conic probes %d  (single-root %d)  branch differs %d\n', nTot, nOneRoot, nDiff);
