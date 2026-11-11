function golden_log(fn, obj, out, errId)
% golden_log  Append one call of an exact-engine entry point to the golden trace (JSON lines).
%   fn    : 'conjQ' | 'biconjQ'
%   obj   : the QuaPol input
%   out   : cell of outputs ({g, dom} for conjQ, {h} for biconjQ), or {} when the call raised
%   errId : '' or the identifier the call raised
% Integers are written with %.0f and the caller's doubles with %.17g, so nothing is rounded on the
% way out; python/cca2/io_matlab.py reads the file back.
    f = getenv('CCA2_GOLDEN_FILE');
    if isempty(f), return, end
    s = struct();
    s.fn = fn;
    s.in = quaPolJson(obj);
    s.err = errId;
    if ~isempty(out)
        if strcmp(fn, 'conjQ')
            s.out = resultJson(out{1});
            s.dom = mat(out{2});
        else
            s.out = resultJson(out{1});
        end
    end
    fid = fopen(f, 'a');
    fprintf(fid, '%s\n', enc(s));
    fclose(fid);
end

function j = quaPolJson(obj)
    j = struct();
    j.nv = obj.nv;  j.ne = obj.ne;  j.nf = obj.nf;
    j.V = matg(obj.V);  j.E = mat(obj.E);  j.f = matg(obj.f);  j.F = mat(obj.F);
    j.exact = obj.isExact();
    if j.exact
        j.fN = mat(obj.fN);  j.fD = mat(obj.fD(:).');  j.VN = mat(obj.VN);  j.VD = mat(obj.VD(:).');
    end
    if ~isempty(obj.fIsConvex), j.fIsConvex = double(obj.fIsConvex); end
end

function j = resultJson(g)
    j = struct();
    j.kind = class(g);
    X = samplePoints();
    switch class(g)
        case 'QuaCon'
            j.nv = g.nv;  j.ne = g.ne;  j.nf = g.nf;
            j.Vname = mat(g.Vname);  j.E = mat(g.E);  j.EcQ = mat(g.EcQ);
            j.fN = mat(g.fN);  j.fD = mat(g.fD(:).');  j.F = mat(g.F);
            fc = cell(1, numel(g.FC));
            for k = 1:numel(g.FC), fc{k} = mat(g.FC{k}); end
            j.FC = ['[' strjoin(fc, ',') ']'];
            [v, idx] = g.eval(X);
            j.evalX = matg(X);  j.evalV = matg(v(:).');  j.evalIdx = mat(idx(:).');
        case 'QuaPol'
            j.pol = quaPolJson(g);
            v = g.eval(X);
            j.evalX = matg(X);  j.evalV = matg(v(:).');
        otherwise
            try
                v = g.eval(X);
                j.evalX = matg(X);  j.evalV = matg(v(:).');
            catch e
                j.evalErr = e.identifier;
            end
            try, j.nf = g.nf; catch, end
    end
end

function X = samplePoints()
    [a, b] = meshgrid([-3.25 -1.5 -0.5 0.125 0.75 1.375 2.5]);
    X = [a(:), b(:)];
end

function s = mat(A)
% an INTEGER matrix as nested JSON arrays, exact
    if isempty(A), s = '[]'; return, end
    rows = cell(1, size(A,1));
    for i = 1:size(A,1)
        rows{i} = ['[' strjoin(arrayfun(@(x) sprintf('%.0f', x), A(i,:), 'UniformOutput', false), ',') ']'];
    end
    s = ['[' strjoin(rows, ',') ']'];
end

function s = matg(A)
% a DOUBLE matrix as nested JSON arrays, round-trip exact (%.17g); Inf/NaN as strings
    if isempty(A), s = '[]'; return, end
    rows = cell(1, size(A,1));
    for i = 1:size(A,1)
        rows{i} = ['[' strjoin(arrayfun(@num, A(i,:), 'UniformOutput', false), ',') ']'];
    end
    s = ['[' strjoin(rows, ',') ']'];
end

function t = num(x)
    if isnan(x), t = '"nan"';
    elseif isinf(x) && x > 0, t = '"inf"';
    elseif isinf(x), t = '"-inf"';
    else, t = sprintf('%.17g', x);
    end
end

function t = enc(s)
% a struct whose char fields are already JSON fragments (from mat/matg) or plain strings
    fn = fieldnames(s);
    parts = cell(1, numel(fn));
    for i = 1:numel(fn)
        v = s.(fn{i});
        if isstruct(v)
            w = enc(v);
        elseif ischar(v)
            if ~isempty(v) && (v(1) == '[' || v(1) == '{'), w = v;
            else, w = ['"' strrep(v, '"', '\"') '"']; end
        elseif islogical(v)
            if v, w = 'true'; else, w = 'false'; end
        else
            w = sprintf('%.17g', v);
        end
        parts{i} = sprintf('"%s":%s', fn{i}, w);
    end
    t = ['{' strjoin(parts, ',') '}'];
end
