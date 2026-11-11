function golden_any(fn, args, outs, errId)
% golden_any  Append one call of ANY traced function to the trace file (JSON lines):
%   {"fn":..., "in":[...], "out":[...], "err":...}
% with every MATLAB value serialised by `ser` below. At most CCA2_TRACE_CAP calls per function
% (default 150) are kept, and a call whose record exceeds 400 kB is dropped.
    f = getenv('CCA2_TRACE_FILE');
    if isempty(f), return, end
    persistent count cap
    if isempty(count)
        count = containers.Map('KeyType', 'char', 'ValueType', 'double');
        cap = str2double(getenv('CCA2_TRACE_CAP'));
        if isnan(cap), cap = 150; end
    end
    if ~isKey(count, fn), count(fn) = 0; end
    if count(fn) >= cap, return, end
    try
        s = ['{"fn":"' fn '","in":' ser(args) ',"out":' ser(outs) ',"err":"' errId '"}'];
    catch e
        s = ['{"fn":"' fn '","serr":"' strrep(e.message, '"', '''') '"}'];
    end
    if numel(s) > 400000, return, end
    count(fn) = count(fn) + 1;
    fid = fopen(f, 'a');
    fprintf(fid, '%s\n', s);
    fclose(fid);
end

function s = ser(x)
    if isa(x, 'bigZ')
        c = x.c;
        v = cell(1, numel(c));
        for i = 1:numel(c), v{i} = ['"' char(c{i}.toString()) '"']; end
        s = ['{"t":"z","sz":' szs(size(c)) ',"v":[' strjoin(v, ',') ']}'];
    elseif isjava(x)
        s = ['{"t":"j","v":"' char(x.toString()) '"}'];
    elseif isnumeric(x) || islogical(x)
        x = double(x);
        v = cell(1, numel(x));
        for i = 1:numel(x), v{i} = num(x(i)); end
        s = ['{"t":"n","sz":' szs(size(x)) ',"v":[' strjoin(v, ',') ']}'];
    elseif ischar(x)
        s = ['"' strrep(strrep(x(:).', '\', '\\'), '"', '\"') '"'];
    elseif isstring(x)
        s = ser(char(x));
    elseif iscell(x)
        v = cell(1, numel(x));
        for i = 1:numel(x), v{i} = ser(x{i}); end
        s = ['{"t":"c","sz":' szs(size(x)) ',"v":[' strjoin(v, ',') ']}'];
    elseif isstruct(x)
        fn = fieldnames(x);
        v = cell(1, numel(x));
        for i = 1:numel(x)
            p = cell(1, numel(fn));
            for k = 1:numel(fn), p{k} = ['"' fn{k} '":' ser(x(i).(fn{k}))]; end
            v{i} = ['{' strjoin(p, ',') '}'];
        end
        s = ['{"t":"s","sz":' szs(size(x)) ',"v":[' strjoin(v, ',') ']}'];
    elseif isa(x, 'function_handle')
        s = ['{"t":"f","v":"' strrep(func2str(x), '"', '''') '"}'];
    elseif isobject(x)
        pn = properties(x);
        p = cell(1, numel(pn));
        for k = 1:numel(pn)
            try, p{k} = ['"' pn{k} '":' ser(x.(pn{k}))];
            catch, p{k} = ['"' pn{k} '":null']; end
        end
        s = ['{"t":"o","cls":"' class(x) '","p":{' strjoin(p, ',') '}}'];
    else
        s = 'null';
    end
end

function s = szs(sz)
    s = ['[' strjoin(arrayfun(@(v) sprintf('%d', v), sz, 'UniformOutput', false), ',') ']'];
end

function t = num(x)
    if isnan(x), t = '"nan"';
    elseif isinf(x) && x > 0, t = '"inf"';
    elseif isinf(x), t = '"-inf"';
    elseif x == fix(x) && abs(x) < 1e15, t = sprintf('%.0f', x);
    else, t = sprintf('%.17g', x);
    end
end
