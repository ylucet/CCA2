classdef bigZ
% bigZ  An array of exact integers in java.math.BigInteger, with MATLAB's arithmetic operators.
%
% WHY IT EXISTS (2026-10-11). The exact pipeline computes integer formulas in doubles and checks
% each result against 2^53 (ratQ.chk). Erratum Example 2 overflowed at one site after another --
% the implicit polynomial, the fibre, the witness quartic, the pulled-back tie -- and every time the
% RAW result was large while its content was huge, so the REDUCED result fit easily. Writing a
% BigInteger copy of each formula by hand is a second spelling of it per site. Instead a formula is
% written once, for doubles, and ratQ.exactFormula re-runs THE SAME CODE on bigZ inputs when the
% double run overflows, then divides out the content and comes back.
%
% Supported: + - .* * (scalar or matrix) .^ ^ (nonnegative integer exponent), unary minus,
% () indexing and assignment, [ , ] and [ ; ], conv, sum, size/numel/length/isempty/end,
% fliplr/flipud/reshape/transpose. Mixed with a double, the double is promoted (it must be an exact
% integer). Anything else raises -- a formula that needs more is not one this class runs.

    properties (SetAccess = private)
        c      % cell array of java.math.BigInteger, the shape of the array
    end

    methods
        function obj = bigZ(A)
            if nargin == 0, obj.c = cell(0,0); return, end
            if isa(A, 'bigZ'), obj = A; return, end
            if iscell(A), obj.c = A; return, end
            if ~isnumeric(A) && ~islogical(A)
                error('bigZ:type', 'bigZ takes numeric input; got %s.', class(A));
            end
            A = double(A);
            if any(A(:) ~= fix(A(:))) || any(~isfinite(A(:)))
                error('bigZ:notInteger', 'bigZ holds exact integers only.');
            end
            if any(abs(A(:)) >= 2^53)
                error('bigZ:inexactInput', 'a double past 2^53 is not an exact integer.');
            end
            obj.c = cell(size(A));
            for i = 1:numel(A), obj.c{i} = java.math.BigInteger.valueOf(int64(A(i))); end
        end

        % ---- shape ---------------------------------------------------------------------------
        function varargout = size(a, varargin), [varargout{1:max(nargout,1)}] = size(a.c, varargin{:}); end
        function n = numel(a, varargin), n = numel(a.c); end
        function n = length(a), n = length(a.c); end
        function tf = isempty(a), tf = isempty(a.c); end
        function e = end(a, k, n)
            sz = size(a.c);
            if n == 1, e = numel(a.c); elseif k <= numel(sz), e = sz(k); else, e = 1; end
        end
        function b = fliplr(a), b = bigZ(fliplr(a.c)); end
        function b = flipud(a), b = bigZ(flipud(a.c)); end
        function b = reshape(a, varargin), b = bigZ(reshape(a.c, varargin{:})); end
        function b = transpose(a), b = bigZ(a.c.'); end
        function b = ctranspose(a), b = bigZ(a.c.'); end

        % ---- indexing ------------------------------------------------------------------------
        function varargout = subsref(a, s)
            switch s(1).type
                case '()'
                    b = bigZ(a.c(s(1).subs{:}));
                    if numel(s) > 1, [varargout{1:nargout}] = subsref(b, s(2:end));
                    else, varargout{1} = b; end
                otherwise
                    [varargout{1:nargout}] = builtin('subsref', a, s);
            end
        end
        function a = subsasgn(a, s, v)
            if ~strcmp(s(1).type, '()') || numel(s) > 1
                error('bigZ:subsasgn', 'only a(i) = v assignment is supported.');
            end
            v = bigZ(v);
            C = subsasgn(a.c, s(1), v.c);           % a 1x1 cell broadcasts, as for any cell
            for i = 1:numel(C)                      % slots created by growth are zero
                if isempty(C{i}), C{i} = java.math.BigInteger.ZERO; end
            end
            a.c = C;
        end
        function n = numArgumentsFromSubscript(~, ~, ~), n = 1; end

        % ---- concatenation -------------------------------------------------------------------
        function b = horzcat(varargin)
            cs = cellfun(@(x) bigZ(x).c, varargin, 'UniformOutput', false);
            b = bigZ(horzcat(cs{:}));
        end
        function b = vertcat(varargin)
            cs = cellfun(@(x) bigZ(x).c, varargin, 'UniformOutput', false);
            b = bigZ(vertcat(cs{:}));
        end

        % ---- arithmetic ----------------------------------------------------------------------
        function r = plus(a, b),  r = bigZ.elementwise(a, b, @(x, y) x.add(y)); end
        function r = minus(a, b), r = bigZ.elementwise(a, b, @(x, y) x.subtract(y)); end
        function r = times(a, b), r = bigZ.elementwise(a, b, @(x, y) x.multiply(y)); end
        function r = uminus(a),   r = bigZ(cellfun(@(x) x.negate(), a.c, 'UniformOutput', false)); end
        function r = uplus(a),    r = a; end

        function r = mtimes(a, b)
            a = bigZ(a);  b = bigZ(b);
            if numel(a.c) == 1 || numel(b.c) == 1, r = times(a, b); return, end
            [m, k] = size(a.c);  [k2, n] = size(b.c);
            if k ~= k2, error('bigZ:dims', 'inner matrix dimensions must agree.'); end
            C = cell(m, n);
            for i = 1:m
                for j = 1:n
                    s = java.math.BigInteger.ZERO;
                    for t = 1:k, s = s.add(a.c{i,t}.multiply(b.c{t,j})); end
                    C{i,j} = s;
                end
            end
            r = bigZ(C);
        end

        function r = power(a, e)
            if isa(e, 'bigZ') || any(e(:) ~= fix(e(:))) || any(e(:) < 0) || ~isscalar(e)
                error('bigZ:power', 'only a nonnegative integer scalar exponent is supported.');
            end
            a = bigZ(a);
            r = bigZ(cellfun(@(x) x.pow(e), a.c, 'UniformOutput', false));
        end
        function r = mpower(a, e)
            if numel(bigZ(a).c) ~= 1, error('bigZ:mpower', 'only a scalar base is supported.'); end
            r = power(a, e);
        end

        function r = conv(a, b)
            a = bigZ(a);  b = bigZ(b);
            na = numel(a.c);  nb = numel(b.c);
            C = repmat({java.math.BigInteger.ZERO}, 1, na + nb - 1);
            for i = 1:na
                for j = 1:nb, C{i+j-1} = C{i+j-1}.add(a.c{i}.multiply(b.c{j})); end
            end
            if iscolumn(a.c) && iscolumn(b.c), C = C.'; end
            r = bigZ(C);
        end

        function r = sum(a, dim)
            if nargin < 2
                if isrow(a.c), dim = 2; else, dim = 1; end
            end
            sz = size(a.c);  sz(dim) = 1;
            C = repmat({java.math.BigInteger.ZERO}, sz);
            for i = 1:numel(a.c)
                sub = cell(1, ndims(a.c));
                [sub{:}] = ind2sub(size(a.c), i);
                sub{dim} = 1;
                k = sub2ind(sz, sub{:});
                C{k} = C{k}.add(a.c{i});
            end
            r = bigZ(C);
        end

        % ---- leaving the type -----------------------------------------------------------------
        function A = double(a)
        % NEAREST doubles -- for display or an evaluator only. ratQ.exactFormula never uses this
        % for its result: it reduces and then requires the reduced value to fit exactly.
            A = cellfun(@(x) x.doubleValue(), a.c);
        end

        function [r, g] = contentFree(a)
        % divide by the POSITIVE gcd of all entries; g is that gcd (BigInteger)
            g = java.math.BigInteger.ZERO;
            for i = 1:numel(a.c), g = g.gcd(a.c{i}); end
            if g.signum() == 0 || g.equals(java.math.BigInteger.ONE), r = a; return, end
            r = bigZ(cellfun(@(x) x.divide(g), a.c, 'UniformOutput', false));
        end

        function A = exactDouble(a, what)
        % the exact doubles, or ratQ:overflow if any entry is past 2^53
            A = zeros(size(a.c));
            for i = 1:numel(a.c)
                v = a.c{i};
                if v.abs().bitLength() > 52
                    error('ratQ:overflow', ...
                        ['%s: reduced in BigInteger, an entry still needs %d bits, past the 53 a ' ...
                         'double holds exactly.'], what, v.abs().bitLength());
                end
                A(i) = v.doubleValue();
            end
        end
    end

    methods (Static, Access = private)
        function r = elementwise(a, b, op)
            a = bigZ(a);  b = bigZ(b);
            if numel(a.c) == 1 && numel(b.c) ~= 1, a = bigZ(repmat(a.c, size(b.c))); end
            if numel(b.c) == 1 && numel(a.c) ~= 1, b = bigZ(repmat(b.c, size(a.c))); end
            if ~isequal(size(a.c), size(b.c))
                error('bigZ:dims', 'array sizes %s and %s do not match.', ...
                      mat2str(size(a.c)), mat2str(size(b.c)));
            end
            r = bigZ(cellfun(op, a.c, b.c, 'UniformOutput', false));
        end
    end
end
