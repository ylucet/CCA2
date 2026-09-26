classdef exactSumLintTest < matlab.unittest.TestCase
% [fast] No `ratQ.chk` around an explicit SUM OF PRODUCTS in production code.
%
% chk tests the value it is given, and a sum of products can be small while a product in it was
% rounded past 2^53 -- the rounded sum then has the wrong sign and passes. That returned WRONG
% answers, silently, from ratQ.feasible2, ratQ.detExact and ratQ.add (DECISIONS 2026-09-25), and the
% same shape was at 80 more sites. A sum of products goes through ratQ.chkB(value, |expression|)
% instead, which refuses once the terms reach 2^53.
%
% WHAT THIS DOES NOT SEE: a matrix or dot product inside chk (`chk(N*s + m)`) hides the same sum
% behind one `*`, and a line scan cannot tell it from a scalar product, which is safe. Those are a
% TODO.md item, not a pass of this test.

    methods (Test)
        function noChkWrapsAnExplicitSumOfProducts(t)
            root = fileparts(mfilename('fullpath'));
            files = dir(fullfile(root, '*.m'));
            bad = {};
            for k = 1:numel(files)
                if endsWith(files(k).name, 'Test.m'), continue, end
                lines = strsplit(fileread(fullfile(root, files(k).name)), newline);
                for i = 1:numel(lines)
                    L = lines{i};
                    c = strfind(L, '%');  if ~isempty(c), L = L(1:c(1)-1); end   % code only
                    if contains(L, 'chkB') || contains(L, 'mulChk'), continue, end
                    for a = chkFirstArgs(L)
                        if ~isempty(regexp(a{1}, '\*.*[-+].*\*', 'once'))
                            bad{end+1} = sprintf('%s:%d', files(k).name, i); %#ok<AGROW>
                            break
                        end
                    end
                end
            end
            t.verifyEmpty(bad, sprintf('chk around a sum of products (use ratQ.chkB): %s', strjoin(bad, ', ')));
        end
    end
end

function args = chkFirstArgs(L)
% The FIRST argument of every `chk(` call on the line, by paren matching, so two calls on one line
% (`[ratQ.chk(-a*Dl,'c8'), ratQ.chk(-b*Dl,'c9')]`) are read as two single products, not one sum.
    args = {};
    for st = strfind(L, 'chk(')
        if st > 1 && any(L(st-1) == 'B'), continue, end
        depth = 0;  j = st + 4;  from = j;
        while j <= numel(L)
            ch = L(j);
            if ch == '(' || ch == '['
                depth = depth + 1;
            elseif ch == ')' || ch == ']'
                if depth == 0, break, end
                depth = depth - 1;
            elseif ch == ',' && depth == 0
                break
            end
            j = j + 1;
        end
        args{end+1} = L(from:j-1); %#ok<AGROW>
    end
end
