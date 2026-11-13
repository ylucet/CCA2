classdef conjQExamplesSweepTest < matlab.unittest.TestCase
% [fast] conj (default engine, conjQ) over EVERY QuaPol example set, each input classified:
%
%   pass     an answer came back
%   REFUSED  the input is OUT OF SCOPE, by one of these identifiers only -- not a gap, never to be
%            fixed: a cubic numerator (the conjugate is for QuaPol), or an inexact QuaPol built by
%            the legacy float pipeline (the exact engine does not compute on rounded data)
%   FAIL     anything else: an in-scope input the engine does not answer
%
% The two counts are kept apart on purpose (decided 2026-10-11): a refusal is correct behaviour, a
% failure is work. A failure outside `knownFail` fails this test, and so does a refusal under any
% other identifier. Values are checked elsewhere (conjQTest, erratumExamplesTest); this pins the
% boundary.

    properties (Constant)
        outOfScope = {'PLQ:op:unsupportedType', 'PLQ:QuaPol:notExact'};
        % in scope, not answered yet -- each one a TODO item
        knownFail = cell(1,0);   % none (examplesNonconvex{2} answered 2026-10-11)
    end

    methods (Test)
        function everyExampleIsAnsweredOrRefusedAsOutOfScope(t)
            [pass, refused, fail] = conjQExamplesSweepTest.sweep();
            t.log(1, sprintf('conj sweep: %d pass, %d refused (out of scope), %d fail', ...
                             numel(pass), numel(refused), numel(fail)));
            failNames = cellfun(@(c) c{1}, fail, 'UniformOutput', false);
            newFail = setdiff(failNames, conjQExamplesSweepTest.knownFail);
            t.verifyEmpty(newFail, ['in-scope inputs conj does not answer: ' strjoin(newFail, ', ')]);
            fixed = setdiff(conjQExamplesSweepTest.knownFail, failNames);
            t.verifyEmpty(fixed, ['now answered -- remove from knownFail: ' strjoin(fixed, ', ')]);
            t.verifyGreaterThan(numel(pass), 0, 'the sweep must reach the engine');
        end
    end

    methods (Static)
        function [pass, refused, fail] = sweep()
        % Each list holds {name, identifier} pairs (identifier '' for a pass).
            pass = {};  refused = {};  fail = {};
            sets = {'examples', 'examples2', 'examplesNonconvex', 'examplesDiscontinuous'};
            for si = 1:numel(sets)
                P = QuaPol.(sets{si})();
                for k = 1:numel(P)
                    name = sprintf('%s{%d}', sets{si}, k);
                    try
                        P{k}.conj();
                        pass{end+1} = {name, ''}; %#ok<AGROW>
                    catch ME
                        if any(strcmp(ME.identifier, conjQExamplesSweepTest.outOfScope))
                            refused{end+1} = {name, ME.identifier}; %#ok<AGROW>
                        else
                            fail{end+1} = {name, ME.identifier}; %#ok<AGROW>
                        end
                    end
                end
            end
        end
    end
end
