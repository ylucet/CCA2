classdef plqStage
% plqStage  Run one stage of the cPLQ pipeline, reusing the previous stage's saved result.
%
% WHY. `testcPLQ/testRectBiconj` was ONE test that ran triangulate -> maximum -> biconjugateF and
% took 3198 s. When it started failing there was nothing to read: no assertion, no stage, just an
% exception at the end of 53 minutes. Splitting it into one test per STAGE, each verifying its own
% output and caching it for the next, turns that into four tests that say WHICH stage broke -- and
% the ones before the break still run in seconds on the next attempt.
%
% HOW. Each stage is keyed by (fixture, stage). `get` loads the cache when it is present and
% newer than every source file the pipeline is built from, and otherwise computes from the
% previous stage and saves. So:
%
%   * a full cold run costs exactly what it costs today, once;
%   * a re-run after an EDIT recomputes -- staleness is decided by file mtime against the
%     repository's own .m files, never by a hand-maintained version number;
%   * a re-run with no edit is instant, which is what makes "run the stage before the failing
%     one" a practical debugging move rather than another 53 minutes.
%
% THE CACHE IS NOT A FIXTURE. It is derived data, it lives under .claude/stagecache/ (which is
% not tracked), and deleting it must change nothing but runtime. Nothing in a test may assert
% against a cached value that the same run did not also verify -- otherwise a stale cache becomes
% a passing test, which is worse than no cache at all.

    methods (Static)
        function d = cacheDir()
            d = fullfile(fileparts(mfilename('fullpath')), '.claude', 'stagecache');
            if ~exist(d, 'dir'), mkdir(d); end
        end

        function t = sourceStamp()
        % Newest mtime over the repository's own .m files. Any edit invalidates every stage.
        % Cheap: a few hundred dir entries, no file is read.
            persistent stamp
            if ~isempty(stamp), t = stamp; return, end
            root = fileparts(mfilename('fullpath'));
            f = dir(fullfile(root, '*.m'));
            t = 0;
            for k = 1:numel(f)
                t = max(t, f(k).datenum);
            end
            stamp = t;
        end

        function out = get(fixture, stage, computeFcn)
        % The result of `stage` for `fixture`, from cache when fresh, else computed and saved.
        % computeFcn takes no arguments and returns the stage's output.
            f = fullfile(plqStage.cacheDir(), sprintf('%s_%s.mat', fixture, stage));
            if exist(f, 'file')
                S = load(f);
                if isfield(S, 'stamp') && S.stamp >= plqStage.sourceStamp()
                    out = S.out;
                    return
                end
            end
            out = computeFcn();
            stamp = plqStage.sourceStamp(); %#ok<NASGU>
            save(f, 'out', 'stamp', '-v7.3');
        end

        function clear(fixture)
        % Drop a fixture's cached stages (all of them when called with no argument).
            d = plqStage.cacheDir();
            if nargin < 1
                pat = '*.mat';
            else
                pat = sprintf('%s_*.mat', fixture);
            end
            f = dir(fullfile(d, pat));
            for k = 1:numel(f), delete(fullfile(d, f(k).name)); end
        end
    end
end
