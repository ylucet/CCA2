classdef guardExerciseTest < matlab.unittest.TestCase
% EVERY NAMED REFUSAL, ACCOUNTED FOR: one constructed input per guard, asserting that the guard
% FIRES and that it REPORTS ITS OWN IDENTIFIER -- plus a pinned list of the ones not yet exercised,
% each with the reason it is not.
%
% WHY THIS EXISTS, and why it is a test rather than a lint. `QuaPar`'s size guard raised
% `MATLAB:error:nonScalarInput` instead of `QuaPar:size`, because `error` does not expand an array
% across format specifiers the way `sprintf` does -- and the real refusal was masked so completely
% that the actual defect could not be read from the message. A guard that cannot report itself is
% worse than no guard. The obvious detector is static, and a static one was built and REFUTED
% (DECISIONS 2026-09-06): 94 flagged out of 302 named `error()` calls, almost all the parser's own
% bug, and MATLAB's `checkcode` has no printf-argument check for `error` either. That entry records
% what a real detector needs -- "not a parser but an EXERCISE" -- and this is it.
%
% WHAT A PASSING TEST HERE PROVES, AND WHAT IT DOES NOT.
%
%   * `verifyError(..., 'pkg:id')` fails if ANY other error comes out, so a guard that cannot
%     report itself is a RED here. That is the defect class this file is for, and it holds whatever
%     the input's provenance.
%   * REACHABILITY is a different claim. A guard fired by an input built to violate the routine's
%     contract proves self-reporting only -- DECISIONS 2026-09-06 records that argument for
%     `maxQuaPar:notDegenerate`. Every exercise below is marked `reachable` (a legal input a user
%     could write) or `contract` (an input built to break an invariant).
%
% THE INVENTORY IS PINNED BOTH WAYS (`everyNamedGuardIsAccountedFor`): a guard with no exercise, no
% other test asserting it and no remainder entry is a red, and so is a remainder entry naming a
% guard that no longer exists. A stale skip list is how a quarantine becomes a permanently passing
% test.
%
% BUCKET: fast. Every input here is closed-form. The guards whose input is a symbolic `region` or
% `plq_1p` are on the remainder list, because constructing one is the symbolic pipeline.

    methods (Test)

        % ---- argument and type validation, all reachable from a legal call ------------------------
        function quaParNargin(t)                                             % reachable
            t.verifyError(@() QuaPar(1, 2, 3), 'QuaPar:nargin');
        end

        function ratPolNargin(t)                                             % reachable
            t.verifyError(@() RatPol(1, 2, 3), 'RatPol:nargin');
        end

        function toQuaParUnsupportedType(t)                                  % reachable
            t.verifyError(@() toQuaPar(5), 'toQuaPar:unsupportedType');
        end

        function conjPieceCPLQType(t)                                        % reachable
            t.verifyError(@() conjPieceCPLQ(5), 'conjPieceCPLQ:type');
        end

        function cornerFaceAlgArguments(t)                                   % reachable
            % The corner candidates complete a curved pair's coverage, so a caller passing the
            % wrong shape would silently lose the DOMAIN rather than a value.
            n = [0 0 0 0 2 0 2 0 0 0];
            t.verifyError(@() cornerFaceAlg(n, 1, [1 0 0], [1 0 1 0 0 -1]), ...
                          'cornerFaceAlg:nargin');
            t.verifyError(@() cornerFaceAlg([1 2 3], 1, [1 0 0], [1 0 1 0 0 -1], 1), ...
                          'cornerFaceAlg:shape');
            t.verifyError(@() cornerFaceAlg(n, 1, [1 0 0], [1 0 1 0 0 -1], 0), ...
                          'cornerFaceAlg:side');
        end

        function envelopeMeshAlgebraicFace(t)                                % reachable
            % A cell covered ONLY by corner faces would need a face whose value is affine with
            % quadratic-irrational coefficients, which a conic subdivision cannot store. Measured
            % never to happen -- a corner's value always ties a quadratic face (DECISIONS.md
            % 2026-09-19) -- so it is built here on purpose: one pair, one corner, nothing else.
            n = [0 0 0 0, 2, 0, 2, 0, 0, 0];
            f = cornerFaceAlg(n, 1, [1 0 1; -1 0 1; 0 1 1; 0 -1 1], [2 0 2 0 0 -3], 1);
            t.verifyNotEmpty(f, 'the fixture must produce a corner, or this proves nothing');
            c = f(1);  c.pairIndex = 1;
            % The corner must share its PAIR with a quadratic face whose region does not reach as
            % far as the wedge: a corner in a pair of its own is dropped by the coverage rule before
            % the guard, and an envelope with no quadratic face at all returns empty before it.
            q = struct('num', [0 0 0 0, 1, 0, 1, 0, 0, 0], 'den', 1, ...
                       'con', [0 0 0 1 0 -1, 1, 0], 'pair', 1);          % only where x >= 1
            env = struct('quad', q, 'bridge', {{c}});
            t.verifyError(@() envelopeMesh(env), 'envelopeMesh:algebraicFace');
        end

        function cellInteriorArguments(t)                                    % reachable
            % The arrangement's emptiness test. A wrong shape here would not throw on its own -- a
            % six-column table simply reads the branch column as missing -- so the shape is pinned,
            % and a CONSTANT row is refused because admitting one made the exact kernel meet the
            % "curve" 1 = 0 and silently stop the pruning (DECISIONS.md 2026-09-18).
            t.verifyError(@() cellInterior([0 0 0 1 0 0, 0]), 'cellInterior:nargin');
            t.verifyError(@() cellInterior([0 0 0 1 0 0], [1 1]), 'cellInterior:shape');
            t.verifyError(@() cellInterior([0 0 0 0 0 1, 0], [1 1]), 'cellInterior:constant');
        end

        function cornerWedgeRowsArguments(t)                                 % reachable
            % The rows a corner's wedge is spelled with; a caller passing something that is not a
            % corner face would otherwise index fields that are not there.
            n = [0 0 0 0 2 0 2 0 0 0];
            f = cornerFaceAlg(n, 1, [1 0 1; -1 0 1; 0 1 1; 0 -1 1], [2 0 2 0 0 -3], 1);
            t.verifyNotEmpty(f, 'the fixture must produce a corner, or this proves nothing');
            t.verifyError(@() cornerWedgeRows(f(1), 1, 1), 'cornerWedgeRows:nargin');
            t.verifyError(@() cornerWedgeRows(struct('a', 1), n, 1, [2 0 2 0 0 -3], 1), ...
                          'cornerWedgeRows:shape');
        end

        function conicMeetConstant(t)                                        % reachable
            % A CONSTANT is not a curve, and a face's region can carry one. Before the guard the
            % line path formed -f/d with d = e = 0 and threw `ratQ:notFinite` from three calls down.
            t.verifyError(@() conicMeet([1 0 0 0 0 -1], [0 0 0 0 0 1]), 'conicMeet:constant');
        end

        function qsFieldArguments(t)                                         % reachable
            % The Q(t) arithmetic both `cornerWedgeRows` and `cellInterior` rest on. A zero
            % denominator or a t term against a field that has no t are the two ways a caller can
            % ask for something that is not an element.
            t.verifyError(@() ratQ.qsRed([1 2]), 'ratQ:qsShape');
            t.verifyError(@() ratQ.qsRed([1 2 0]), 'ratQ:qsZeroDen');
            t.verifyError(@() ratQ.qsMul([1 1 1], [1 0 1], [1 0]), 'ratQ:qsField');
            t.verifyError(@() ratQ.qsSign([1 1 1], [1 0], [-1 1]), 'ratQ:qsField');
        end

        function bridgeMemberArguments(t)                                    % reachable
            % The membership predicate is called from `eval` and from every test that asks "is x in
            % this face", so its arguments are checked where they arrive.
            pr = struct('d',1,'N1',eye(2),'m1',[0;0],'N2',2*eye(2),'m2',[0;0], ...
                        'Dl',-eye(2),'dl',[0;0],'kp',0,'P',zeros(0,3));
            t.verifyError(@() bridgeMember(pr), 'bridgeMember:nargin');
            t.verifyError(@() bridgeMember(struct('d',1), [0 0]), 'bridgeMember:pair');
            t.verifyError(@() bridgeMember(pr, [0 0 0]), 'bridgeMember:shape');
        end

        function assembleQuaConCellsConShape(t)                              % reachable
            % `con` is 7 columns (conic + side) or 8 (plus the BRANCH that names one line of a
            % degenerate conic). Anything else is a caller that has not been told about the column.
            t.verifyError(@() assembleQuaConCells(struct('num', zeros(1,10), 'den', 1, ...
                                                          'con', [0 0 0 1 0 0])), ...
                          'assembleQuaConCells:conShape');
        end

        function conjConvexPolygonNargin(t)                                  % reachable
            t.verifyError(@() conjConvexPolygon([0 0; 1 0; 1 1], 0, 0, eye(2)), ...
                          'conjConvexPolygon:nargin');
        end

        function conjConvexPolygonShape(t)                                   % reachable
            t.verifyError(@() conjConvexPolygon(zeros(0,2), 0, 0, eye(2), [0;0], 0), ...
                          'conjConvexPolygon:shape');
        end

        function xyFrameNotIndefinite(t)                                     % reachable
            % xyFrame's whole contract is the INDEFINITE case; a PSD Q is the ordinary convex one,
            % handled elsewhere, so refusing it is the guard doing its job.
            t.verifyError(@() xyFrame(eye(2), [0;0], 0), 'xyFrame:notIndefinite');
        end

        function sweepMaxQuaParCurvedSplitSeed(t)                            % reachable
            % A sweep with no seed is unreproducible, and DESIGN 0.1 forbids quoting numbers that
            % cannot be re-derived -- so the sweep refuses to run rather than produce them.
            t.verifyError(@() sweepMaxQuaParCurvedSplit(), 'sweepMaxQuaParCurvedSplit:seed');
        end

        function sweepQuaParEvalAtVerticesSeed(t)                            % reachable
            t.verifyError(@() sweepQuaParEvalAtVertices(), 'sweepQuaParEvalAtVertices:seed');
        end

        function ratQBadQuadSignInput(t)                                     % reachable
            % A conic is 6 wide; three coefficients is not one.
            t.verifyError(@() ratQ.quadSignOnPoly([1 2 3], [1 0 0]), 'ratQ:badQuadSignInput');
        end

        % ---- shape and geometry guards -----------------------------------------------------------
        function quaPolToPlqCurvedEdge(t)                                    % reachable
            % A QuaPar with a NONZERO Ec row is a legal QuaPar -- a curved edge -- and
            % quaPolToPlq's translation covers straight edges only. Boundedness was never the right
            % test here; curvature is, per the guard's own note.
            t.verifyError(@() quaPolToPlq(curvedQuaPar()), 'quaPolToPlq:curvedEdge');
        end

        function maxQuaParNotFullDomain(t)                                   % reachable
            % `F` names the face on each side of each edge, and outside a bounded domain there is
            % none -- so an ordinary triangle QuaPar carries zeros in `F` and is +infinity off the
            % triangle. maxQuaPar's contract is two everywhere-finite functions, so this refusal is
            % reachable from perfectly ordinary input.
            g = triangleQuaPar();
            t.verifyError(@() maxQuaPar(g, g), 'maxQuaPar:notFullDomain');
        end

        function ratPolToPlqZeroDenominator(t)                               % reachable
            V = [0 0; 1 0; 0 1];  E = [1 2 1; 2 3 1; 3 1 1];  F = repmat([1 0], 3, 1);
            r = RatPol(V, E, [0 0 0 0, 2, 0, 2, 0, 0, 0], F, [0 0 0]);
            t.verifyError(@() ratPolToPlq(r), 'ratPolToPlq:zeroDenominator');
        end

        function addQuaParNoOverlap(t)                                       % reachable
            % Two functions on disjoint triangles: the sum's domain is the intersection, which is
            % empty, so there is nothing to return and saying so is the only right answer.
            f = triangleQuaPar();  g = triangleQuaPar(10);
            t.verifyError(@() addQuaPar(f, g), 'addQuaPar:noOverlap');
        end

        function clipArcByConicDegenerateArc(t)                              % contract
            % An arc whose two endpoints coincide has no parameter interval to clip.
            t.verifyError(@() clipArcByConic([1 0 0 0 -1 0], [0 0], [0 0], [0 0 0 1 0 -1]), ...
                          'clipArcByConic:degenerateArc');
        end

        % ---- the two `:cubic` guards are SHADOWED, and that is a measurement, not a gap ----------
        function cubicIsRefusedBeforeThePerFaceGuard(t)
        % `quaPolToPlq:cubic` and `ratPolToPlq:cubic` CANNOT FIRE. Both routines call
        % `assertOperable` at the top (quaPolToPlq.m:52, ratPolToPlq.m:30), which refuses degree > 2
        % with `PLQ:op:unsupportedType`, and the per-face cubic check sits further down (line 76 and
        % line 38). So a cubic face is always refused by the OUTER guard and the inner one is dead
        % code -- measured here rather than asserted, by sending a cubic in and reading which
        % identifier comes out. Recorded in SUPPORT_MATRIX.md; the two inner guards are on the
        % remainder list as `shadowed`.
            V = [0 0; 1 0; 0 1];  E = [1 2 1; 2 3 1; 3 1 1];  F = repmat([1 0], 3, 1);
            cubic = [1 0 0 0, 2, 0, 2, 0, 0, 0];
            t.verifyError(@() quaPolToPlq(QuaPol(V, E, cubic, F)), 'PLQ:op:unsupportedType');
            t.verifyError(@() ratPolToPlq(RatPol(V, E, cubic, F)), 'PLQ:op:unsupportedType');
        end

    end

    methods (Test)
        function everyNamedGuardIsAccountedFor(t)
        % THE PIN. Every `error('pkg:id', ...)` identifier in the project root is exercised above,
        % asserted by some other test, or listed here with its reason class:
        %
        %   symbolic-input : the input is a symbolic `region` / `plq_1p`, so constructing one is the
        %                    symbolic pipeline -- a slow-bucket exercise, not a fast one
        %   loop-safety    : a bound on a loop whose termination is argued at the site; firing it
        %                    means constructing an input that argument says cannot exist
        %   internal       : an internal invariant, so any input that fires it is out of contract
        %                    and the firing would prove self-reporting only, never reachability
        %   contract       : reachability needs an IN-CONTRACT input the pipeline alone can build
        %                    (DECISIONS 2026-09-06 for maxQuaPar:notDegenerate)
        %   shadowed       : an earlier guard refuses the same condition first -- measured, see
        %                    cubicIsRefusedBeforeThePerFaceGuard
        %   todo           : constructible, not yet constructed
        %
        % The count is not the point; the NAMES are, so the failure message can say which guard
        % arrived without an exercise.
            remainder = { ...
                'addQuaPar:degenerateConic',                'todo'; ...
                'addQuaPar:internal',                       'internal'; ...
                'addQuaPar:notImplemented',                 'todo'; ...
                'addQuaPol:internal',                       'internal'; ...
                'bridgeArcCut:witnessOnCut',                'internal'; ...
                'bridgeFaceAlg:branchOverflow',             'todo'; ...
                'bridgeFaceRat:internal',                   'internal'; ...
                'bridgeFaceRat:noWitness',                  'todo'; ...
                'bridgeFaceAlg:branchUndecided',            'internal'; ...
                'bridgeFaceAlg:emptyFibre',                 'internal'; ...
                'bridgeFaceAlg:noStationary',               'internal'; ...
                'bridgeFaceAlg:noZ',                        'internal'; ...
                'bridgeFaceAlg:valueNotARoot',              'internal'; ...
                'bridgeFaceAlg:witnessOverflow',            'todo'; ...
                'clipArcByConic:tooManyCrossings',          'todo'; ...
                'ratQ:parallelPair',                        'internal'; ...
                'conicMeet:noShear',                        'loop-safety'; ...
                'conjConvexOverPiece:noCell',               'symbolic-input'; ...
                'conjConvexOverPiece:nonAffineFacet',       'symbolic-input'; ...
                'conjConvexOverPiece:notConvex',            'symbolic-input'; ...
                'conjConvexPolygon:internal',               'internal'; ...
                'conjConvexPolygon:unreachableFace',        'internal'; ...
                'conjPieceCPLQ:internal',                   'internal'; ...
                'convEnvCPLQ:degenerateFace',               'todo'; ...
                'convEnvCPLQ:internal',                     'internal'; ...
                'convEnvUnbounded:convexAlongRay',          'symbolic-input'; ...
                'convEnvUnbounded:degenerate',              'symbolic-input'; ...
                'convEnvUnbounded:noVertex',                'symbolic-input'; ...
                'convEnvUnbounded:notHalfStrip',            'symbolic-input'; ...
                'convEnvUnbounded:notPointed',              'symbolic-input'; ...
                'convEnvUnbounded:notTriangle',             'symbolic-input'; ...
                'convEnvUnbounded:notWedge',                'symbolic-input'; ...
                'envelopeEval:faceKind',                    'internal'; ...
                'exactQ:internal',                          'internal'; ...
                'exactQ:zeroDenominator',                   'todo'; ...
                'fanUnboundedFace:baseVertex',              'symbolic-input'; ...
                'fanUnboundedFace:noVertex',                'symbolic-input'; ...
                'fanUnboundedFace:nonAffineFacet',          'symbolic-input'; ...
                'maxQuaPar:assemblyDisagreesWithPieces',    'internal'; ...
                'maxQuaPar:notDegenerate',                  'contract'; ...
                'maxQuaPar:pieceInvariant',                 'internal'; ...
                'plqCheck:halfStripRegion',                 'symbolic-input'; ...
                'plqCheck:triRegion',                       'symbolic-input'; ...
                'plqCheck:wedgeRegion',                     'symbolic-input'; ...
                'plq_1p:pieceVars',                         'symbolic-input'; ...
                'quaPolToPlq:cubic',                        'shadowed'; ...
                'ratPolToPlq:cubic',                        'shadowed'; ...
                'ratQ:divideNoConverge',                    'loop-safety'; ...
                'ratQ:gcdNoConverge',                       'loop-safety'; ...
                'ratQ:isolateNoConverge',                   'loop-safety'; ...
                'ratQ:premNoConverge',                      'loop-safety'; ...
                'ratQ:pseudoRemNotIntegral',                'internal'; ...
                'ratQ:signAtNoConverge',                    'loop-safety' };

            ids      = guardIdentifiers();
            asserted = union(assertedIdentifiers(), exercisedHere());

            unhandled = setdiff(ids, union(asserted, remainder(:,1).'));
            t.verifyEmpty(unhandled, sprintf(['these named guards have no exercise, no other ' ...
                'test asserting them, and no remainder entry: %s'], strjoin(unhandled, ', ')));

            stale = setdiff(remainder(:,1).', ids);
            t.verifyEmpty(stale, sprintf(['these remainder entries name no guard that exists ' ...
                '(renamed, deleted, or a typo): %s'], strjoin(stale, ', ')));
        end
    end
end

% ================================================================================================
% ---- the inventory, read from the sources rather than listed ------------------------------------
function ids = guardIdentifiers()
% Every `error('pkg:id'` identifier in the project root. Read from the files, because a hard-coded
% inventory is the thing that goes stale.
    ids = scanFor('error\(''([A-Za-z_]\w*:\w+)''', sourceFiles());
end

function ids = assertedIdentifiers()
% Every `pkg:id` string appearing in some OTHER test file. A proxy -- the string could sit in a
% comment -- and deliberately the generous direction: this pin exists to catch a guard with NO
% coverage at all, so a false "covered" is one missed catch while a false "uncovered" would be
% noise, and noise is what gets a pin disabled.
    ids = scanFor('''([A-Za-z_]\w*:\w+)''', testFiles());
end

function ids = exercisedHere()
% What this file itself asserts, read the same way, so adding an exercise above is all it takes.
    % PER STATEMENT, and the two false starts are worth the three lines they cost. Stopping the
    % match at the first `)` matched nothing, because the call's own function handle carries
    % parentheses; stopping it at the first `;` then missed exactly the three calls whose input
    % contains a column vector `[0;0]`. So the statement is reassembled from its continuation lines
    % and the LAST identifier in it is the one asserted.
    ids = {};
    lines = strsplit(fileread([mfilename('fullpath') '.m']), newline);
    stmt = '';
    for i = 1:numel(lines)
        stmt = [stmt, strtrim(lines{i})]; %#ok<AGROW>
        if endsWith(strtrim(lines{i}), '...'), stmt = extractBefore(stmt, strlength(stmt)-2); continue, end
        if contains(stmt, 'verifyError')
            m = regexp(stmt, '''([A-Za-z_]\w*:\w+)''', 'tokens');
            if ~isempty(m), ids{end+1} = m{end}{1}; end %#ok<AGROW>
        end
        stmt = '';
    end
    ids = unique(ids);
end

function ids = scanFor(pat, files)
    ids = {};
    for i = 1:numel(files)
        txt = fileread(files{i});
        m = regexp(txt, pat, 'tokens');
        for k = 1:numel(m), ids{end+1} = m{k}{1}; end %#ok<AGROW>
    end
    ids = unique(ids);
end

function f = sourceFiles()
    f = listM();
    f = f(~contains(lower(f), 'test'));
end

function f = testFiles()
    f = listM();
    f = f(contains(lower(f), 'test') & ~contains(f, 'guardExerciseTest.m'));
end

function f = listM()
    d = dir(fullfile(fileparts(mfilename('fullpath')), '*.m'));
    f = fullfile({d.folder}, {d.name});
end

% ---- fixtures -----------------------------------------------------------------------------------
function q = curvedQuaPar()
% A one-face QuaPar on a triangle whose first edge carries a nonzero conic row -- a curved edge.
% The shape .claude/guardWiredProbe.m uses, and a legal QuaPar: a parabolic edge satisfies
% b^2 - 4ac = 0 and assertParabolicEdges accepts it.
    V  = [0 0; 1 0; 0 1];
    E  = [1 2 1; 2 3 1; 3 1 1];
    Ec = zeros(3,6);   Ec(1,:) = [1 0 0 0 -1 0];          % y = x^2, a parabola
    q  = QuaPar(V, E, Ec, [0 0 0 0, 2, 0, 2, 0, 0, 0], repmat([1 0], 3, 1));
end

function q = triangleQuaPar(shift)
% |x|^2 on a triangle, optionally translated -- an ordinary bounded-domain QuaPar.
    if nargin < 1, shift = 0; end
    V  = [0 0; 1 0; 0 1] + shift;
    E  = [1 2 1; 2 3 1; 3 1 1];
    q  = QuaPar(V, E, zeros(3,6), [0 0 0 0, 2, 0, 2, 0, 0, 0], repmat([1 0], 3, 1));
end
