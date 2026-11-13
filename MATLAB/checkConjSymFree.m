function R = checkConjSymFree(verbose)
% checkConjSymFree  MEASURE how often `conj` reaches the symbolic Case C, and why.
%
% objective: "conj is sym-free except as a fallback" is a claim about a RATE, and a rate has to be
%   counted. This runs `conj('cplq')` over a fixture family spanning every shape the dispatch
%   distinguishes -- full-domain quadratics, single triangles of each Hessian class, bounded
%   polygons, and unbounded multi-face domains -- and reports, per fixture, which route ran and
%   what the fallback reason was when the numeric one declined.
%
% [input]  verbose : (optional) print a row per fixture (default true)
% [output] R : struct array with fields name, kind, symbolic, why, secs
%
% WHY IT IS A CHECK AND NOT A TEST. It asserts nothing about values; it is the instrument that says
% where the remaining symbolic calls are, so the work can be aimed. The value assertions live in
% conjCPLQTest, and the definition oracle lives in plqCheck.
%
% READING THE RESULT. `symbolic` is true exactly when the returned object is a QuaParCPLQ, which is
% cPLQ's own symbolic form and the only thing Case C can produce. `why` is the identifier the
% numeric route raised, recorded through conjCPLQ's own CCA2_CONJ_FALLBACK global.

    if nargin < 1, verbose = true; end
    F = fixtures();
    R = struct('name', {}, 'kind', {}, 'symbolic', {}, 'why', {}, 'secs', {});
    global CCA2_CONJ_FALLBACK %#ok<GVMIS>
    for i = 1:numel(F)
        CCA2_CONJ_FALLBACK = {};
        t0 = tic;
        why = ''; kd = '';
        try
            g = F(i).q.conj('cplq');
            kd = g.kind();
        catch ME
            kd = ['ERROR:' ME.identifier];
        end
        secs = toc(t0);
        if ~isempty(CCA2_CONJ_FALLBACK), why = strjoin(unique(CCA2_CONJ_FALLBACK), ','); end
        R(end+1) = struct('name', F(i).name, 'kind', kd, ...
            'symbolic', strcmp(kd, 'QuaParCPLQ'), 'why', why, 'secs', secs); %#ok<AGROW>
        if verbose
            fprintf('%-46s %-14s sym=%d %7.2fs  %s\n', F(i).name, kd, R(end).symbolic, secs, why);
        end
    end
    CCA2_CONJ_FALLBACK = [];
    if verbose
        n = numel(R); ns = sum([R.symbolic]);
        fprintf('\nSYMBOLIC %d of %d fixtures (%.0f%%)\n', ns, n, 100*ns/max(n,1));
        ws = {R([R.symbolic]).why};
        if ~isempty(ws)
            u = unique(ws);
            for k = 1:numel(u)
                fprintf('   %3d  %s\n', sum(strcmp(ws, u{k})), u{k});
            end
        end
    end
end

% ================================================================================================
function F = fixtures()
% objective: one fixture per SHAPE the conj dispatch distinguishes. Kept small and explicit rather
%   than generated: the point is to cover the branches, and a generated family hides which branch
%   a row is exercising.
    F = struct('name', {}, 'q', {});

    % ---- Case A: full-domain quadratics ------------------------------------------------------
    F(end+1) = mk('A: energy (x^2+y^2)/2 on R^2',        QuaPol([1 0 1 0 0 0]));
    F(end+1) = mk('A: general PD quadratic on R^2',      QuaPol([2 1 3 -1 4 5]));

    % ---- Case B: a single bounded triangle, one per Hessian class ------------------------------
    T = [0 0; 1 0; 0 1];
    F(end+1) = mk('B: affine on a triangle',             tri(T, [0 0 0 1 2 3]));
    F(end+1) = mk('B: convex PD on a triangle',          tri(T, [2 0 2 0 0 0]));
    F(end+1) = mk('B: rank-1 PSD on a triangle',         tri(T, [1 0 0 0 0 0]));
    F(end+1) = mk('B: concave on a triangle',            tri(T, [-2 0 -2 0 0 0]));
    F(end+1) = mk('B: indefinite xy on a triangle',      tri(T, [0 1 0 0 0 0]));
    F(end+1) = mk('B: indefinite x^2-y^2 on a triangle', tri(T, [1 0 -1 0 0 0]));
    F(end+1) = mk('B: A.3 one-convex-edge triangle',     tri([1 1; 0 0; 2 0], [0 1 0 0 0 0]));
    F(end+1) = mk('B: A.5 three-convex-edge triangle',   tri([5/2 3/2; 0 0; 1/2 1], [0 1 0 0 0 0]));

    % ---- Case B2: bounded polygons -------------------------------------------------------------
    S = [0 0; 1 0; 1 1; 0 1];
    F(end+1) = mk('B2: convex quadratic on a square',    poly4(S, [2 0 2 0 0 0]));
    F(end+1) = mk('B2: concave on a square',             poly4(S, [-2 0 -2 0 0 0]));
    F(end+1) = mk('B2: xy on a square (two triangles)',  twoTri(S, [0 1 0 0 0 0]));
    F(end+1) = mk('B2: xy on a general quadrilateral',   twoTri([0 0; 2 0; 5/2 1; 1/2 1], [0 1 0 0 0 0]));

    % ---- multi-piece: the case that forces a non-adjacent comparison ----------------------------
    F(end+1) = mk('C?: two DIFFERENT quadratics on a square', twoTriTwoF(S, [2 0 2 0 0 0], [0 1 0 0 0 0]));

    % ---- unbounded ------------------------------------------------------------------------------
    F(end+1) = mk('U: 4 convex quadratic wedges (unbounded)', fourConvexWedges());
    F(end+1) = mk('U: max(0,x,y) as three AFFINE wedges',     maxZeroXY());
end

function q = fourConvexWedges()
% Four unbounded faces round the origin, each carrying a CONVEX quadratic. Taken from
% conjCPLQTest/multiFaceUnboundedConvexFacesConjugateExactly, which is the test that pins the
% values; here it is the ROUTE that is being measured.
    V = [0 0; -1 0; 0 1; 1 0; 0 -1];
    E = [1 2 0; 1 3 0; 1 4 0; 1 5 0];
    f = [1 0 1 0 0 0; 1 0 2 0 0 0; 2 0 2 0 0 0; 2 0 1 0 0 0];
    Fa = [1 2; 2 3; 3 4; 4 1];
    q = QuaPol(V, E, f, Fa);
end

function q = maxZeroXY()
% max(0, x, y) as three unbounded wedges: convex, its own biconjugate, and every face AFFINE.
% Included because an affine unbounded face is a DIFFERENT gap from a convex quadratic one --
% conjConvexPolygon needs a positive definite A, and the conjugate of an affine function over an
% unbounded polygon is a support function, i.e. an INDICATOR, not a quadratic.
    V = [0 0; 0 -1; -1 0; 1 1];
    E = [1 2 0; 1 3 0; 1 4 0];
    Fa = [2 1; 1 3; 3 2];
    f = [0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0];
    q = QuaPol(V, E, f, Fa);
end

function s = mk(name, q), s = struct('name', name, 'q', q); end

function q = tri(V, f6)
    E = [1 2 1; 2 3 1; 3 1 1];
    F = [1 0; 1 0; 1 0];
    q = QuaPol(V, E, f6, F);
end

function q = poly4(V, f6)
    E = [1 2 1; 2 3 1; 3 4 1; 4 1 1];
    F = [1 0; 1 0; 1 0; 1 0];
    q = QuaPol(V, E, f6, F);
end

function q = twoTri(V, f6)
% the same quadratic on both halves -- Step 0 merges these back into one face
    E = [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1];
    F = [1 0; 1 0; 2 1; 2 0; 2 0];
    q = QuaPol(V, E, [f6; f6], F);
end

function q = twoTriTwoF(V, fa, fb)
% two DIFFERENT quadratics: Step 0 cannot merge them, so Step 3 must compare two pieces
    E = [1 2 1; 2 3 1; 1 3 1; 3 4 1; 4 1 1];
    F = [1 0; 1 0; 2 1; 2 0; 2 0];
    q = QuaPol(V, E, [fa; fb], F);
end

