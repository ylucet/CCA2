function sweepCaseEnumeration(seed, nPerCase)
% sweepCaseEnumeration  The Step-1/Step-2 case enumeration, measured over RANDOM instances of each
%   case rather than one hand-picked example per case.
%
% [input]  seed     : rng seed (default 20260802). Seeded so the table is reproducible --
%                     SUPPORT_MATRIX.md section 0.1's rule.
%          nPerCase : how many instances to carry through Step 2 per case (default 12).
%
% WHY THIS EXISTS. The first version of this enumeration used ONE triangle and ONE quadratic per
% case, and that is not enough: the consequences tabulated (face counts, presence of a parabolic
% edge) vary across geometric configurations even when the CLASSIFICATION does not. It was wrong
% in a way that mattered -- for a rank-1 PSD quadratic it reported "5 faces, no parabolic edge",
% but conjPieceCPLQ documents 6 faces with three parabolic edge strips generically and 5 in a
% DEGENERATE tie (a triangle edge exactly perpendicular to A's nonzero eigenvector). The single
% example picked was x^2 on an axis-aligned right triangle, i.e. exactly that degenerate tie, and
% its answer was then reported as the case's answer.
%
% The CLASSIFICATION is exhaustive by construction and needs no sampling: in 2D a symmetric Q is
% PSD, NSD, or indefinite (a trichotomy), and an indefinite one reduces to the bilinear frame
% where the convex-edge count nCE is 0, 1, 2 or 3. What is sampled is only what each class DOES.
%
% Instances are drawn at random and BUCKETED BY MEASURED CLASS rather than constructed to order,
% because nCE is a property of the triangle in the bilinear frame and cannot be prescribed
% directly. That also means the bucket sizes themselves report how generic each case is.

    if nargin < 1 || isempty(seed),     seed = 20260802; end
    if nargin < 2 || isempty(nPerCase), nPerCase = 12;   end
    rng(seed);

    names = {'affine','convex PD','convex rank-1 PSD','concave', ...
             'indefinite 0CE','indefinite 1CE','indefinite 2CE','indefinite 3CE'};
    inst = cell(1, numel(names));
    for i = 1:numel(inst), inst{i} = {}; end

    % ---- draw instances until every bucket has enough (or we run out of tries) ---------------
    tries = 0;
    while tries < 4000 && any(cellfun(@numel, inst) < nPerCase)
        tries = tries + 1;
        V = randTriangle();
        which = mod(tries, 4);
        switch which
            case 0, Q = randPSD();          % convex PD
            case 1, Q = randRank1PSD();     % convex, rank 1
            case 2, Q = -randPSD();         % concave
            case 3, Q = randIndefinite();   % indefinite -> bucketed by nCE
        end
        if mod(tries, 17) == 0, Q = zeros(2); end       % affine, occasionally
        L = randn(2,1); c = randn();
        f6 = [Q(1,1) Q(1,2) Q(2,2) L(1) L(2) c];
        b = bucketOf(Q, V);
        if b == 0, continue, end
        if numel(inst{b}) >= nPerCase, continue, end
        inst{b}{end+1} = {V, f6};
    end

    % ---- run Step 1 and Step 2 on each ------------------------------------------------------
    E = [1 2 1; 2 3 1; 3 1 1]; F = [1 0; 1 0; 1 0];
    fprintf('\nsweepCaseEnumeration(seed=%d, nPerCase=%d)\n', seed, nPerCase);
    fprintf('%-20s %5s | %-26s | %-34s\n', 'CASE', 'n', 'STEP 1 envelope', 'STEP 2 conjugate');
    fprintf('%s\n', repmat('-', 1, 100));
    for b = 1:numel(names)
        s1 = containers.Map('KeyType','char','ValueType','double');
        s2 = containers.Map('KeyType','char','ValueType','double');
        for t = 1:numel(inst{b})
            V = inst{b}{t}{1}; f6 = inst{b}{t}{2};
            q = QuaPol(V, E, f6, F);
            k1 = 'ERROR';
            try
                env = convEnvCPLQ(q);
                kinds = {};
                for k = 1:env.nf
                    if any(abs(env.den(k,1:2)) > sqrt(eps)),        kinds{end+1} = 'quad/lin'; %#ok<AGROW>
                    elseif any(abs(env.f(k,5:7)) > sqrt(eps)),      kinds{end+1} = 'quadratic'; %#ok<AGROW>
                    else,                                          kinds{end+1} = 'affine'; %#ok<AGROW>
                    end
                end
                k1 = sprintf('%df: %s', env.nf, strjoin(unique(kinds),'+'));
            catch ME
                k1 = sprintf('ERR %s', tailId(ME.identifier));
            end
            s1 = bump(s1, k1);

            k2 = 'ERROR';
            try
                g = q.conj('cplq');
                if isa(g,'QuaParCPLQ')
                    k2 = 'symbolic fallback';
                else
                    nPar = 0;
                    if ~isempty(g.Ec), nPar = sum(any(g.Ec ~= 0, 2)); end
                    k2 = sprintf('%df, %d arc', g.nf, nPar);
                end
            catch ME
                k2 = sprintf('ERR %s', tailId(ME.identifier));
            end
            s2 = bump(s2, k2);
        end
        fprintf('%-20s %5d | %-26s | %-34s\n', names{b}, numel(inst{b}), fmtMap(s1), fmtMap(s2));
    end
end

% ============================================================================================
function b = bucketOf(Q, V)
% 1 affine, 2 convex PD, 3 convex rank-1 PSD, 4 concave, 5..8 indefinite with nCE = 0..3.
    tol = 1e-9 * max(1, max(abs(Q(:))));
    ev = sort(eig((Q+Q')/2));
    if max(abs(ev)) <= tol,                         b = 1; return, end
    if ev(1) >= -tol
        if ev(1) > tol, b = 2; else, b = 3; end
        return
    end
    if ev(2) <= tol,                                b = 4; return, end
    % indefinite: count convex edges in the bilinear frame, exactly as convEnvCPLQ does
    M = bilinearFrame(Q);
    b = 5 + numConvexEdges((M*V')');
end

function M = bilinearFrame(Q)
% Same construction as convEnvCPLQ's file-local of this name: 1/2 x'Qx = u1*u2.
    [R, Lam] = eig(Q); lam = diag(Lam);
    [lp, ip] = max(lam); [ln, in] = min(lam);
    r1 = R(:,ip); r2 = R(:,in);
    a = sqrt(lp/2)*r1 + sqrt(-ln/2)*r2;
    bb = sqrt(lp/2)*r1 - sqrt(-ln/2)*r2;
    M = [a'; bb'];
end

function n = numConvexEdges(V)
% Edges of positive finite slope -- convEnvCPLQ's classifyConvexEdges, count only.
    tol = sqrt(eps); ed = [1 2; 2 3; 3 1]; n = 0;
    for t = 1:3
        vi = V(ed(t,1),:); vj = V(ed(t,2),:);
        dx = vj(1)-vi(1); dy = vj(2)-vi(2);
        if abs(dx) < tol, continue, end
        if dy/dx > tol, n = n + 1; end
    end
end

function V = randTriangle()
    for k = 1:50
        V = 4*randn(3,2);
        a = (V(2,1)-V(1,1))*(V(3,2)-V(1,2)) - (V(2,2)-V(1,2))*(V(3,1)-V(1,1));
        if abs(a) > 0.35, if a < 0, V = V([1 3 2],:); end, return, end
    end
end

function Q = randPSD()
    A = randn(2); Q = A'*A + 0.35*eye(2);
end

function Q = randRank1PSD()
    th = 2*pi*rand(); u = [cos(th); sin(th)]; Q = (0.5 + 2*rand())*(u*u');
end

function Q = randIndefinite()
    th = 2*pi*rand(); R = [cos(th) -sin(th); sin(th) cos(th)];
    Q = R*diag([0.5+2*rand(), -(0.5+2*rand())])*R';
    Q = (Q+Q')/2;
end

function m = bump(m, k)
    if isKey(m,k), m(k) = m(k)+1; else, m(k) = 1; end
end

function s = fmtMap(m)
    ks = keys(m); vs = values(m); parts = {};
    [~, ord] = sort(cell2mat(vs), 'descend');
    for i = ord, parts{end+1} = sprintf('%s x%d', ks{i}, vs{i}); end %#ok<AGROW>
    s = strjoin(parts, ', ');
end

function s = tailId(id)
    p = strsplit(id, ':'); s = p{end};
end
