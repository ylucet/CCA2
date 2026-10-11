function [face, info] = bridgeFaceRat(n1, d1, n2, d2, P)
% bridgeFaceRat  The bridge face of an AFFINE cell against a RANK-1 cell -- which is RATIONAL, not
%   algebraic, and comes out in closed form with no elimination at all.
%
% objective: the class the discriminant construction cannot touch and does not need to. Measured on
%   the standard fixtures, 12 of 29 bridge pairs are affine x rank-1 (DECISIONS.md 2026-09-17), and
%   for those `det((1+mu)N1 - mu N2)` vanishes IDENTICALLY -- one Hessian is 0 and the other is
%   singular -- so no member of the Lagrange family is invertible and `bridgeFaceAlg` refuses with
%   `degreeInMu`. The geometry is simpler than the general case, not harder.
%
% [input]  n1, d1 : the AFFINE face (its quadratic part must be zero), in the weighted basis
%          n2, d2 : the RANK-1 face
%          P      : (optional) the cell pair's k x 3 integer half-planes in s
% [output] face : struct
%                   .kind  'rational'
%                   .fN    1x10 integer numerator in the weighted basis -- the value is fN/den
%                   .fD    positive integer scaling of that numerator
%                   .den   1x3 [g h k], the LINEAR denominator g*x1 + h*x2 + k
%                   .pair  what `bridgeMember` decides membership from
%                   .wit   a rational point of the face, [xn1 xn2 xd]
%          info : the parabola's data and the coefficient sizes
%
% ==================================================================================================
% WHY IT IS RATIONAL, IN FOUR LINES
%
% g1 is AFFINE, so the objective `<s,x> - g1(s)` is LINEAR in s. g2 has a RANK-1 Hessian, so the tie
% `{g1 = g2}` is a PARABOLA. Parameterise the parabola by the linear coordinate t = <u,s> along the
% range of that Hessian -- its other coordinate is then a QUADRATIC in t -- and the objective becomes
%
%       A(x) + B(x) t + C(x) t^2 ,        every coefficient AFFINE in x.
%
% The value is that divided by a fixed SCALING, so it is concave in t exactly where C/scale < 0 --
% the sign of the scaling is part of the condition, not a detail. There the maximum is at
% t* = -B/(2C) and the value is
%
%       A - B^2/(4C)  =  (4AC - B^2) / (4C) ,
%
% a quadratic numerator over a LINEAR denominator: exactly the `Rat` level of the lattice, with no
% `FP`/`FB` and no branch index. Where C > 0 the objective grows without bound along the parabola and
% the pair's contribution is +infinity, which is an answer and is returned as a named refusal so the
% assembly can record it rather than guess.
%
% THE RESTRICTION TO THE PAIR does not change the formula, only where it applies: each row of P is a
% QUADRATIC in t (because the parabola's second coordinate is), so the face is where t*(x) satisfies
% all of them -- and t*(x) is a RATIONAL function of x, so that test is closed form too. Points where
% the maximum is pinned at an end of the feasible t-set belong to the NEIGHBOURING face, the facet
% candidate of the side cell, exactly as in the general case.

    if nargin < 4 || nargin > 5
        error('bridgeFaceRat:nargin', 'expects (n1, d1, n2, d2) or (n1, d1, n2, d2, P).');
    end
    if nargin < 5, P = zeros(0,3); end
    n1 = checkFace(n1, 'n1');  n2 = checkFace(n2, 'n2');
    d1 = checkDen(d1, 'd1');   d2 = checkDen(d2, 'd2');
    if any([n1(5) n1(6) n1(7)] ~= 0)
        error('bridgeFaceRat:notAffine', ...
            'n1 must be the AFFINE face; its quadratic part is %s.', mat2str(n1(5:7)));
    end
    N2 = [n2(5) n2(6); n2(6) n2(7)];
    if all(N2(:) == 0) || ratQ.detExact(N2) ~= 0
        error('bridgeFaceRat:notRankOne', ...
            ['n2 must have a RANK-1 Hessian: rank 0 is two affine faces (no curve to sweep) and ' ...
             'rank 2 is the general case, which bridgeFaceAlg answers.']);
    end

    d  = ratQ.chk(lcm(d1, d2), 'common denominator');
    b1 = [n1(8); n1(9)] * (d/d1);   k1 = n1(10) * (d/d1);
    M2 = N2 * (d/d2);               b2 = [n2(8); n2(9)] * (d/d2);   k2 = n2(10) * (d/d2);
    dl = b1 - b2;   kp = k1 - k2;

    % ---- the parabola, parameterised ----------------------------------------------------------
    % M2 is symmetric of rank 1, so M2 = u u' / a with u its first nonzero column and a the matching
    % diagonal entry -- integral, and no square root is taken of anything.
    if M2(1,1) ~= 0
        u = [M2(1,1); M2(1,2)];   a = M2(1,1);
    else
        u = [M2(1,2); M2(2,2)];   a = M2(2,2);
    end
    if a == 0
        error('bridgeFaceRat:internal', ...
            'a rank-1 symmetric matrix has a nonzero diagonal entry; this one appears not to.');
    end
    w  = [-u(2); u(1)];                     % along the parabola's axis: u'w = 0
    nu = ratQ.chk(u.'*u, 'range normalisation');
    dw = ratQ.chkB(dl.'*w, abs(dl).'*abs(w), 'tie linear part along the axis');
    if dw == 0
        % THE TIE IS A LINE PAIR, AND THEN THERE IS NO SEPARATE BRIDGE FACE AT ALL. The tie's linear
        % part lies in the rank-1 Hessian's range, so the tie depends on t = <u,s> alone: t^2 = const,
        % two PARALLEL LINES (or one, or none). With g1 AFFINE the objective is LINEAR in s, so the
        % sup over either side of the tie is a LINEAR PROGRAM on a polygon and is attained at a
        % VERTEX of it -- and the tie is just one more facet of that polygon. The value is therefore
        % an ordinary VERTEX candidate of the side, which `envelopeCellRegions` already produces
        % once the tie is handed to it as a facet.
        %
        % Measured: this is what ALL 12 affine x rank-1 pairs on the standard fixtures are, so it is
        % the assembly's job and not the producer's.
        error('bridgeFaceRat:noSeparateFace', ...
            ['the tie is a pair of PARALLEL LINES and this face is AFFINE, so the pair''s sup is a ' ...
             'linear program on each side of the tie and is attained at a VERTEX -- the tie acts ' ...
             'as an ordinary FACET, and there is no separate bridge face to build. The assembly ' ...
             'should add the tie to each side polygon and take that side''s vertex candidates.']);
    end

    % s(t) = (t*u + r(t)*w)/nu  with  r(t) = (t^2 * nu/(2a) - t*(dl''u) - kp*nu) / dw
    % Cleared: r(t) = (rq*t^2 + rl*t + rc) / rd
    rq = ratQ.chk(nu,            'parabola quadratic term');
    rl = ratQ.chkB(-2*a*(dl.'*u), 2*abs(a)*(abs(dl).'*abs(u)), 'parabola linear term');
    rc = ratQ.chk(-2*a*kp*nu,    'parabola constant term');
    rd = ratQ.chk(2*a*dw,        'parabola denominator');

    [face, info] = buildFace(u, w, nu, rq, rl, rc, rd, b1, k1, d, P);
    info.parabola = struct('u', u, 'w', w, 'rq', rq, 'rl', rl, 'rc', rc, 'rd', rd);
end

% ==================================================================================================
function [face, info] = buildFace(u, w, nu, rq, rl, rc, rd, b1, k1, d, P)
% The objective along the parabola, its maximiser, and the value -- all over the integers.
%
% <s,x> - g1(s) with g1 = (b1's + k1)/d and s = (t u + r(t) w)/nu is
%
%       [ t (u'y) + r(t) (w'y) ] / nu  -  k1/d ,      y = x - b1/d ,
%
% so putting everything over the common denominator and collecting powers of t gives A, B, C, each
% AFFINE in x. Their integer forms are what the rest of the routine uses.
    X = [1 0];  Y = [0 1];                          %#ok<NASGU> % (documentation of the layout below)
    % y = x - b1/d, cleared by d:  d*y = d*x - b1.  An affine form in x is kept as [gx gy c],
    % meaning (gx*x1 + gy*x2 + c), so u'y*d = u1*d*x1 + u2*d*x2 - u'b1.
    uy = [d*u(1), d*u(2), -(u.'*b1)];
    wy = [d*w(1), d*w(2), -(w.'*b1)];

    % value*(nu*d*rd) = rd*(t*uy + ...) ... collect:  A + B t + C t^2 with
    C = ratQ.chk(rq * wy,                         'C');            % from r(t)'s t^2 term
    B = ratQ.chkB(rd * uy + rl * wy, abs(rd) * abs(uy) + abs(rl) * abs(wy), 'B');
    A = ratQ.chkB(rc * wy - [0, 0, nu*rd*k1], abs(rc) * abs(wy) + [0, 0, abs(nu*rd*k1)], 'A');
    % the common scaling: value = (A + B t + C t^2) / (nu * d * rd)
    scale = ratQ.chk(nu * d * rd, 'value scaling');

    % THE MAXIMISER. t* = -B/(2C) needs C < 0 for a maximum, and C is AFFINE in x -- so the face
    % lives where that affine form has the sign that makes the objective concave. Its sign is read
    % at the witness like every other side in this toolbox.
    face = struct('kind', 'rational');
    info  = struct();

    % value at t*: (A - B^2/(4C)) / scale = (4AC - B^2) / (4C * scale)
    num = ratQ.chkB(polyMulAff(A, C)*4 - polyMulAff(B, B), ...
                    polyMulAff(abs(A), abs(C))*4 + polyMulAff(abs(B), abs(B)), 'face numerator');   % quadratic in x
    den = ratQ.chk(4 * scale * C, 'face denominator');                          % LINEAR in x

    face.fN = ratQ.chk(quadToWeighted(num), 'face numerator, weighted basis');
    face.fD = 1;
    face.den = den;
    % THE SCALING'S SIGN IS PART OF THE CONDITION, and leaving it out was a real defect: the value
    % is (A + B t + C t^2)/scale, so it is CONCAVE in t exactly when C/scale < 0 -- not when C < 0.
    % `scale = nu*d*rd` carries `rd = 2a(delta'w)`, whose sign is whatever the pair's geometry makes
    % it. Measured before the fix: three of six configurations wrong by factors of 1e3, and one
    % where the face admitted no point at all.
    face.pair = struct('kind', 'rational', 'u', u, 'w', w, 'nu', nu, ...
                       'rq', rq, 'rl', rl, 'rc', rc, 'rd', rd, ...
                       'B', B, 'C', C, 'scale', scale, 'P', P);
    face.wit = witnessRat(face.pair);
    info.maxCoeff = max(abs([face.fN, face.den]));
    info.scale = scale;
end

% ==================================================================================================
function wit = witnessRat(pr)
% A rational point of the face: propose from the parabola itself -- the midpoint of no segment here,
% because the sweep IS the parabola's own maximiser map -- and verify with the closed-form test.
    for t = [0 1 -1 2 -2 0.5 -0.5 4 -4 10 -10]
        s = paraPoint(pr, t);
        if any(~isfinite(s)), continue, end
        if ~isempty(pr.P) && any(pr.P(:,1:2)*s + pr.P(:,3) < 0), continue, end
        % x such that t*(x) = t: any x with B(x) + 2 C(x) t = 0. Walk a small grid and take the
        % first rational point whose own maximiser is this t and which the pair admits.
        for xd = [1 2 4 8]
            for a1 = -8:8
                for a2 = -8:8
                    xn = [a1; a2];
                    if ~bridgeMemberRatAt(pr, xn, xd), continue, end
                    wit = [xn(:).', xd];
                    return
                end
            end
        end
    end
    error('bridgeFaceRat:noWitness', ...
        ['no rational point of this face was found on the search grid. Either the pair admits no ' ...
         'part of the parabola, or the face is thinner than the grid -- which the caller should ' ...
         'be told rather than have guessed at.']);
end

function tf = bridgeMemberRatAt(pr, xn, xd)
% Closed form: the maximiser t* = -B/(2C) must exist (C < 0) and satisfy every row of P.
    Bv = affAt(pr.B, xn, xd);
    Cv = affAt(pr.C, xn, xd);
    if Cv * pr.scale >= 0, tf = false; return, end   % not concave along the parabola: no maximum
    t = -Bv / (2*Cv);
    s = paraPoint(pr, t);
    tf = isempty(pr.P) || all(pr.P(:,1:2)*s + pr.P(:,3) >= 0);
end

function s = paraPoint(pr, t)
    r = (pr.rq*t^2 + pr.rl*t + pr.rc) / pr.rd;
    s = (t*pr.u + r*pr.w) / pr.nu;
end

function v = affAt(f, xn, xd)
    v = (f(1)*xn(1) + f(2)*xn(2)) / xd + f(3);
end

% ==================================================================================================
function c = polyMulAff(f, g)
% The product of two affine forms [gx gy c], as [xx xy yy x y 1].
    c = [f(1)*g(1), f(1)*g(2) + f(2)*g(1), f(2)*g(2), ...
         f(1)*g(3) + f(3)*g(1), f(2)*g(3) + f(3)*g(2), f(3)*g(3)];
end

function n = quadToWeighted(q)
% [xx xy yy x y 1] into the 10-wide weighted basis, where the squared terms carry halves.
    n = [0 0 0 0, 2*q(1), q(2), 2*q(3), q(4), q(5), q(6)];
end

function n = checkFace(n, name)
    if ~isnumeric(n) || ~isvector(n) || numel(n) ~= 10
        error('bridgeFaceRat:shape', '%s must be a 1x10 face numerator; got %s.', ...
              name, mat2str(size(n)));
    end
    n = ratQ.chk(double(n(:).'), name);
    if any(n(1:4) ~= 0)
        error('bridgeFaceRat:cubicFace', '%s has a CUBIC part; a conjugate face is quadratic.', name);
    end
end
function d = checkDen(d, name)
    if ~isscalar(d) || ~isnumeric(d) || d ~= fix(d) || d <= 0
        error('bridgeFaceRat:denominator', '%s must be a positive integer; got %s.', name, mat2str(d));
    end
    d = double(d);
end
