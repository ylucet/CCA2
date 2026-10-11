function [tf, info] = bridgeMember(pr, X)
% bridgeMember  Does a point belong to a BRIDGE face? Decided from the GEOMETRY, not from sign
%   conditions on curves.
%
% objective: the membership test the convex envelope's bridge face actually needs. Its boundary
%   includes a SEGMENT -- a bounded piece of a line -- and an H-form cell can only say "this side of
%   this whole curve", which is why two constructions over sign conditions were built and measured
%   wrong (DECISIONS.md 2026-09-17, and bridgeArcCut's header). This asks the defining question
%   instead, and the defining question has a bounded answer built in.
%
% [input]  pr : the cell pair's own data, as `bridgeFaceAlg` returns it in `face.pair`:
%                 .d, .N1, .m1, .N2, .m2   the two faces over one integer denominator
%                 .Dl, .dl, .kp            the tie's data (Delta, delta, kappa), same denominator
%                 .P                       k x 3 half-planes of the pair; EMPTY for the uncut sweep
%          X  : n x 2 double
% [output] tf   : n x 1 logical
%          info : struct with `mu` (n x 1, the sweeping parameter where one was found, NaN where
%                 not) and `base` (n x 2, the tie point the segment through X starts from)
%
% ==================================================================================================
% THE DEFINING QUESTION
%
% A bridge face is swept: a tie point s carries the segment [grad g1(s), grad g2(s)], and the face is
% the union of those segments over the tie points INSIDE the cell pair. So
%
%       x is in the face  <=>  there is a mu in [-1, 0] with  Phi(x, mu) = 0        (the segment
%                              through x exists, at parameter mu)                    reaches x)
%                              AND its base point s(mu) lies in P .                 (the pair owns
%                                                                                    that tie point)
%
% Both halves are bounded conditions -- mu ranges over a closed interval and P is an intersection of
% half-planes IN THE DUAL PLANE, where everything is straight -- which is exactly what the H-form in
% x could not express. Nothing here is a curve in x at all.
%
% EXACT AT A RATIONAL POINT, CACHED IN DOUBLES HERE. `bridgeMemberExact` answers the same question
% over the integers: Phi's roots are isolated by Sturm and each row of P is the SIGN OF A POLYNOMIAL
% AT A ROOT, decided by `ratQ.signAt`. This routine is its double-precision counterpart, for `eval`
% and for plotting -- the same split QuaCon.eval already documents, and `ratQ.evalFace` before it.
%
% WHAT THIS DOES NOT GIVE, said plainly: a BOUNDARY. A predicate answers "is x inside" and nothing
% else, so figures, edge lists and the paper's topology tables still want the arcs built explicitly.
% That is a separate item and it is not on the path to a correct `biconj`.

    if nargin ~= 2
        error('bridgeMember:nargin', 'expects (pr, X).');
    end
    if isstruct(pr) && isfield(pr, 'kind') && strcmp(pr.kind, 'rational')
        % A RATIONAL face (bridgeFaceRat: affine against rank-1, the tie a parabola) needs no
        % fibre at all: the maximiser along the parabola is t* = -B/(2C), it exists exactly where
        % the objective is concave (C/scale < 0, the scaling's sign included), and the point is
        % owned when that maximiser lies in the pair. Same closed form bridgeFaceRat uses for
        % its witness; until 2026-10-10 this routine refused such a face by name.
        if ~isnumeric(X) || size(X,2) ~= 2
            error('bridgeMember:shape', 'X must be n x 2; got %s.', mat2str(size(X)));
        end
        Bv = X*pr.B(1:2).' + pr.B(3);
        Cv = X*pr.C(1:2).' + pr.C(3);
        t = -Bv ./ (2*Cv);
        r = (pr.rq*t.^2 + pr.rl*t + pr.rc) / pr.rd;
        S = (t*pr.u(:).' + r*pr.w(:).') / pr.nu;
        tf = Cv * pr.scale < 0;
        if ~isempty(pr.P)
            tol = sqrt(eps) * max(1, max(abs(S), [], 2));
            tf = tf & all(S*pr.P(:,1:2).' + pr.P(:,3).' >= -tol, 2);
        end
        info = struct('mu', nan(size(X,1),1), 'base', nan(size(X,1),2));
        return
    end
    if ~isstruct(pr) || ~all(isfield(pr, {'d','N1','m1','N2','m2','Dl','dl','kp','P'}))
        error('bridgeMember:pair', ...
            ['pr must carry the cell pair''s own data (d, N1, m1, N2, m2, Dl, dl, kp, P) -- it is ' ...
             'what `bridgeFaceAlg` puts in `face.pair`.']);
    end
    if ~isnumeric(X) || size(X,2) ~= 2
        error('bridgeMember:shape', 'X must be n x 2; got %s.', mat2str(size(X)));
    end

    n = size(X,1);
    tf = false(n,1);  mus = nan(n,1);  base = nan(n,2);
    for i = 1:n
        [Phi, Dt, U1, U2] = bridgePhi(pr.N1, pr.m1, pr.N2, pr.m2, pr.Dl, pr.dl, pr.kp, pr.d, ...
                                      X(i,:).', 1, false);
        p = fliplr(Phi);                                   % `roots` wants highest power first
        nz = find(abs(p) > eps*max(abs(p)), 1);
        if isempty(nz), continue, end
        p = p(nz:end);
        if numel(p) < 2, continue, end
        r = roots(p);
        r = real(r(abs(imag(r)) <= 1e-9*max(1, abs(r))));
        r = r(r >= -1 - 1e-12 & r <= 1e-12);
        dScale = max(abs(Dt));
        for j = 1:numel(r)
            mu = min(max(r(j), -1), 0);
            D  = polyval(fliplr(Dt), mu);
            % A ROOT WHERE det(M_mu) VANISHES IS AN ARTEFACT OF CLEARING IT. Phi was formed by
            % multiplying through by det(M_mu)^2, so every zero of that determinant is a root
            % whatever x is, and the base point read there is meaningless -- measured on a VERTEX
            % cell, where det(M_mu) = (1+mu)^2 det(N1) vanishes at mu = -1 exactly: the base came
            % back as 1e7 and the face claimed points it does not own. The threshold is RELATIVE to
            % the determinant's own coefficients, because an absolute one means a different thing
            % at every scale. The exact path refuses the same root by the sign of det being 0.
            if abs(D) <= 1e-9 * max(1, dScale), continue, end
            s = [polyval(fliplr(U1), mu); polyval(fliplr(U2), mu)] / D;
            if ~isempty(pr.P)
                % THE TOLERANCE IS ON THE PAIR, NOT ON THE FACE, and that is the cheap direction:
                % a point just outside a facet is swept from a tie point a neighbouring pair owns,
                % and that pair's own face gives the same value there -- the envelope is continuous
                % across a cut. Being slightly generous here therefore cannot produce a wrong
                % VALUE, while being strict would leave a hairline crack between two faces.
                slack = pr.P(:,1:2)*s + pr.P(:,3);
                if any(slack < -1e-9 * max(1, max(abs(s)))), continue, end
            end
            tf(i) = true;  mus(i) = mu;  base(i,:) = s.';
            break
        end
    end
    if nargout > 1
        info = struct('mu', mus, 'base', base);
    end
end
