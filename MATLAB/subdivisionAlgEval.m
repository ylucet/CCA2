function [val, who] = subdivisionAlgEval(subdivision, X)
% subdivisionAlgEval  Evaluate a subdivision whose edges may be ALGEBRAIC -- the counterpart of `QuaCon.eval` for
%   what `envelopeSubdivisionAlg` builds, until the output type can store such an edge.
%
% [input]  subdivision : as `envelopeSubdivisionAlg` returns it -- including the tie rows `.ii` of a minimum over
%                 several pairs, value_1 >= value_2
%          X    : n x 2
% [output] val  : n x 1, +Inf where no face owns the point -- which IS the envelope there, a hole in
%                 the subdivision being how +infinity is expressed
%          who  : n x 1 face index, 0 where none owns it
%
% EVALUATION IS IN DOUBLES, deliberately, exactly as `QuaCon.eval` and `envelopeEval` are: the
% operations that must leave Q are the PREDICATES that decide the subdivision, and reading a finished
% function at a point is not one of them. The algebraic rows are therefore read through `algBranch`
% here; `algCurveSide` is what decides them when the subdivision is being BUILT.
    if nargin ~= 2
        error('subdivisionAlgEval:nargin', 'expects (subdivision, X).');
    end
    if ~isnumeric(X) || size(X,2) ~= 2
        error('subdivisionAlgEval:points', 'X must be n x 2; got %s.', mat2str(size(X)));
    end
    n = size(X,1);
    val = inf(n,1);  who = zeros(n,1);
    for i = 1:n
        x = X(i,:);
        for k = 1:numel(subdivision)
            f = subdivision(k);
            if ~conHolds(f.con, x), continue, end
            if ~algHolds(f.alg, x), continue, end
            if isfield(f, 'ii') && ~iiHolds(f.ii, x), continue, end
            if isfield(f, 'tri') && ~triHolds(f.tri, f.tfs, x), continue, end
            v = faceValue(f, x);
            if isnan(v), continue, end
            if isinf(val(i)) || v > val(i), val(i) = v;  who(i) = k; end
        end
    end
end

function tf = conHolds(con, x)
    tf = true;
    for t = 1:size(con,1)
        c = con(t,1:6);
        if all(c(1:5) == 0)
            if con(t,7)*c(6) < -1e-9, tf = false; return, end
            continue
        end
        br = 0;
        if size(con,2) >= 8, br = con(t,8); end
        if br ~= 0
            v = ratQ.pairLineEval(c, br, x);
        else
            v = c(1)*x(1)^2 + c(2)*x(1)*x(2) + c(3)*x(2)^2 + c(4)*x(1) + c(5)*x(2) + c(6);
        end
        if con(t,7)*v < -1e-9, tf = false; return, end
    end
end

function tf = algHolds(alg, x)
% side * (q(x) - value_FB(x)) >= 0, read in doubles.
    tf = true;
    for t = 1:numel(alg)
        a = alg(t);
        v = algBranch(a.FP, x, a.FB);
        if ~isfinite(v), tf = false; return, end
        q = ratQ.evalFace(a.num, a.den, x);
        if a.side * (q - v) < -1e-9, tf = false; return, end
    end
end

function tf = iiHolds(ii, x)
% value_1(x) >= value_2(x) for two irrational faces (a minimum over pairs), read in doubles.
    tf = true;
    for t = 1:numel(ii)
        v1 = algBranch(ii(t).FP1, x, ii(t).FB1);  v2 = algBranch(ii(t).FP2, x, ii(t).FB2);
        if ~isfinite(v1) || ~isfinite(v2) || v1 - v2 < -1e-9, tf = false; return, end
    end
end

function tf = triHolds(tri, tfs, x)
% side * s_e(x) >= 0 for each row [t e side] -- the e-th side of triangle t, read in doubles from
% its corner cache and oriented +1 toward the opposite corner (triangleFaceSide decides it exactly).
    tf = true;
    for r = 1:size(tri,1)
        T = tfs(tri(r,1)).T;  i = tri(r,2);  j = mod(i,3) + 1;  o = 6 - i - j;
        e = T(j,:) - T(i,:);
        cr = @(p) e(1)*(p(2) - T(i,2)) - e(2)*(p(1) - T(i,1));
        v = cr(x) * sign(cr(T(o,:)));
        if tri(r,3) * v < -1e-9 * max(1, norm(e)) * max(1, norm(x)), tf = false; return, end
    end
end

function v = faceValue(f, x)
    if strcmp(f.kind, 'triangle')
        v = x * f.tf.s - f.tf.c;
    elseif strcmp(f.kind, 'bridge')
        v = algBranch(f.FP, x, f.FB);
    else
        v = ratQ.evalFace(f.num, f.den, x);
    end
end
