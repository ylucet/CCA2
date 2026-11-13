function h = maxQ(g1, g2)
% maxQ  The pointwise maximum of two QuaCons, computed EXACTLY.
%
% objective: h(s) = max(g1(s), g2(s)), with every coefficient and every cell boundary exact.
%   This is Step 3 of [COAP]/[JOGO] -- the fold that combines the per-piece conjugates -- and it is
%   the reason a conjugate needs conic edges at all.
%
% [input]  g1, g2 : QuaCon
% [output] h      : QuaCon
%
% ------------------------------------------------------------------------------------------------
% WHY THE MAX IS WHERE THE CONICS COME FROM
%
% On the overlay cell where g1's face i and g2's face j both apply, the winner changes across
% {g_i = g_j}, whose quadratic part is (H_i - H_j)/2 and whose discriminant is -det(H_i - H_j).
% [COAP]/[JOGO] Theorem 6 makes that vanish when the two pieces compared are ADJACENT pieces of f
% (CONJ_FIELD_PROOF.md Theorem 3), but Step 3 compares ALL pairs, so from three pieces up some
% compared pair is non-adjacent and the boundary is a genuine ELLIPSE. That is precisely why QuaCon
% exists, and this routine is the thing that produces one.
%
% THE ALGORITHM IS SHORT BECAUSE THE H-FORM MAKES IT SHORT. An overlay cell is the CONCATENATION of
% two constraint lists -- there is no geometric intersection to compute, no vertex merging, no
% arrangement. Splitting by the winner is one more constraint on each half. Compare that with the
% V-form: maxQuaPar.m is 4654 lines, and DECISIONS.md's record of it is a sequence of defects about
% clipping, orphaned arcs, collapsed cells and lost far-field pieces -- all of them consequences of
% deciding geometry from stored coordinates.
%
% THE SIGNED DIFFERENCE, NOT ratQ.diffConic. That routine returns the CURVE {g1 = g2} in canonical
% form, which deliberately discards the scale and hence the sign -- exactly what a curve does not
% need and what a comparison cannot do without. Here the difference is built directly so that
% "g1 >= g2" is a statement about a specific integer row, and only then canonicalised, with the
% side carried through the normalisation.
%
% WHAT IS INCOMPLETE. assembleQuaConCells' two stated gaps apply: a cell that is empty only because
% of a CURVED constraint is not detected (so nf is an upper bound), and corners involving a curved
% edge are not named. Both need Phase 2c's exact degree-4 sign kernel. Neither can produce a wrong
% VALUE -- eval reads the exact face functions and the exact sign conditions.

    if ~isa(g1, 'QuaCon') || ~isa(g2, 'QuaCon')
        error('PLQ:maxQ:input', 'maxQ takes two QuaCons; got a %s and a %s.', ...
            class(g1), class(g2));
    end

    cells = struct('num', {}, 'den', {}, 'con', {});
    for i = 1:g1.nf
        ci = rowsOf(g1, i);
        for j = 1:g2.nf
            cj = rowsOf(g2, j);
            base = [ci; cj];                       % the overlay cell: both cells' conditions

            n1 = g1.fN(i,:);  d1 = g1.fD(i);
            n2 = g2.fN(j,:);  d2 = g2.fD(j);

            [dn, dd] = ratQ.sub(n1, d1, n2, d2);   % g1 - g2, exactly
            if all(dn == 0)
                % The SAME function on this overlay cell. Splitting it would produce two cells
                % carrying one function separated by a boundary that is not there.
                cells(end+1) = struct('num', n1, 'den', d1, 'con', base); %#ok<AGROW>
                continue
            end

            % R * m(s) = 2*dd*(g1 - g2)(s), and dd > 0, so sign(R*m) IS sign(g1 - g2).
            R = ratQ.chk([dn(5), 2*dn(6), dn(7), 2*dn(8), 2*dn(9), 2*dn(10)], 'difference');
            Rc = ratQ.conic(R);
            % ratQ.conic may have negated the row; the side has to follow, or each half would
            % silently take the OTHER function -- the same hazard sgnOf guards in conjQ.
            flip = 1;
            nz = find(R ~= 0, 1);
            if R(nz) < 0, flip = -1; end

            % g1 wins where R >= 0, g2 where R <= 0. The shared boundary belongs to both, which is
            % correct rather than sloppy: the two functions are EQUAL there.
            cells(end+1) = struct('num', n1, 'den', d1, ...
                                  'con', [base; Rc,  flip]); %#ok<AGROW>
            cells(end+1) = struct('num', n2, 'den', d2, ...
                                  'con', [base; Rc, -flip]); %#ok<AGROW>
        end
    end

    h = assembleQuaConCells(cells);
end

function rows = rowsOf(g, k)
% objective: face k's constraints as explicit [a b c d e f sign] rows, so that two QuaCons with
%            different edge numberings can be overlaid without translating indices.
    fc = g.FC{k};
    rows = zeros(size(fc,1), 7);
    for r = 1:size(fc,1)
        rows(r,:) = [g.EcQ(fc(r,1), :), fc(r,2)];
    end
end
