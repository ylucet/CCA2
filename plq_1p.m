classdef plq_1p
    properties
        f;
        d;
        envelope = functionNDomain.empty();
        conjfia = [];
        conjugates = functionNDomain.empty();
        maxConjugate = functionNDomain.empty();
        % Set by convexEnvelope when this piece's quadratic is INDEFINITE but not already x*y:
        % struct('M',M,'a',a,'c0',c0) recording the change of variables q(Mz) = z1*z2 + a'z + c0.
        % conjugate() then works in the z-frame and substitutes back via f*(s) = h*(M's-a) - c0.
        % Empty means "no frame change was needed", which is every input cPLQ itself ever had.
        frame = [];
        % For printouts
        lCE = false
        lConj = false
        lMConj = false
        lPrintEta = false %true
    end




    methods % creation & print
         function obj = plq_1p(d,f)
            % put checks for type of f and d
            if nargin > 0
              obj.f = f;
              obj.d = d;
            end 
         end

         function print(obj)
         
           disp("Domain")
           obj.d.print
           fprintf("\n")
           disp("Function")
          
           obj.f.f = simplifyFraction(obj.f.f);
           obj.f.print
           fprintf("\n\n\n")
           
           disp("Convex Envelope")
           for j=1: size(obj.envelope,2)
             if obj.lCE  
                 disp('Function')  
                
                 obj.envelope(j).f.f = simplifyFraction(obj.envelope(j).f.f);
                 obj.envelope(j).f.print
                 disp('Domain')
                 obj.envelope(j).d.print
             end
             if obj.lConj
                 if (size(obj.conjfia,1) > 0)
                     obj.conjugates(obj.conjfia(j):obj.conjfia(j+1)-1).printL
                 end
             end
           end
           fprintf("\n\n\n\n\n")
           if ~obj.lMConj
               return
           end
           disp("Maximum conjugate")
           obj.maxConjugate.printL
         end

         function printLatex(obj)
         
           disp("Domain")
           %obj.d.printLatex
           obj.d.polygon.printLatex
           fprintf("\n")
           disp("Function")
           
           obj.f.printLatex
           fprintf("\n\n\n")
              
           disp("Convex Envelope")
           disp(" ")
           for j=1:size(obj.envelope,2) 
             disp('Function')  
             obj.envelope(j).f.f = simplifyFraction(obj.envelope(j).f.f);
             obj.envelope(j).f.printLatex
             disp('Domain')
             disp(" ")
             obj.envelope(j).d.printLatex
             if (size(obj.conjfia,1) > 0)
                 obj.conjugates(obj.conjfia(j):obj.conjfia(j+1)-1).printLLatex
             end
           end
           fprintf("\n\n\n\n\n")
           disp("Maximum conjugate")
           obj.maxConjugate.printLLatex
         end


         
         function Mprint(obj)
           fprintf("display(inequal(("); 
           obj.d.polygon.printMaple
           
           fprintf("),x=-5..5,y=-5..5,color=[red,blue,yellow,green],nolines)) \n")
           if obj.lCE  
             
           disp("Convex Envelope")
           fprintf("display(inequal((");
           for j=1:size(obj.envelope,2)-1 
             obj.envelope(j).d.printMaple
             fprintf(",");
           end
           obj.envelope(size(obj.envelope,2)).d.printMaple
           fprintf("),x=-5..5,y=-5..5,color=[red,blue,yellow,green],nolines)) \n")
           end
           if obj.lConj
             
           

           for j=1:size(obj.envelope,2) 
           if (size(obj.conjfia,1) > 0)
                  obj.conjugates(obj.conjfia(j):obj.conjfia(j+1)-1).printM
           end
           end
           end
           if ~ obj.lMConj
             return
           end
           disp('max in piece')
           obj.maxConjugate.printM;
         end

         function plotMaxConjugateDomain(obj)
             figure;
             colors = ['b', 'r', 'g', 'm', 'c', 'y'];
             n = 0
             f = obj.maxConjugate (1).f
             c = colors(mod(n,6)+1)
             for i =1:size(obj.maxConjugate,2)
                if (f.f ~= obj.maxConjugate (i).f.f)
                  n = n + 1
                  c = colors(mod(n,6)+1)
                  f = obj.maxConjugate (i).f
                end
                obj.maxConjugate (i).d.plot;
                textR = "R"+num2str(i);
                textR="";
                obj.maxConjugate (i).d.plotRegionC(textR,c);
             end
         end

         function plotDomain(obj)
             figure;
             obj.d.plot;
         end
    end

    methods
        function ps = triangulate (obj, ps)

            d = obj.d;
            vars = d.polygon.vars;

            % An UNBOUNDED face cannot go down the fan below: that fan reads vx/vy as
            % coordinates and rebuilds each triangle through the BOUNDED domain(t,x,y)
            % constructor, and an unbounded region stores its directions as vertices flagged
            % +/-intmax. Feeding those through produced envelopes carrying 2147483647 and
            % intmax^2 = 4611686014132420609 (measured -- see quaPolToPlq's header). Such a face
            % is covered by fanUnboundedFace instead, which works from the HALF-PLANES and emits
            % triangles, half-strips and wedges. The cover is a cover, not a partition, which is
            % all plq.maximum needs: a sup over a union is the max of the sups.
            try
                [~, rkind] = d.polygon.recessionRays;
            catch
                rkind = 'bounded';      % non-affine facet: not a case this branch handles
            end
            if ~strcmp(rkind, 'bounded')
                ds = fanUnboundedFace(d.polygon, vars(1), vars(2));
                for i = 1:numel(ds)
                    ps = [ps, plq_1p(ds{i}, obj.f)]; %#ok<AGROW>
                end
                return
            end

            % A CONVEX (or affine) quadratic needs NO TRIANGULATION -- and triangulating it is
            % expensive and slightly wrong. co(q|P) = q on any convex P, so Step 1 has nothing to
            % do (convexEnvelope1's own convex branch says exactly this), and Step 2 conjugates
            % it by the KKT active set -- vertex, edge and interior cells -- which
            % conjConvexOverPiece computes for a region with ANY number of affine facets, not
            % just three. Splitting the polygon first only forces Step 3 to glue back together
            % what was never broken.
            %
            % MEASURED on f = (x^2+y^2)/2 over the unit square, which is the QPLIB diagonal-term
            % shape:
            %     via triangulate   2 triangles -> 7 + 7 cells -> Step 3 -> 9 cells, 456 s,
            %                       worst error 1.6e-4 against the closed form
            %     direct            9 cells, 4.2 s, EXACT at 10 of 10
            % The 9 cells are the product structure the answer actually has: f separates and the
            % box is a product, so f*(s) = g(s1) + g(s2) with g the 1-D conjugate, three pieces
            % each. Triangulation destroys that structure and Step 3 pays to rebuild it.
            if any(strcmp(obj.quadKind, {'convex', 'affine'})) || obj.separableParts
                ps = [ps, obj];
                return
            end

            if d.polygon.nv == 3
                ps = obj.appendTriangle(ps, [d.polygon.vx(:), d.polygon.vy(:)], vars);
                return
            end
            vx = d.polygon.vx;
            vy = d.polygon.vy;

            % (mvx,mvy) : vertex in up left corner
            mvy = max(vy);
            ind = [];
            for i = 1:d.polygon.nv
                if vy(i) ~= mvy
                    ind = [ind,i];
                end
            end
            vx(ind) = [];
            mvx = min(vx);
            vx = d.polygon.vx;
            for i = 1:d.polygon.nv
                if vx(i) == mvx & vy(i) == mvy
                    break;
                end
            end

            start = i;
            

            for i = start+1:d.polygon.nv-1
                triangle = [start,i,i+1];
                 t(:,1) = vx(triangle);
                 t(:,2) = vy(triangle);
                 ps = obj.appendTriangle(ps, t, vars);


            end
            if start ~= 1 & start ~= d.polygon.nv
                triangle = [start,d.polygon.nv,1];
                 t(:,1) = vx(triangle);
                 t(:,2) = vy(triangle);
                 ps = obj.appendTriangle(ps, t, vars);
            end
            for i = 1: start-2
                triangle = [start,i,i+1];
                 t(:,1) = vx(triangle);
                 t(:,2) = vy(triangle);
                 ps = obj.appendTriangle(ps, t, vars);
            end

        end

        function ps = appendTriangle(obj, ps, t, vars)
        % Append the triangle `t` (3x2 of vertex coordinates) as one piece -- or, when cPLQ's own
        % closed form is not the convex ENVELOPE there, as the sub-triangles on which it is.
        %
        % [COAP] Appendix A.4 shows the two-convex-edge form is tight only over a sub-region, and
        % A.5's three-convex-edge case has no such form at all -- so cPLQ's Step 1 returns a
        % MINORANT for the first and NOTHING for the second. Splitting the DOMAIN, rather than the
        % envelope, leaves every sub-piece on a path Step 2 already has.
        % splitTightTriangleSym has the derivation, the measurements and the reason it is symbolic.
        %
        % GATED on the piece being exactly x*y. That is what cPLQ's closed forms -- and the split
        % -- are written for; convexEnvelope sends any other indefinite quadratic through xyFrame,
        % and plq_1p.isCanonicalXY is the same test that gate uses. Everything else keeps the
        % vertices it arrives with, exactly.
        %
        % ON BY DEFAULT since 2026-08-18. `CCA2_NO_A45_SPLIT` opts OUT, for comparing against the
        % old behaviour; `CCA2_A45_SPLIT` is still honoured so the tests that set it keep working.
        %
        % WHY IT WAS OPT-IN, and what changed. A.4's cevian foot is IRRATIONAL, so a split
        % sub-triangle has SURD coordinates and every symbolic operation downstream works in a
        % quadratic extension instead of the rationals. The split itself is cheap -- 20 ms on the
        % no-split path, 0.3 s (A.4) / 1.2 s (A.5) when it fires -- but what followed was not:
        % `testcPLQ` went 1542 s to 4728 s, AND `testcPLQ/testRectBiconj` ERRORED. Trading a
        % documented, loud failure on one domain shape for a new one on another is not a trade to
        % make silently, so it stayed off.
        %
        % BOTH OBJECTIONS ARE GONE, measured 2026-08-18 after Step 2 was made exact (three double
        % leaks: `domain.mE`/`cE`, `region.limitOfFAtVertices`, `plq_1p.quadPartsOf` +
        % `conjConvexOverPiece`) and Step 3's merge was repaired:
        %   * `testcPLQ` 8 passed / 0 failed in 2273 s, against 4728 s -- and testRectBiconj is
        %     one of the eight. That exception was a casualty of the double leaks, not a defect
        %     of the split; nothing in the test or the split was changed to fix it.
        %   * assembling f* for the general quadrilateral: 86 cells and 73 min -> 60 cells and
        %     43 min.
        % 2273 s against 1542 s off is a bucket question, not a blocker -- testcPLQ is already in
        % the slow bucket and finishes well inside its timeout. The standing rule is in
        % DECISIONS.md (2026-08-17): a computation has to be CORRECT even when it is slow, and a
        % slow correct path gets its test moved down a bucket rather than traded away.
        %
        % What the split buys is a general convex quadrilateral coming out EXACT -- 10 of 10
        % probe points, 8 of 8 through the full assembly -- where the default used to raise
        % MATLAB:badsubscript on a 3-convex-edge triangle and return a MINORANT in place of the
        % envelope on a 2-convex-edge one.
            if ~isempty(getenv('CCA2_NO_A45_SPLIT')) || ~obj.isBilinear
                ps = [ps, plq_1p(domain(t, vars(1), vars(2)), obj.f)];
                return
            end
            sub = splitTightTriangleSym(t);
            if numel(sub) < 2
                ps = [ps, plq_1p(domain(t, vars(1), vars(2)), obj.f)];
                return
            end
            for k = 1:numel(sub)
                ps = [ps, plq_1p(domain(sub{k}, vars(1), vars(2)), obj.f)]; %#ok<AGROW>
            end
        end

        function tf = isBilinear(obj)
        % Is this piece's function EXACTLY x*y -- no linear part, no constant? Tested rather than
        % assumed, because the piece reaching triangulate can also carry a RATIONAL function fed
        % back in by ratPolToPlq, whose quadParts throws rather than returning a constant matrix.
            tf = false;
            try
                [Q, L, c] = obj.quadParts;
                tf = plq_1p.isCanonicalXY(Q, L, c);
            catch
                tf = false;
            end
        end
    end

    methods (Static)
        function tf = isCanonicalXY(Q, L, c)
        % Is q ALREADY exactly x*y? cPLQ's closed forms assume that literally -- not "indefinite",
        % not "a multiple of x*y", but q = x*y with no linear or constant part -- so anything
        % else goes through xyFrame. Checking here rather than assuming keeps cPLQ's own inputs
        % (every one of which is x*y) on the untouched path, so the existing suite is unaffected.
            tol = 1e-9;
            tf = norm(double(Q) - [0 1; 1 0]) <= tol && norm(double(L(:))) <= tol && ...
                 abs(double(c)) <= tol;
        end

        function [isAff, a, b, c] = affineParts(sf, vars)
        % Is the symbolicFunction sf affine in vars, and if so what are its coefficients?
        %
        % Coefficients come out by EVALUATION at (0,0), (1,0), (0,1), which recovers an affine
        % expression exactly however it happens to be written -- the same reasoning
        % region.linearForm records. Two guards, both of them live cases rather than caution:
        % a RATIONAL envelope (cPLQ's nCE==1 form) can divide by zero at those points, and a
        % quadratic envelope must be reported non-affine rather than silently truncated to its
        % linear part. The identity is confirmed symbolically, not assumed from the three
        % samples, so a quadratic that happens to agree at all three is still rejected.
            a = sym(0); b = sym(0); c = sym(0); isAff = false;
            try
                if ~sf.isPolynomial
                    return                      % rational: not affine, and unsafe to evaluate
                end
                e = sf.f;
                c = subs(e, vars, [0 0]);
                a = subs(e, vars, [1 0]) - c;
                b = subs(e, vars, [0 1]) - c;
                if ~all(isfinite(double([a b c])))
                    a = sym(0); b = sym(0); c = sym(0);
                    return
                end
                isAff = isAlways(simplify(e - (a*vars(1) + b*vars(2) + c)) == 0, ...
                                 'Unknown', 'false');
            catch
                a = sym(0); b = sym(0); c = sym(0); isAff = false;
            end
            if ~isAff
                a = sym(0); b = sym(0); c = sym(0);
            end
        end
    end

    methods % convex envelope
        function vars = pieceVars(obj)
        % The piece's two variables, in order. obj.f.getVars reports only the variables that
        % actually OCCUR in f, so a piece carrying a function of one variable -- q = -x^2, say,
        % which is a perfectly ordinary face function -- returns a single symbol and every
        % `vars(2)` downstream errors with "Index exceeds the number of array elements". The
        % domain always carries both, in the order the piece was built with, so it is the
        % authority on which plane this piece lives in.
            vars = obj.f.getVars;
            if numel(vars) == 2
                return
            end
            if ~isempty(obj.d) && ~isempty(obj.d.polygon) && numel(obj.d.polygon.vars) == 2
                vars = obj.d.polygon.vars;
                return
            end
            error('plq_1p:pieceVars', ...
                'cannot determine this piece''s two variables from either its function or its domain.');
        end

        function [Q, L, c] = quadParts(obj)
        % This piece's quadratic as q(x) = 1/2 x'Qx + L'x + c. Read by differentiation rather
        % than by matching monomials, so it does not care how the expression is written.
            vars = obj.pieceVars;
            q = obj.f.f;
            Q = double(hessian(q, vars));
            g = double(subs(gradient(q, vars), vars, [0 0]));
            L = g(:);
            c = double(subs(q, vars, [0 0]));
        end

        function k = quadKind(obj)
        % 'convex' | 'concave' | 'indefinite' | 'affine' | 'other', by the SIGNS OF THE
        % EIGENVALUES of Q -- never by cPLQ's nCE, which counts edges of positive finite SLOPE
        % and therefore only classifies f = x*y correctly. This is the dispatch the rest of
        % Step 1 keys on.
        %
        % 'other' is not a defensive afterthought: this piece's f is NOT always a quadratic.
        % ratPolToPlq feeds Step 1's envelope back in as a plq_1p, and for a 2-convex-edge
        % triangle that envelope is RATIONAL -- hessian/double then throw rather than returning a
        % constant matrix. Every caller must let 'other' fall through to cPLQ's own nCE branches,
        % which is what handles a rational face today (conjCPLQTest's
        % indefiniteTriangleTwoConvexEdgesSplitViaCPLQStep2 is exactly that path, and it broke
        % when this method assumed a quadratic).
            k = 'other';
            try
                Q = obj.quadParts;
                if any(~isfinite(Q(:))), return, end
            catch
                return
            end
            tol = 1e-9 * max(1, max(abs(Q(:))));
            ev = eig((Q+Q')/2);
            if max(abs(ev)) <= tol,            k = 'affine';
            elseif min(ev) >= -tol,            k = 'convex';
            elseif max(ev) <=  tol,            k = 'concave';
            else,                              k = 'indefinite';
            end
        end

        function [ok, ax, dx, ay, dy, c, lo, hi] = separableParts(obj)
        % Is this piece SEPARABLE OVER A PRODUCT -- f(x,y) = f1(x) + f2(y) + c with the domain an
        % axis-aligned BOX? Both halves are required and neither implies the other: a separable f
        % over a triangle does not separate (the domain couples the variables), and a box under a
        % cross term does not either.
        %
        % When it holds, conjSeparableOverBox computes f* as two 1-D conjugates and a product,
        % with no 2-D region arithmetic at all. See that file for the derivation and the numbers.
        %
        % Everything here is a decision, so it is taken numerically and conservatively: anything
        % unreadable answers false, and the caller falls back to the general path.
            ok = false; ax = 0; dx = 0; ay = 0; dy = 0; c = 0; lo = [0 0]; hi = [0 0];
            try
                vars = obj.pieceVars;
                [Q, L, cc] = obj.quadParts;
            catch
                return
            end
            if any(~isfinite(Q(:))) || any(~isfinite(L(:))) || ~isfinite(cc)
                return
            end
            % NO CROSS TERM. Q is the Hessian, so Q(1,2) is the coefficient of x*y.
            if abs(Q(1,2)) > 1.0d-12 * max(1, max(abs(Q(:))))
                return
            end
            % THE DOMAIN MUST BE A BOX: every facet affine and axis-aligned, two bounds per axis.
            d = obj.d.polygon;
            if isempty(d), return, end
            [A, b, lin] = d.linearForm;
            if ~all(lin), return, end
            loB = [-inf -inf]; hiB = [inf inf];
            for k = 1:size(A,1)
                r = A(k,:);
                nr = norm(r);
                if nr <= 1.0d-12, continue, end
                r = r / nr; bk = b(k) / nr;
                if abs(r(2)) <= 1.0d-9                      % a bound on x
                    if r(1) > 0, hiB(1) = min(hiB(1), bk/r(1));
                    else,        loB(1) = max(loB(1), bk/r(1));
                    end
                elseif abs(r(1)) <= 1.0d-9                  % a bound on y
                    if r(2) > 0, hiB(2) = min(hiB(2), bk/r(2));
                    else,        loB(2) = max(loB(2), bk/r(2));
                    end
                else
                    return                                  % an oblique facet: not a box
                end
            end
            if any(~isfinite([loB hiB])) || any(hiB <= loB)
                return                                      % unbounded or empty in some axis
            end
            ax = Q(1,1)/2; ay = Q(2,2)/2;                   % q = 1/2 x'Qx + L'x + c
            dx = L(1);     dy = L(2);
            c  = cc;
            lo = loB; hi = hiB;
            ok = true;
        end

        function obj = convexEnvelope(obj)
          vars = obj.pieceVars;
          x = vars(1);
          y = vars(2);

          % An INDEFINITE quadratic that is not already x*y is moved into the frame where it IS
          % x*y, and everything downstream happens there; see xyFrame.m. Without this, Step 1's
          % closed forms -- which never reference obj.f -- compute the envelope of x*y whatever
          % the face carries, and Case C returns wrong values for any f other than x*y (measured:
          % f*(0.3,0.4) = 0.4 where the truth is 0.125). The envelope is deliberately left EMPTY
          % here: conjugate() rebuilds this piece in the z-frame, envelope and all, so computing
          % one in the wrong frame now would only be discarded -- and could raise on the way.
          if strcmp(obj.quadKind, 'indefinite')
              [Q, L, c] = obj.quadParts;
              if ~plq_1p.isCanonicalXY(Q, L, c)
                  [M, a, c0] = xyFrame(Q, L, c);
                  obj.frame = struct('M', M, 'a', a, 'c0', c0);
                  obj.lCE = true;
                  return
              end
          end

          obj = convexEnvelope1 (obj,x,y);
        end

        function obj = convexEnvelope1 (obj,x,y)
            % a=sym('a');
            % b=sym('b');
              
            % etaV : eta functions corresponding to set obj.d.V
            % etaE : eta functions corresponding to set obj.d.E
            % etaR : domain of etaE - stored as etaR(i,1:3) : [function,lb,ub] => lb <= function <= ub 
            
            % [etaV, etaE, etaR] =  getEtaFunctions (obj,x,y);
            % 
            % if obj.lPrintEta
            %   obj.d.V
            %   obj.d.E
            %   disp('etaV')
            %   etaV.printL
            %   disp('etaE')
            %   etaE.printL
            %   disp('etaR')
            %   etaR.printL
            % end

            % CLASSIFY BY THE QUADRATIC, not by nCE. nCE counts edges of positive finite SLOPE,
            % which detects "q is convex along this edge" only when q is x*y -- and by the time
            % control reaches here an indefinite q HAS been reduced to x*y (convexEnvelope sends
            % anything else through xyFrame), so the nCE branches below are exactly as valid as
            % they are in cPLQ. What nCE cannot classify at all is a convex or concave q, which
            % cPLQ never sees and CCA2 does:
            %   convex / affine -> co(q|P) = q on a convex P, so there is no envelope to compute.
            %                      Step 2 conjugates the quadratic directly (conjConvexOverPiece).
            %   concave         -> co(q|P) is affine; convEnvUnbounded builds it from the actual
            %                      values of q, for a triangle, a wedge or a half-strip alike.
            kind = obj.quadKind;
            if strcmp(kind, 'convex') || strcmp(kind, 'affine')
                obj.envelope = [obj.envelope, functionNDomain(obj.f, obj.d.polygon)];
                obj.lCE = true;
                return
            end

            % An UNBOUNDED piece -- a wedge or a half-strip out of fanUnboundedFace -- has no
            % nCE classification to speak of: obj.d.nE is computed by walking consecutive
            % vx/vy pairs, and for such a piece some of those are intmax direction markers, not
            % points. Its envelope is built from the half-planes instead. convEnvUnbounded also
            % owns the -inf decision (region.quadUnboundedBelow), which has no counterpart in
            % the bounded branches below because a bounded face never needs one.
            try
                [~, rkind] = obj.d.polygon.recessionRays;
            catch
                rkind = 'bounded';      % non-affine facet: not a case this branch handles
            end
            if strcmp(kind, 'concave') || ~strcmp(rkind, 'bounded')
                expr = convEnvUnbounded(obj.d.polygon, obj.f.f, [x y]);
                obj.envelope = [obj.envelope, ...
                                functionNDomain(symbolicFunction(expr), obj.d.polygon)];
                obj.lCE = true;
                return
            end

            nCE = obj.d.nE;  %size(etaE,1)


            if nCE == 0
              % HISTORY: this was a closed form in the vertex COORDINATES alone --
              %   a = (x1*y1*y2 - x1*y1*y3 - ...)/(x1*y2 - x2*y1 - ...), and likewise b and c --
              % which never referenced obj.f. That is the affine interpolant through
              % (xi, yi, xi*yi), i.e. it silently computed the envelope of x*y whatever the
              % piece's actual function was; pinned by test, it returned the identical envelope
              % for x*y, x^2-y^2, (x^2+y^2)/2 and 3xy+7x-2y+5. Harmless in cPLQ, whose every
              % caller passes f = x*y, but not in CCA2, where quaPolToPlq builds a general
              % quadratic per face. convEnvUnbounded interpolates the actual values of obj.f and
              % reduces to exactly the old formula when obj.f is x*y. It also tests the
              % "no convex edge" precondition as d'Qd <= 0 per edge, which is what nCE == 0
              % means for a general quadratic (nCE itself tests the SLOPE, correct only for x*y).
              expr = convEnvUnbounded(obj.d.polygon, obj.f.f, [x y]);
              f = symbolicFunction(expr);
              obj.envelope = [obj.envelope,functionNDomain(f, obj.d.polygon)];

              obj.lCE = true;
              return
            end
            if nCE == 1
              
              x2 = obj.d.polygon.vx(2);  
              x3 = obj.d.polygon.vx(3);  
              y1 = obj.d.polygon.vy(1);  
              y2 = obj.d.polygon.vy(2);  
              y3 = obj.d.polygon.vy(3);  

              m =  obj.d.mE(1);
              q =  obj.d.cE(1);
              x1 = obj.d.polygon.vx(obj.d.V(1))  ;
              y1 = obj.d.polygon.vy(obj.d.V(1))  ;

              a = -m*y1;
              b = q;
              c = x1;
              d = -q*y1 + m*x1*y1;
              e = -q*x1 - x1*y1;
              f0 = q*x1*y1;
              g = -m;
              h = 1;
              k = -y1 + m*x1;
              f = symbolicFunction(a*x^2+b*x*y+c*y^2+d*x+e*y+f0,g*x+h*y+k);
              obj.envelope = [obj.envelope,functionNDomain(f, obj.d.polygon)];

              obj.lCE = true;
              return
            end

            if nCE == 2
                m_h = sym('m_h');
                 m_w = sym('m_w');
                 q_h = sym('q_h');
                 q_w = sym('q_w');
                % m_h =  obj.d.mE(1);
                % q_h =  obj.d.cE(1);
                % m_w =  obj.d.mE(2);
                % q_w =  obj.d.cE(2);
                
                a =  (m_h*m_w)/(m_h + m_w + 2*sqrt(m_h*m_w));
                b =  (2*sqrt(m_h*m_w))/(m_h + m_w + 2*sqrt(m_h*m_w));
                c =  (1)/(m_h + m_w + 2*sqrt(m_h*m_w));
                d =  ((m_h*q_w + m_w*q_h))/(m_h + m_w + 2*sqrt(m_h*m_w));
                e =  ((- (q_h + q_w)))/(m_h + m_w + 2*sqrt(m_h*m_w));
                f0 =  (q_h*q_w)/(m_h + m_w + 2*sqrt(m_h*m_w));

                
                f = symbolicFunction(a*x^2+b*x*y+c*y^2+d*x+e*y+f0);
                f = f.subsF([m_h,m_w],obj.d.mE(1:2));
                f = f.subsF([q_h,q_w],obj.d.cE(1:2));
                
                obj.envelope = [obj.envelope,functionNDomain(f, obj.d.polygon)];
                
                obj.lCE = true;
                return
            end
           
             
        end
           % disp('out')

         

        function [etaV, etaE, etaR] = getEtaFunctions (obj,x,y)

            a=sym('a');
            b=sym('b');
            eta = obj.f - symbolicFunction(a*x+b*y);
            
            % Eta for Edges
            for i = 1:obj.d.nE
                xv1 = obj.d.polygon.vx(obj.d.E(i,1));
                yv1 = obj.d.polygon.vy(obj.d.E(i,1));
                
                xv2 = obj.d.polygon.vx(obj.d.E(i,2));
                yv2 = obj.d.polygon.vy(obj.d.E(i,2));
                
                edgey = obj.d.mE(i) * x + obj.d.cE(i);
                etaT = eta.subsVarsPartial([y],[edgey]);
                df = etaT.dfdx(x);

                f0 = obj.f.subsVarsPartial([y],[edgey]);
                df0 = f0.dfdx(x);
                df1 = df0.subsF([x],[xv1]);
                df2 = df0.subsF([x],[xv2]);
                if df1 < df2
                    etaE(i,1) = eta.subsVarsPartial([x,y],[xv1,yv1]);
                    etaE(i,3) = eta.subsVarsPartial([x,y],[xv2,yv2]);
                    etaR(i,2) = df1; 
                    etaR(i,3) = df2; 
                else
                    etaE(i,3) = eta.subsVarsPartial([x,y],[xv1,yv1]);
                    etaE(i,1) = eta.subsVarsPartial([x,y],[xv2,yv2]);
                    etaR(i,3) = df1; 
                    etaR(i,2) = df2; 
                end
                etaR(i,1) = symbolicFunction(a+obj.d.mE(i)*b);
                etaE(i,2) =  symbolicFunction((-(a+obj.d.mE(i)*b-obj.d.cE(i))^2/(4*obj.d.mE(i)))-b*obj.d.cE(i))
               
                
                xp = df.solve(x);
                etaE(i,2) = etaT.subsVarsPartial([x],[xp]);
                etaR(i,1) = symbolicFunction(a+obj.d.mE(i)*b);
                
            end

            % Eta for Vertices
            for i = 1:obj.d.nV
                xv = obj.d.polygon.vx(obj.d.V(i));
                yv = obj.d.polygon.vy(obj.d.V(i));
                etaV(i) = eta.subsVarsPartial([x,y],[xv,yv]);
            end
            if obj.d.nV == 0
                etaV=symbolicFunction.empty();
            end
            if obj.d.nE == 0
                etaE=symbolicFunction.empty();
                etaR=symbolicFunction.empty();
            end
            
        end

    end

    methods % conjugate
        function obj = conjugate (obj)
            % SEPARABLE OVER A BOX -- taken FIRST, because it makes every step below unnecessary.
            % f = f1(x) + f2(y) + c on a product domain conjugates to f1*(s1) + f2*(s2) - c, two
            % 1-D problems in closed form. No envelope, no normal cones, no region arithmetic.
            % This is checked before the frame change on purpose: an indefinite DIAGONAL quadratic
            % such as x^2 - y^2 is separable as it stands, and rotating it into the x*y frame
            % would destroy both the separability and the box. See conjSeparableOverBox.m.
            [okSep, axS, dxS, ayS, dyS, cS, loS, hiS] = obj.separableParts;
            if okSep
                obj.conjugates = conjSeparableOverBox(axS, dxS, ayS, dyS, cS, loS, hiS, ...
                                                      [sym('s_1'), sym('s_2')]);
                obj.conjfia = [1, numel(obj.conjugates)+1];
                obj.lConj = true;
                return
            end

            % FRAME CHANGE. When convexEnvelope decided this piece's quadratic is indefinite but
            % not x*y, the whole computation is redone in the z-frame where it IS x*y -- domain
            % and all -- and the resulting conjugate is read back through
            %       f*(s) = h*(M's - a) - c0
            % (see xyFrame.m for the derivation). Doing it here rather than in convexEnvelope
            % keeps the transform in ONE place: the z-frame piece is an ordinary plq_1p carrying
            % exactly x*y, so it goes down the untouched cPLQ path and needs no special casing.
            if ~isempty(obj.frame)
                % pvars must be a LOCAL: obj.pieceVars(1) parses as the method call
                % pieceVars(obj,1), not as indexing the returned pair, and errors with
                % MATLAB:TooManyInputs.
                pvars = obj.pieceVars;
                objT = plq_1p(transformDomain(obj.d, inv(obj.frame.M), pvars), ...
                              symbolicFunction(pvars(1) * pvars(2)));
                objT = objT.convexEnvelope;
                objT = objT.conjugate;
                obj.envelope   = objT.envelope;      % kept in the z-frame, for printouts only
                obj.conjugates = substituteFrame(objT.conjugates, obj.frame);
                obj.conjfia = [1, size(obj.conjugates,2)+1];
                obj.lConj = true;
                return
            end

            obj.conjfia(1) = 1;
            for i=1:max(1 , size(obj.envelope,2)) %For triangles where convex envelope is not computed
              obj = obj.conjugateFunction(i);
              obj.conjfia(i+1) = size(obj.conjugates,2)+1;
            end
        end

        function obj = conjugateFunction (obj,i)
            
            nCE = obj.d.nE;
            vars = obj.pieceVars;
            s1 = sym('s_1');
            s2 = sym('s_2');
            dualVars = [s1,s2];

            % DISPATCH. nCE = obj.d.nE counts edges of positive finite SLOPE, computed by
            % walking consecutive vx/vy pairs. For an UNBOUNDED piece some of those pairs are
            % +/-intmax direction markers, so the slopes -- and hence nCE -- are meaningless,
            % and the piece would fall into whichever of the nCE==1/2 branches the garbage
            % happened to select. (Measured: a wedge with a convex q landed in a quadratic
            % branch and produced conjugate pieces carrying 2147483647*s_2 and
            % intmax^2 = 4611686014132420609, max error 1.15e18.) What actually decides which
            % branch is correct is the ENVELOPE: an affine envelope is conjugated by the
            % support-function construction below, whatever the domain's shape or boundedness.
            %
            % A BOUNDED piece with nCE ~= 0 is left exactly as it was: it is not probed, not
            % classified, and not routed differently. That matters, because the probe below is
            % an evaluation and the nCE==1 envelope is RATIONAL -- evaluating it at (0,0) can
            % divide by zero, which is a live case in conjCPLQTest, not a hypothetical.
            % Boundedness IS tested for every piece, including nCE==0. nCE is computed from
            % slopes between consecutive vx/vy entries, so for an unbounded piece it is garbage
            % and may perfectly well come out 0 -- keying the test on nCE~=0 would let exactly
            % those pieces slip back onto the vertex-list path this branch exists to avoid.
            % It is the affine PROBE, not this, that must stay off the bounded nCE~=0 pieces.
            envUnbounded = false;
            try
                [~, rkind] = obj.envelope(i).d.recessionRays;
                envUnbounded = ~strcmp(rkind, 'bounded');
            catch
                envUnbounded = false;   % non-affine facet: not a case this branch handles
            end

            % The affine probe runs for EVERY piece now. It used to be restricted to
            % nCE == 0 or unbounded, because evaluating cPLQ's RATIONAL nCE==1 envelope at
            % (0,0) divides by zero -- but affineParts guards that itself (it returns early
            % unless the envelope is a polynomial), so the restriction is no longer buying
            % anything and it was actively wrong: a CONCAVE q has an affine envelope while nCE,
            % which tests edge slopes, can perfectly well report 1. Measured on q = -(x^2+y^2)/2
            % over conv{(0,0),(1,0),(1,1)}: envelope -x/2 - y/2, correct, but nCE = 1 sent it to
            % the rank-1 branch, which returned f*(0.3,0.4) = -0.3 where the truth is 1.7.
            % cPLQ's own inputs are unaffected: for f = x*y the nCE==1 envelope is rational and
            % the nCE==2 envelope is a genuine quadratic, so neither is reported affine.
            [envIsAffine, ea, eb, ec] = plq_1p.affineParts(obj.envelope(i).f, vars);
            % A CURVED CONVEX envelope -- which is what a convex face gets, since co q = q --
            % is conjugated by the KKT active-set decomposition, vertex/edge/interior cells.
            % This is the branch cPLQ has no counterpart for: its Step 1 always hands Step 2 an
            % affine or rank-1-PSD envelope, so a curved convex one never arises there, and
            % conjugateOfPiecePoly correspondingly returns only the vertex cells for one. It
            % covers the bounded triangle and the unbounded wedge/half-strip with the same code,
            % because a ray contributes its direction exactly as a bounded edge does.
            envKind = envelopeKind(obj.envelope(i).f, vars);
            if ~envIsAffine && strcmp(envKind, 'convex')
                [eQ, eL, ec2] = quadPartsOf(obj.envelope(i).f, vars);
                obj.conjugates = [obj.conjugates, ...
                    conjConvexOverPiece(obj.envelope(i).d, eQ, eL, ec2, dualVars)];
                obj.lConj = true;
                return
            end

            if envUnbounded && ~envIsAffine
                error('plq_1p:conjugateFunction:unboundedNonAffine', ...
                    ['this piece is unbounded and its convex envelope %s is neither affine nor ' ...
                     'convex, so neither the support-function construction nor the active-set ' ...
                     'one applies.'], char(obj.envelope(i).f.f));
            end
            if nCE == 0 && ~envIsAffine
                error('plq_1p:conjugateFunction:nonAffineEnvelope', ...
                    ['this piece has no convex edge but its envelope %s is not affine, so the ' ...
                     'support-function construction below does not apply.'], ...
                    char(obj.envelope(i).f.f));
            end

            % An UNBOUNDED piece must have its vertex list put in BOUNDARY CYCLIC ORDER first.
            % getNormalConeVertex, below, walks consecutive vertices j, j+1 and wraps at nv, so
            % it is meaningless on an unordered list -- and getVertices does not order the
            % box-clip vertices it appends for an unbounded region. Unordered, the first
            % quadrant's list came out (0,0), (INF,INF), (0,INF), (INF,0), which makes (0,0)
            % adjacent to (INF,INF) and yields the cone {s1+s2 <= 0, s1 <= 0} instead of
            % {s1 <= 0, s2 <= 0} -- reporting f*(-10,5) = 0 where the truth is +inf.
            %
            % poly2orderUnbounded is the routine for exactly this, and it used to throw here;
            % see its own HISTORY note and region.getEdges'. Ordered, the list is
            % (0,0), (INF,0), (INF,INF), (0,INF) and every FINITE vertex gets its true normal
            % cone. The at-infinity vertices still get intmax-laden cones, which is why the
            % emit loop below skips them: they are box-clip artefacts, not vertices of the sup.
            if envUnbounded
                obj.envelope(i).d = obj.envelope(i).d.poly2orderUnbounded;
            end

            if nCE == 0 || envIsAffine


                NCV = obj.envelope(i).d.getNormalConeVertex(s1, s2);
                [subdV,undV] =  obj.envelope(i).getSubdiffVertexT1 (NCV, dualVars);

                % HISTORY, two defects in one place, both invisible on cPLQ's own inputs.
                %
                % (1) a,b,c used to be RECOMPUTED here by the same closed form in the vertex
                % coordinates that convexEnvelope1's nCE==0 branch used -- the affine
                % interpolant through (xi,yi,xi*yi). So this branch conjugated the envelope of
                % x*y no matter what obj.f was, and it did so even when convexEnvelope1 had
                % just produced a different envelope. The envelope is right there in
                % obj.envelope(i).f; read it. Coefficients come out by EVALUATION at (0,0),
                % (1,0), (0,1), which recovers an affine expression exactly however it happens
                % to be written (the same reasoning region.linearForm records), and stays
                % symbolic so the exact-arithmetic chain is not broken.
                %
                % (2) the loop ran over every one of obj.envelope(i).d.nv vertices. For an
                % UNBOUNDED face some of those entries are the +/-intmax markers standing for
                % recession DIRECTIONS, not points, and each one contributed a bogus affine
                % piece s1*2147483647 + ... and a constant intmax^2 = 4611686014132420609.
                % That is next-step 1(b) of the session handoff, measured there as a max error
                % of 1.15e18. A direction is not a vertex of the sup: it constrains the DOMAIN
                % of the conjugate (via the normal cone), it does not add a piece to it.
                a = ea; b = eb; c = ec;
                for j = 1:obj.envelope(i).d.nv
                  x1 = obj.envelope(i).d.vx(j);
                  y1 = obj.envelope(i).d.vy(j);
                  if abs(x1) == intmax || abs(y1) == intmax
                      continue                  % a recession direction, not a vertex: no piece
                  end

                  conjf = symbolicFunction(s1 * x1 + s2 * y1  - (a*x1+b*y1+c));
                  conjd = region(subdV(j,:), dualVars);
                  obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                end

              obj.lConj = true;
              return
            end
            
            if nCE == 1
                
              m =  obj.d.mE(1);
              q =  obj.d.cE(1);
              x1 = obj.d.polygon.vx(obj.d.V(1))  ;
              y1 = obj.d.polygon.vy(obj.d.V(1))  ;

              a = -1;
              b = -2*m;
              d = 2*q+4*m*x1;
              c = -m^2;
              e = -(2*m*q - 4*m*y1);
              f = -(q^2 + 4*m*x1*y1);
             

             
              crs = a*s1^2+b*s1*s2 + c* s2^2 + d*s1 + e*s2 + f;   % nonconvex part


              NCV = obj.envelope(i).d.getNormalConeVertex(s1, s2);
              [subdV,undV] =  obj.envelope(i).getSubdiffVertexT1 (NCV, dualVars);
              subdV = obj.envelope(i).getSubDiffVertexSpT1(subdV, undV, -crs);


              NCE = obj.envelope(i).d.getNormalConeEdge(s1, s2);
              [subdE,unR] = obj.envelope(i).getSubdiffVertexT2 (NCE, dualVars);
              
              %obj.envelope(i).d.print

              edgeNo = obj.envelope(i).d.getEdgeNos(vars);
             
              [subdE, unR, crs] = obj.envelope(i).getSubDiffEdgeT1(subdE, edgeNo, undV, crs, dualVars);
              
              
              %expr = obj.envelope(i).conjugateExprVerticesT1 (dualVars, undV )

                x1 = obj.d.polygon.vx(1);  
                x2 = obj.d.polygon.vx(2);  
                x3 = obj.d.polygon.vx(3);  
                y1 = obj.d.polygon.vy(1);  
                y2 = obj.d.polygon.vy(2);  
                y3 = obj.d.polygon.vy(3);  
                a = ((x1*y1*y2 - x1*y1*y3 - x2*y1*y2 + x2*y2*y3 + x3*y1*y3 - x3*y2*y3))/((x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
                b = ((x1*x2*y2 - x1*x2*y1 + x1*x3*y1 - x1*x3*y3 - x2*x3*y2 + x2*x3*y3))/((x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
                c = ((x1*x2*y1*y3 - x1*x3*y1*y2 - x1*x2*y2*y3 + x2*x3*y1*y2 + x1*x3*y2*y3 - x2*x3*y1*y3))/((x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2));
                for j = 1:obj.envelope(i).d.nv
                  if undV(j)
                      if j == obj.envelope(i).d.nv
                        e0 = 1;
                      else    
                        e0 = j+1;
                      end  
                      x1 = obj.d.polygon.vx(j);   % same as obj.envelope(i).d.vx(1) for triangles
                      y1 = obj.d.polygon.vy(j);  
                  
                      conjf = symbolicFunction(s1 * x1 + s2 * y1  - (a*x1+b*y1+c));
                      conjd = region([subdE(e0,1:2),-subdE(e0,3)], dualVars);
                      obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                    %%%%%%%%%%%%%
                      r = subdV(e0,:);
                      r(2) = -r(2);
                       conjd = region(r, dualVars);
                       obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                       %%%%%%%%%%%%%
                       if j == 1
                        e0 = obj.envelope(i).d.nv;
                       else    
                        e0 = j-1;
                       end  
                       r = subdV(e0,:);
                       r(1) = -r(1);
                       conjd = region(r, dualVars);
                       obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                  else



                      x1 = obj.d.polygon.vx(j);   % same as obj.envelope(i).d.vx(1) for triangles
                      y1 = obj.d.polygon.vy(j);  
                        
                      conjf = symbolicFunction(s1 * x1 + s2 * y1  - (a*x1+b*y1+c));
                      conjd = region(subdV(j,:), dualVars);
                      obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                  end
                end
                


               for j = 1:obj.envelope(i).d.nv

                  if unR(j)
                       a = 1/(4*m);
                       b = 1/2;
                       c = m/4;
                       d = -q/(2*m);
                       e =  q/2;
                       f = q^2/(4*m);
                       conjf = symbolicFunction(a*s1^2+b*s1*s2+c*s2^2+d*s1+e*s2+f);
                       % obj.envelope(1).f
                       % edge = vars(2) - m*vars(1) -q
                       % edge
                       % conjugateExpr (edge, obj.envelope(1).f.f,vars(1),vars(2))
                       % pause
                       conjd = region(subdE(j,:), dualVars);
                       obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                  
                  end
               end
               obj.lConj = true;
               return
            end
            if nCE == 2
                
              NCV = obj.envelope(i).d.getNormalConeVertex(s1, s2);
              [subdV,undV] =  obj.envelope(i).getSubdiffVertexT1 (NCV, dualVars);
              

               [coef,terms] = coeffs(obj.envelope(i).f.f,obj.envelope(i).d.vars);
               x = obj.d.polygon.vars(1);
               y = obj.d.polygon.vars(2);
               for j = 1:obj.envelope(i).d.nv
                 x1 = obj.d.polygon.vx(j); 
                 y1 = obj.d.polygon.vy(j);  
                 % change this for direct computation
                 conjf = symbolicFunction(s1 * x1 + s2 * y1)  - obj.envelope(i).f.subsF([x,y],[x1,y1]);
                 %subdV(j,:)
                 conjd = region(subdV(j,:), dualVars);
                 %conjd.print
                 if isempty(conjd)
                     continue
                 end
                 obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
                end
                %pause
              % HISTORY: a..f are read off the envelope polynomial by MATCHING monomials, so a
              % monomial that is ABSENT left its coefficient undefined and `grad` below then
              % errored with "Unrecognized function or variable 'd'". That is not exotic: the
              % envelope of a triangle with no linear part -- e.g. CCA2's Step 1 output for the
              % 3-convex-edge triangle conv{(0,0),(1,1),(3,2)}, whose first face is a pure
              % quadratic form alpha*xy + beta*x^2 + gamma*y^2 -- has no x, y or constant term at
              % all. An absent monomial has coefficient zero, which is exactly what the loop below
              % means to leave in place, so seed all six here.
              a = 0; b = 0; c = 0; d = 0; e = 0; f = 0; %#ok<NASGU> (c,f used only by the
              % commented-out diagnostics below; seeded anyway so the set stays consistent)
              for j = 1:size(terms,2)
                  if isAlways (terms(j) == x^2)
                      a = coef(j);
                  end
                  if isAlways (terms(j) == x*y)
                      b = coef(j);
                  end
                  if isAlways (terms(j) == y^2)
                      c = coef(j);
                  end
                  if isAlways (terms(j) == x)
                      d = coef(j);
                  end
                  if isAlways (terms(j) == y)
                      e = coef(j);
                  end
                  if isAlways (terms(j) == 1)
                      f = coef(j);
                  end
              end
              %obj.envelope(i).d.ineqs(1)
              grad = b*s1 - 2*a*s2 - b*d + 2*a*e;

              NCE = obj.envelope(i).d.getNormalConeEdge(s1, s2);
              [subdE,unR] = obj.envelope(i).getSubdiffVertexT2 (NCE, dualVars); 
              
              % check signs and then put gradient ineq properly
              % simplifyFraction(grad)
              % simplifyFraction(grad)
              % simplifyFraction(subdE(3,1))
              % simplifyFraction(subdE(3,2))
              % pause
              % obj.envelope(i).d.ineqs(1)
              % obj.envelope(i).d.ineqs(2)
              % obj.envelope(i).d.ineqs(3)
              [c0,t0] = (coeffs(simplifyFraction(grad),[s1,s2]));
              %double(c0)
              %pause
              % HISTORY: subdE(:,3) -- which half-plane of `grad` belongs to which edge --
              % used to be the hard-coded pattern subdE(1,3)=grad, subdE(2,3)=-grad,
              % subdE(3,3)=grad. That is a statement about one particular edge ORDERING, not
              % about the geometry, so it silently handed edges the wrong half-plane for any
              % other ordering. It is now derived per edge inside the loop below; see the
              % comment there.
 %            for j = 1:obj.envelope(i).d.nv
 % %                 j
 %  %                obj.envelope(i).d.ineqs(j)
 %                  d1(j) = 2*a*obj.envelope(i).d.vx(j) + b* obj.envelope(i).d.vy(j) +d;
 %                  d2(j) = b*obj.envelope(i).d.vx(j) + 2*c* obj.envelope(i).d.vy(j) +e;
 % 
 % 
 % 
 %            end
 %            double(d1)
 %            double(d2)
              for j = 1:obj.envelope(i).d.nv
 %                 j
                    coeffs0 = obj.envelope(i).d.ineqs(j).getLinearCoeffs(obj.envelope(i).d.vars);

                    % if coeffs0(2)/coeffs0(1) > 0
                    %     subdE(j,3) = -grad;
                    % else
                    %     subdE(j,3) = grad;
                    % end

                    % if j < obj.envelope(i).d.nv
                    %    double(subs(grad,[s1,s2],[(d1(j)+d1(j+1))/2,(d2(j)+d2(j+1))/2]))
                    % else
                    %    double(subs(grad,[s1,s2],[(d1(j)+d1(1))/2,(d2(j)+d2(1))/2]))
                    % end
                    % pause
  %                obj.envelope(i).d.ineqs(j)
                  
                  c0 = obj.envelope(i).d.ineqs(j).getLinearCoeffs ([x,y]);

                  % Which half-plane of `grad` this edge owns, from the geometry rather than
                  % from j. In this branch the envelope is a rank-1 PSD quadratic, so its
                  % Hessian [2a b; b 2c] has kernel kdir = [b, -2a], and
                  %   grad = <s - L, kdir>   (L = [d,e] the linear part),
                  % i.e. the objective s'x - f(x) is AFFINE along kdir with slope grad. The sup
                  % is therefore pushed along +kdir when grad > 0, so it can only land on an
                  % edge whose OUTWARD normal nj = [c0(1),c0(2)] satisfies <kdir,nj> > 0; that
                  % edge owns {grad >= 0}, i.e. the inequality -grad <= 0, and an edge with
                  % <kdir,nj> < 0 owns {grad <= 0}.
                  % An edge PARALLEL to kdir (<kdir,nj> = 0) is a level edge of f: the sup is
                  % attained in its relative interior only on the measure-zero line {grad = 0},
                  % where its two endpoints give the same value, so the (closed) vertex regions
                  % already cover it and the edge contributes no region at all.
                  kdir = double([b, -2*a]);
                  nj   = double([c0(1), c0(2)]);
                  tj   = kdir*nj';
                  if abs(tj) <= 1.0d-10*max(1, norm(kdir)*norm(nj))
                      continue
                  elseif tj > 0
                      subdE(j,3) = -grad;
                  else
                      subdE(j,3) = grad;
                  end

                  m = -c0(1)/c0(2);   % put check for zero d
                  q = -c0(3)/c0(2);
% conj = (mh^2*y^2 + 2*mh*qh*y + 2*mh*x*y + qh^2 - 2*qh*x + x^2)/(4*mh)


                  av = 1/(4*m);
                   bv = (2*m)/(4*m);
                   cv = (m^2)/(4*m);
                   dv = (-(2*q))/(4*m);
                   ev = ((2*m*q))/(4*m);
                   fv = (q^2)/(4*m);

                   if m > 0
                    conjf = symbolicFunction(av*s1^2+bv*s1*s2+cv*s2^2+dv*s1+ev*s2+fv);
                   else
                    conjf =  symbolicFunction( conjugateExpr(obj.envelope(i).d.ineqs(j).f,obj.envelope(i).f.f,x,y));
                    conjf = conjf.subsF([x,y],[s1,s2]);
                   end
                   %obj.envelope(i).d.ineqs(j).f
                   %obj.envelope(i).f.f
                   %conjf =  symbolicFunction( conjugateExpr(obj.envelope(i).d.ineqs(j).f,obj.envelope(i).f.f,x,y))
                   %pause
                   %conjf = conjf.subsF([x,y],[s1,s2]);
                   %pause
                   conjd = region(subdE(j,:), dualVars);
                   
                   conjd = conjd.simplifyUnboundedRegion;
                  % conjd.print
                  % pause
                    if isempty(conjd)
                     continue
                 end
                   obj.conjugates = [obj.conjugates,functionNDomain([conjf],conjd)];
              end
              %pause
            end  

               obj.lConj = true;
              return
              
            end
        end

    

    methods % max

        function obj = maximumConjugate(obj)
          for k = obj.conjfia(1):obj.conjfia(2)-1 
             obj.maxConjugate(k) = obj.conjugates(k);
          end
          for i = 2:size(obj.envelope,2) %size(obj.conjfia,2)-1
              obj.maxConjugate = obj.maxConjugate * obj.conjugates(obj.conjfia(i):obj.conjfia(i+1)-1);
              %obj.maxConjugate.printM2
              %obj.maxConjugate.printL
             obj.maxConjugate = obj.maxConjugate.maximumP(true);
          end
          obj.lMConj = true;
        end
    end


end


% ==================================================================================================
function [Q, L, c] = quadPartsOf(sf, vars)
% The quadratic parts of a symbolicFunction, by differentiation rather than monomial matching.
%
% SYMBOLIC, not double. These three go straight into conjConvexOverPiece, which builds a whole
% conjugate cell out of them; taking a double here means two cells that share a facet can end up
% carrying two DIFFERENT doubles of the same exact number, and then region.merge's facet test
% cannot match them. Measured 2026-08-17: two cells sharing a facet carried
% 659536895553805/562949953421312 and 5276295164430439/4503599627370496 -- both `4 - 2*sqrt(2)`,
% one ULP apart. DECISIONS.md has it.
    q = sf.f;
    Q = hessian(q, vars);
    g = subs(gradient(q, vars), vars, [0 0]);
    L = g(:);
    c = subs(q, vars, [0 0]);
end

function k = envelopeKind(sf, vars)
% 'convex' | 'concave' | 'indefinite' | 'affine' | 'other' for an ENVELOPE expression. 'other'
% covers anything not a polynomial of degree <= 2 -- notably cPLQ's nCE==1 envelope, which is
% RATIONAL and must not be handed to a quadratic classifier.
    k = 'other';
    try
        if ~sf.isPolynomial, return, end
        Q = double(hessian(sf.f, vars));
        if any(~isfinite(Q(:))), return, end
    catch
        return
    end
    tol = 1e-9 * max(1, max(abs(Q(:))));
    ev = eig((Q+Q')/2);
    if max(abs(ev)) <= tol,   k = 'affine';
    elseif min(ev) >= -tol,   k = 'convex';
    elseif max(ev) <=  tol,   k = 'concave';
    else,                     k = 'indefinite';
    end
end
