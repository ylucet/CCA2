classdef AlgCon < Alg & Conic
   % AlgCon: the ROOT of the type lattice, the ONE place the mesh data is declared, and the type
   % `biconj` returns when a face of the convex envelope has no formula.
   %
   % It takes over the first two roles from RatCon, whose header stated them and still explains
   % the diamond, the constructor protocol and the no-stored-flag rule. What follows is only what
   % CHANGED and why.
   %
   % THE FUNCTION AXIS GAINED A TOP LEVEL. It is now a three-level chain, most special first:
   %
   %       Qua      <      Rat       <      Alg
   %       den = 1         numerator over   a face is cut out IMPLICITLY, by an integer
   %                       a linear den     polynomial p(x,y,z) = 0 plus a BRANCH INDEX
   %                                        naming which real root in z is the value
   %
   % so the lattice is a 3 x 3 product of function level against subdivision level
   % (Pol < Par < Conic), and this class is its top cell:
   %
   %                 AlgCon  (algebraic on conic)      < Alg & Conic   <- this class
   %                    |
   %                 RatCon  (rational on conic)       < AlgCon & Rat
   %                    |
   %                 RatPar  (rational on parabolic)   < RatCon & Par
   %                 /    \
   %          RatPol        QuaPar                     < RatPar & Pol / RatPar & Qua
   %                 \    /
   %                 QuaPol                            < RatPol & QuaPar
   %
   %       QuaCon      quadratic on a CONIC subdivision            < RatCon & Qua
   %       QuaParCPLQ  a QuaPar still in cPLQ's own symbolic form  < RatPar & Qua
   %
   % WHY THIS ROOT IS CONCRETE, where RatCon was abstract. The concrete cells are named
   % <Function><Subdivision> -- QuaPol, RatPol, QuaPar, QuaCon -- and at (Alg, Conic) the root of
   % the lattice and the concrete cell ARE THE SAME CELL, because `biconj` actually produces one.
   % One class is therefore honest where two names would be a spelling. Probed before committing
   % (`.claude/probeAlg/probeRun2.m`): an abstract class MAY inherit a concrete one, and the root's
   % own constructor does NOT clobber the leaf, so the CONSTRUCTOR PROTOCOL below is unaffected.
   %
   % WHY THE LEVEL EXISTS AT ALL -- see Alg.m, which carries the measurement. In one line: the
   % convex envelope's bridge face is affine PLUS A NORM, which neither Qua nor Rat can hold, but
   % sqrt(Q) satisfies z^2 = Q, so written implicitly the same face is a quadric with RATIONAL
   % coefficients and every stored number stays an integer.
   %
   % THE PAYLOAD IS ONE THING PER FACE, NOT FOUR KINDS. The 2026-08-23 plan gave a face `fkind` in
   % {QUAD, RAT, RULED, AFFINE}. The implicit form SUBSUMES the AFFINE case -- a dual vertex p of
   % degree 4 satisfies its own minimal polynomial over Q, so eliminating it leaves a rational
   % relation in (x,y,z) -- so a face carries `FP{k}` (the polynomial) and `FB(k)` (the branch), or
   % neither, in which case `fN`/`fD` (or `f`/`den`) hold it exactly as before. DECISIONS.md
   % 2026-09-13, and that entry is also why the type is AlgCon and not `AlgAlg`: the bridge REGION
   % was measured to be an intersection of sign conditions on CONICS, so only the function axis
   % grows.
   %
   % PINNED, LIKE EVERY OTHER AXIS-VARYING PROPERTY. `den` is general on a Rat and pinned to 1 on a
   % Qua; `Ec` is general on a Par and pinned to 0 on a Pol; `FP`/`FB` are general on an Alg and
   % PINNED TO EMPTY on a Rat. All three are enforced by `set.` validators that read the object's
   % own TRAITS, from one definition site, so no type can be made to lie about itself -- including
   % types that do not exist yet.
   %
   % WHAT IS *NOT* CHANGED. Adding a level ABOVE the previous top changes no existing type's
   % behaviour: every object that was a RatCon before is one now, is additionally an Alg (a
   % rational function IS algebraic), and cannot carry a payload it could not carry before. The
   % four legacy types and QuaCon are untouched; `Rat < Alg` is the whole of the edit they see.

   properties
        % Everything the lattice shares, declared here ONCE -- the nine mesh properties AND the two
        % axis-varying ones (`den`, `Ec`). No child declares any property of its own; see the
        % "DATA LIVES HERE, ONCE" note above for why that is forced rather than merely tidy.
        nv {mustBeInteger,mustBeNonnegative}%number of vertices
        ne {mustBeInteger,mustBeNonnegative}%number of edges
        nf {mustBeInteger,mustBeNonnegative}%number of faces
        V (:,2){mustBeNumeric} % nv x 2 matrix storing unique vertices
        E (:,3){mustBeInteger,mustBeNonnegative} % ne x 3 matrix storing edge indices where Edge j is [V(E(j,1)), V(E(j,2))]
            %and E(j,3) = 1 for a segment, 0 for a ray, and 2 for a CLOSED conic edge -- a whole
            %ellipse, which has NO endpoints, so E(j,1:2) must be 0 there. See
            %AlgCon.assertEdgeFlags for the convention and why the third value exists; a closed
            %edge bounds the bridge face domain (DECISIONS.md 2026-09-06). Nothing in the toolbox
            %produces one yet, and every consumer whose edge handling assumes two endpoints calls
            %AlgCon.assertNoClosedEdges so the refusal is named rather than silent.
        f (:,10){mustBeNumeric} % nf x 10 matrix storing the coefficients of the cubic
            %f(k,:)=c means at the point (x,y) in the face k  the function value
            % is C*([x^3; x^2 y; x y^2; y^3; x^2; x y; y^2; x; y; 1] .* [1/6; 1/2; 1/2; 1/6; 1/2; 1; 1/2; 1; 1; 1]
            % the coefficients vector [1/6; 1/2; 1/2; 1/6; 1/2; 1; 1/2; 1; 1; 1] is required to easily manipulate
            % Hessians; it is the same reason we work with 0.5 x^2 instead of x^2, i.e. so that diff(0.5 x^2,1) = 1
            % For RatPol this is the NUMERATOR; see RatPol's own `den`.
        F (:,2){mustBeInteger,mustBeNonnegative} % ne x 2 matrix storing indicating function indices. F(j)=[k1, k2] means
            %the quadratic on the left (resp. right) of edge j has coefficients f(k1) (resp. f(k2)).
            %an index value of 0 indicates the function has value +infinity on that side
        % Special cases: functions with domain dimension less than 2
            %the needle function is stored as V = [0 0];E = [];f = [1 2 3];F = [];%V has 1 row, f has 1 row (and should
            %       be a constant but any cubic is accepted
            %the edge chain function
            %
        P %{cell} %cell array of size nf representing an adjacency list. Each
            %element P{i} is an array of indices k with 1<= abs(k) <= ne(i)
            %where ne(i) is the number of edges in face i.
            %If k>0, then the face lies on the right of the edge; otherwise k<0.
            %The edge indices are ordered to obtain a unique
            %representation. If the face is bounded, the smallest index
            %is on the left, then go clockwise. If the face is
            %unbounded, the first index is the unbounded edge with the
            %smallest index, then go clockwise (the last index is the
            %remaining unbounded edge).
        dom {struct}% struct that stores the domain information (set on which the function is finite); PLC assumes the domain is CONNECTED
            %   dim: nonnegative integer storing the dimension of the domain
            %       dim=0 means the domain is a single vertex (only row in V); PLC is a needle function
            %       dim=1 has 3 cases: segment (2 vertexes, 1 edge), ray (2 vertexes, 1 edge), and boundary of a slice
            %       (3 vertexes, 2 edges). COULD ALSO BE A CHAIN OF EDGES OR ANY GRAPH WITH FUNCTION ONLY DEFINED
            %       ON EDGES; NOT TESTED SINCE THE CLASS FOCUSES MOSTLY ON CONVEX FUNCTIONS
            %   P: if nonempty, the domain has dimension two; in that case P stores the indexes of the face
            %       representing the domain with the INVERSE convention as P above to obtain a unique representation.
            %   isConvex: boolean. true if the domain is convex

        % ---- a CALLER'S ASSERTION, not a computed fact ---------------------------------------
        fIsConvex {mustBeNumericOrLogical} = []   % [] = unknown (the default), true = the CALLER
            % asserts the function is convex, false = asserts it is not. It is an assertion about
            % f AS AN EXTENDED-REAL FUNCTION, i.e. including the indicator of its domain.
            %
            % WHY IT IS TRUSTED. Convexity of a SINGLE quadratic piece is free to determine (the
            % eigenvalues of a 2x2 Hessian), and the operators do determine it -- no flag needed.
            % Convexity of a MULTI-PIECE function is not: it needs per-piece convexity AND the
            % gradient jump across every shared edge to be consistent (isFaceConvex + isEdgeConvex,
            % the latter still marked untested). The flag exists for exactly that case, and it
            % unlocks the largest short-circuit in the toolbox -- biconj of a convex f IS f, with
            % nothing to compute.
            %
            % WHAT IS STILL CHECKED. The free NECESSARY condition: every piece's Hessian must be
            % positive semidefinite. A flag contradicting that is refused LOUDLY rather than
            % trusted, because the failure it would otherwise cause is silent -- biconj returning
            % a non-convex f as its own convex envelope. See ALGORITHM.md.
        % ---- the two AXIS-VARYING properties -------------------------------------------------
        % Each is general on its own axis and pinned to a trivial value on the specialization; the
        % pinning is enforced by the set validators below, keyed on the object's own traits.
        den (:,3){mustBeNumeric} % nf x 3 per-face linear DENOMINATORS: den(k,:)=[g h k0] gives the
            % denominator g*x + h*y + k0 of the rational function on face k, so the value on face k
            % is (f(k,:) in the cubic basis) / den(k,:). PINNED to [0 0 1] (i.e. 1) on every face
            % when the object is a Qua -- a polynomial is exactly a rational with unit denominator.
        Ec (:,6){mustBeNumeric} % ne x 6 per-edge conic: Ec(j,:)=[a b c d e f] for the curve
            % a x^2 + b xy + c y^2 + d x + e y + f = 0, which must be a PARABOLA (b^2-4ac=0).
            % An all-zero row means edge j is the straight line through its endpoints. PINNED to
            % all-zero when the object is a Pol -- a polyhedral subdivision is exactly a parabolic
            % one with no curvature.

        % ---- the ALGEBRAIC payload: general on an Alg, PINNED TO EMPTY on a Rat --------------
        EcA (:,1) struct = emptyEcA()   % ne_alg x 1 ALGEBRAIC edges, the fourth subdivision
            % level (see Crv.m). Entry j has fields FP, FB, num, den, and stands for the curve
            %       { x : q(x) = value_FB(x) },   q = num . basis / den,
            % where value_FB is the FB-th real root of FP(x,y,.) -- the boundary between a BRIDGE
            % face and a quadratic one. The polynomial FP(x,y,q(x,y)) that vanishes there is
            % rational but names the PAIR of branches, so the branch index is part of the edge and
            % not decoration; `algCurveSide` decides a side exactly, by counting roots.
            % PINNED to empty unless the object is a Crv -- a conic subdivision has no such edge.
            % An FC row selects this table with a fourth column of 1; see FC's own note.
        FP (:,1) cell = {}   % nf x 1. FP{k} is face k's defining polynomial p_k(x,y,z), as ratQ's
            % dense polyN array in three variables: FP{k}(i+1,j+1,m+1) is the INTEGER coefficient
            % of x^i y^j z^m. EMPTY where face k has a formula, in which case fN/fD (or f/den)
            % hold it and nothing here applies. Measured degrees run from (2,2,2) for an isotropic
            % bridge to (6,6,4) for an anisotropic one, with integer coefficients under 1e9
            % (DECISIONS.md 2026-09-06), so dense storage is at most a few hundred numbers a face.
        FB (:,1){mustBeInteger} = []   % nf x 1 BRANCH INDEX: which real root of FP{k}(x,y,.) in z
            % is the function's value, counted WITH multiplicity and 1-based from an END: POSITIVE
            % counts up from the smallest real root, NEGATIVE down from the largest (-1 is the
            % largest). 0 exactly where FP{k} is empty. See algBranch.m for why an index is a
            % legitimate NAME for an irrational value -- the argument QuaCon.Vname already rests on
            % one dimension down -- and for why the END matters: a colliding pair of roots can
            % LEAVE the real line inside a face, which shifts every index above it and none below.
            % A bridge face is a SUPREMUM, so its value is the largest real root and it is stored
            % as -1; measured over 4758 region points, its ascending index takes two different
            % values on ONE face, so the ascending spelling could not name it. DECISIONS.md
            % 2026-09-17. Positive indices keep their original meaning unchanged.
        src = []             % PROVENANCE ONLY, never the definition of a value: the QuaCon f* this
            % envelope was assembled from, so a face can say which dual cell it came from -- the
            % paper's classification table wants exactly that, and so does every debugging session.
            % A value is ALWAYS recoverable from FP/FB alone and nothing reads src to evaluate one.
        FR (:,2){mustBeInteger,mustBeNonnegative} = []   % PROVENANCE ONLY: nf x 2, the [dual
            % vertex, adjacent face] of `src` that face k came from, or 0 where unknown.
   end

   % ---------------------------------------------------------------------------------------
   % CONSTRUCTOR PROTOCOL -- leaf-only initialization. READ BEFORE EDITING ANY CONSTRUCTOR.
   %
   % Rule: every constructor in this hierarchy must have a NO-ARGUMENT path that writes NOTHING
   % and returns immediately, and all state-writing must happen in the LEAF constructor after its
   % superclass constructors have run.
   %
   % Why: with multiple inheritance MATLAB re-runs a shared base constructor ONCE PER INHERITANCE
   % PATH. For the planned lattice (QuaPol < RatPol & QuaPar, both < RatPar) constructing a QuaPol
   % invokes RatPar's constructor TWICE:
   %       RatPar -> RatPol -> RatPar -> QuaPar -> QuaPol
   % That is not merely wasteful: the SECOND base invocation OVERWRITES whatever the first
   % inheritance path already wrote. Verified directly -- with a mid class writing a field after
   % its own super-call, the constructed leaf ended up with the BASE's value, silently discarding
   % the mid class's:
   %       RatPol alone : f = set-by-RatPol
   %       QuaPol       : f = set-by-RatPar        <- clobbered
   % MATLAB offers no virtual-inheritance escape (a child cannot redefine an inherited property:
   % MATLAB:class:RedefinedProperty), so the protocol above is the fix, not a workaround.
   %
   % Consequence for callers: `QuaPol()`, `QuaPar()`, `RatPol()` with no arguments now return a
   % blank object instead of raising an error. That is also what MATLAB itself needs for
   % `ClassName.empty`, object arrays and `load`, so it is a gain rather than a concession.
   % ---------------------------------------------------------------------------------------

   methods (Static)
      function E = assertEdgeFlags(E, who)
      % objective: validate the segment/ray/CLOSED flag in column 3 of an edge table.
      % [input]  E : ne x 3 edge table; who : (optional) name used in the error message
      % [output] E : unchanged, so this can wrap an expression
      %
      % THE THIRD VALUE, AND WHY IT HAD TO EXIST (2026-09-06). Column 3 has always been 1 for a
      % SEGMENT and 0 for a RAY, and both spellings name their endpoints in columns 1 and 2. A
      % CLOSED conic edge -- a whole ellipse -- has NO endpoints, and the bridge face's own domain
      % is bounded by two of them: the isotropic bridge region is an ANNULUS between two concentric
      % circles (DECISIONS.md 2026-09-06). There is nothing to put in columns 1 and 2, and there is
      % no way to NAME such a point either, because `Vname` names a vertex as the meeting of TWO
      % edge conics and a lone closed curve meets none.
      %
      % So the convention is:
      %
      %       E(j,3) = 1   segment, endpoints E(j,1:2)
      %       E(j,3) = 0   ray, apex E(j,1) and direction through E(j,2)
      %       E(j,3) = 2   CLOSED, no endpoints, and E(j,1:2) MUST be 0
      %
      % ONE MORE VALUE IN A COLUMN THAT ALREADY EXISTED, which is the same move `FC`'s side column
      % made when it gained 0 for "ON the curve" to let a face be thin -- QuaCon.m's header records
      % that choice and the reason: it costs one value rather than a new class.
      %
      % WHY THIS VALIDATOR IS THE WHOLE POINT. Every consumer in the toolbox reads the flag as
      % `E(j,3) == 0` (a ray) or `~= 0` (assumed a segment), so a 2 arriving unannounced is read as
      % a SEGMENT and the routine then indexes vertex 0. Nothing today produces a 2, so this fires
      % on nothing and costs nothing -- it exists so that when a bridge face is first stored, the
      % failure is a named refusal at the boundary instead of a silent wrong answer deep inside.
         if nargin < 2, who = 'edge table'; end
         if isempty(E), return, end
         bad = ~ismember(E(:,3), [0 1 2]);
         if any(bad)
            error('RatCon:badEdgeFlag', ...
               ['%s: edge %d has flag %g in column 3. The only values are 1 (segment), ' ...
                '0 (ray) and 2 (CLOSED, no endpoints).'], who, find(bad,1), E(find(bad,1),3));
         end
         cl = E(:,3) == 2;
         if any(cl & any(E(:,1:2) ~= 0, 2))
            j = find(cl & any(E(:,1:2) ~= 0, 2), 1);
            error('RatCon:closedEdgeHasEndpoints', ...
               ['%s: edge %d is CLOSED (flag 2) but names endpoints [%g %g]. A closed conic has ' ...
                'no endpoints and no vertex can be named on it, so columns 1 and 2 must be 0.'], ...
               who, j, E(j,1), E(j,2));
         end
      end

      function assertNoClosedEdges(E, who)
      % objective: refuse, BY NAME, an edge table containing a CLOSED edge.
      % [input]  E : ne x 3 edge table; who : the routine that cannot handle one
      %
      % Called at the entry of every routine whose edge handling is written as "ray or segment".
      % Those routines are not wrong -- a closed edge is genuinely new geometry and each will need
      % its own answer -- but they must not treat flag 2 as a segment and index vertex 0. One
      % definition site, so the message and the convention stay in one place.
         if isempty(E), return, end
         if any(E(:,3) == 2)
            error('RatCon:closedEdge', ...
               ['%s cannot yet handle a CLOSED edge (E(:,3) = 2): its edge handling assumes ' ...
                'every edge names two endpoints. A closed conic bounds the bridge face domain ' ...
                '(DECISIONS.md 2026-09-06) and each consumer needs its own answer for one. ' ...
                'Edge(s): %s.'], who, mat2str(find(E(:,3) == 2).'));
         end
      end
   end
   methods
       % NOTE ON THE ERROR IDENTIFIERS BELOW. `set.den` and `set.Ec` still raise `RatPar:...`
       % and the edge-flag statics still raise `RatCon:...`, although all of them now live on
       % AlgCon. That is deliberate and it is the same reason twice: an error identifier is a
       % PUBLIC, CATCHABLE name, the conditions they report are unchanged, and renaming one would
       % break every caller and test that pins it for no gain (`RatParTest` pins two of them).
       % The identifiers that are NEW here -- `AlgCon:FPMustBeEmpty`, `AlgCon:FBMustBeEmpty` --
       % name this class, because nothing pinned them before.
       function obj = set.den(obj, val)
       % Enforce the Qua axis from the single place `den` is defined: a Qua IS a Rat whose
       % denominator is 1, so a QuaPar/QuaPol/QuaParCPLQ must never carry a non-trivial one.
       % Reading the trait (rather than the class name) is what lets one validator serve the whole
       % lattice, including types added later.
            if ~isempty(val) && isa(obj,'Qua')
                trivial = repmat([0 0 1], size(val,1), 1);
                if any(val ~= trivial, 'all')
                    error('RatPar:denMustBeTrivial', ...
                        ['a %s is quadratic (isa Qua), so its denominator is identically 1: ' ...
                         'den must be [0 0 1] on every face.'], class(obj));
                end
            end
            obj.den = val;
       end

       function obj = set.Ec(obj, val)
       % Enforce the Pol axis from the single place `Ec` is defined: a Pol IS a Par with no
       % curvature, so a RatPol/QuaPol must never carry a nonzero conic. See set.den.
            if ~isempty(val) && isa(obj,'Pol') && any(val ~= 0, 'all')
                error('RatPar:EcMustBeZero', ...
                    ['a %s is polyhedral (isa Pol), so every edge is straight: ' ...
                     'Ec must be all zero.'], class(obj));
            end
            obj.Ec = val;
       end

       function obj = set.EcA(obj, val)
       % Enforce the Crv axis from the single place `EcA` is defined: only a Crv may carry an edge
       % that is not a conic. Reading the TRAIT rather than the class name is the same argument
       % set.den, set.Ec and set.FP rest on, and it covers types that do not exist yet.
            if ~isempty(val) && ~isa(obj,'Crv')
                error('AlgCon:EcAMustBeEmpty', ...
                    ['a %s has a CONIC subdivision (it is not a Crv), so every edge is a conic ' ...
                     'and none is an arbitrary algebraic curve: EcA must be empty. An edge that ' ...
                     'needs a polynomial belongs on a Crv -- see Crv.m.'], class(obj));
            end
            obj.EcA = val;
       end

       function obj = set.FP(obj, val)
       % Enforce the Alg axis from the single place `FP` is defined: every face of a Rat has a
       % FORMULA, so a rational type must never carry an implicit one. Reading the TRAIT rather
       % than the class name is what lets one validator serve the whole lattice -- the same
       % argument set.den and set.Ec above rest on, and it covers types that do not exist yet.
            if ~isempty(val) && isa(obj,'Rat')
                error('AlgCon:FPMustBeEmpty', ...
                    ['a %s is rational (isa Rat), so every face has a formula and none is cut ' ...
                     'out implicitly: FP must be empty. A face that needs a polynomial belongs ' ...
                     'on an AlgCon.'], class(obj));
            end
            obj.FP = val;
       end

       function obj = set.FB(obj, val)
       % The branch index is meaningless without a polynomial to index into, so it is pinned on
       % the same axis and for the same reason. All-zero is allowed everywhere: that is what a
       % face with a formula carries, and a producer filling FB before FP should not be stopped.
            if ~isempty(val) && any(val ~= 0) && isa(obj,'Rat')
                error('AlgCon:FBMustBeEmpty', ...
                    ['a %s is rational (isa Rat), so no face names a root and there is no branch ' ...
                     'to select: FB must be empty or all zero.'], class(obj));
            end
            obj.FB = val;
       end

       function k = kind(obj)
       % objective: which concrete type this RatPar actually is -- the "type flag" callers switch
       %            on after receiving a statically-typed RatPar from conj/biconj/convEnv/...
       % [output] k : one of 'QuaPol' | 'RatPol' | 'QuaPar' | 'QuaParCPLQ'
       %
       % Derived from the object's real class rather than stored, so it can never disagree with
       % what the object IS -- see this class's header. Prefer kind() over class() in user code:
       % class() also returns the name, but kind() is the documented, stable contract (class names
       % are an implementation detail that a future refactor may change).
            k = class(obj);
       end

       function tf = isMeshed(obj)
       % objective: whether this object carries an explicit V/E/F mesh (true for QuaPol, RatPol
       %            and QuaPar) or is still held in cPLQ's symbolic form with the mesh not yet
       %            reconstructed (false for QuaParCPLQ, whose inherited mesh properties are empty
       %            precisely because reconstructing V/E/Ec/F/P from it is still open -- see
       %            QuaParCPLQ.m and DESIGN.md II.5.1).
       % [output] tf : logical
            tf = ~isempty(obj.f);
       end
   end

   methods (Access = protected)
       function obj = setPinnedDefaults(obj)
       % Fill in whichever axis-varying property this type pins, so that EVERY constructed object
       % is fully formed with respect to the lattice. A QuaPol really is a RatPol with unit
       % denominators and a QuaPar with zero conics, so it should carry them: code that legitimately
       % reads `p.den` after an isa(p,'RatPol') test (e.g. conjPieceCPLQ's rational-piece guard)
       % must not trip over an empty matrix just because the object happens to be the pinned
       % specialization.
       %
       % Called by the LEAF constructor only, after it has written the mesh -- see this file's
       % CONSTRUCTOR PROTOCOL note. Leaves whatever is already set untouched, so a RatPol keeps its
       % genuine denominators and a QuaPar keeps its genuine conics.
            if isempty(obj.den), obj.den = repmat([0 0 1], size(obj.f,1), 1); end
            if isempty(obj.Ec),  obj.Ec  = zeros(size(obj.E,1), 6); end
       end
   end
end

% ==================================================================================================
function e = emptyEcA()
% The empty ALGEBRAIC edge table, with its fields declared so that every EcA -- empty or not -- has
% the same shape. A default of `struct([])` would have no fields, and a caller appending to it would
% silently create whatever it wrote first.
    e = struct('FP', {}, 'FB', {}, 'num', {}, 'den', {});
end
