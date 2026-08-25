classdef (Abstract) RatCon < Rat & Conic
   % RatCon: the common PARENT of every function type CCA2's operators return, and the ONE
   % place the mesh data is declared.
   %
   % PURPOSE -- one static return type. Before this class existed, `conj` returned a DIFFERENT
   % class depending on the shape of its input (QuaPol for a full-domain strictly convex
   % quadratic, QuaPar for a single bounded triangle, QuaParCPLQ for a general bounded domain), so
   % no caller could write `g = f.conj()` without dispatching on class first. Every one of those
   % is now a RatPar, so `conj`/`biconj`/`convEnv`/... are statically typed: they return a RatPar,
   % and `kind(g)` says what it actually is. See RETURN_TYPE.md.
   %
   % TWO AXES, NOT ONE LIST. A piecewise function is built from two independent choices: the
   % FUNCTION type (rational, or its polynomial specialization) and the SUBDIVISION type (parabolic,
   % or its polyhedral specialization). Those axes are modelled as property-less abstract markers --
   % Rat with Qua < Rat, and Par with Pol < Par -- so either axis can be queried on its own:
   %
   %       isa(f,'Qua')   the denominator is pinned to 1
   %       isa(f,'Pol')   every edge conic is pinned to 0
   %
   % The four named types are the four cells of that 2x2 grid, and the inheritance is the grid's own
   % partial order -- a genuine diamond, because that is what a product lattice IS:
   %
   %                 RatCon  (rational on conic)        < Rat & Conic   <- this class
   %                    |
   %                 RatPar  (rational on parabolic)     < RatCon & Par
   %                 /    \
   %      (rational on    (quadratic on
   %       polyhedral)     parabolic)
   %          RatPol        QuaPar                       < RatPar & Pol / RatPar & Qua
   %                 \    /
   %                 QuaPol  (quadratic on polyhedral)   < RatPol & QuaPar
   %
   %       QuaCon      quadratic on a CONIC subdivision              < RatCon & Qua
   %       QuaParCPLQ  a QuaPar still in cPLQ's own symbolic form    < RatPar & Qua
   %
   % The subdivision axis has THREE levels, Pol < Par < Conic, because a conjugate's edges are only
   % parabolic when the two pieces compared are adjacent (CONJ_FIELD_PROOF.md Theorem 3); Step 3
   % compares non-adjacent pieces too, and those give genuine ellipses. See Conic.m. Con was added
   % ABOVE the previous top, so no existing type's behaviour changed: RatPar and everything under
   % it still is-a Par and still cannot hold a non-parabolic conic.
   %
   % QuaPol IS-A RatPol (a polynomial is a rational with unit denominator) and IS-A QuaPar (a
   % polyhedral subdivision is a parabolic one with no curvature), so both edges of the diamond are
   % mathematically honest and `isa` reports the full lattice.
   %
   % DATA LIVES HERE, ONCE. `den` and `Ec` are declared in this class alone -- NOT overridden in the
   % children -- because MATLAB makes a property defined in two superclasses fatal AND unresolvable
   % (MATLAB:class:conflictingSuperClassProperty; a child cannot redefine it either:
   % MATLAB:class:RedefinedProperty). Had RatPol declared `den` and QuaPar declared its own pinned
   % `den`, QuaPol could not exist. Instead each type's pinned values are enforced by the
   % `set.den`/`set.Ec` validators below, which read the object's own TRAITS -- so one definition
   % site serves the whole lattice and no type can be made to lie about itself.
   %
   % QuaParCPLQ deliberately does NOT inherit from QuaPar: it carries the same mathematical type
   % (marked Qua, and Par via this class) but not the mesh, so making `isa(x,'QuaPar')` true for it
   % would invite callers to read V/E/F that it does not have. `isMeshed()` is the intended test.
   %
   % NO STORED TYPE FLAG, BY DESIGN. `kind()` is a METHOD returning the object's actual class name,
   % not a property set at construction: a stored flag can drift out of sync with the real class
   % (after a copy, a subclass change, or a constructor path that forgets to set it), and a flag
   % that disagrees with the object is worse than no flag. Deriving it from `class(obj)` makes that
   % impossible. It also keeps MATLAB's own dispatch as the single source of truth -- which the
   % operator code already depends on heavily (e.g. maxQuaPar/infConv/moreau/proxAverage compose
   % with QuaParCPLQ only because MATLAB dispatches on the operand's actual class).
   %
   % VALUE SEMANTICS. RatPar is a VALUE class (not a handle), matching all four children: every
   % operator returns a new object and never mutates its argument. Do not change this without
   % auditing every operator, which all rely on copy-on-assign.
   %
   % SCOPE. Introducing RatPar does NOT mean implementing every possible RatPar. CCA2's goal is
   % QuaPol conjugate/biconjugate; the RatPols and QuaPars that actually arise from it are very
   % special cases (e.g. a parabolic edge only ever occurs surrounded by two parallel rays, and
   % hyperbolic edges never arise at all -- see [COAP]/[JOGO] and SUPPORT_MATRIX.md section 0).
   % RatPar exists here as a COMMON TYPE, not as a general rational-cubic-on-parabolic engine.

   properties
        % Everything the lattice shares, declared here ONCE -- the nine mesh properties AND the two
        % axis-varying ones (`den`, `Ec`). No child declares any property of its own; see the
        % "DATA LIVES HERE, ONCE" note above for why that is forced rather than merely tidy.
        nv {mustBeInteger,mustBeNonnegative}%number of vertices
        ne {mustBeInteger,mustBeNonnegative}%number of edges
        nf {mustBeInteger,mustBeNonnegative}%number of faces
        V (:,2){mustBeNumeric} % nv x 2 matrix storing unique vertices
        E (:,3){mustBeInteger,mustBeNonnegative} % ne x 3 matrix storing edge indices where Edge j is [V(E(j,1)), V(E(j,2))]
            %and E(j,3) = 1 for a segment and 0 for a ray.
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

   methods
       % NOTE ON THE ERROR IDENTIFIERS BELOW. They still read `RatPar:...` although the
       % validators now live on RatCon. That is deliberate: an error identifier is a public,
       % catchable name, and the two conditions it reports are unchanged. Renaming it would break
       % every caller and test that pins it for no gain.
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
