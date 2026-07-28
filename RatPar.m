classdef (Abstract) RatPar
   % RatPar: the common PARENT of every function type CCA2's operators return.
   %
   % PURPOSE -- one static return type. Before this class existed, `conj` returned a DIFFERENT
   % class depending on the shape of its input (QuaPoly for a full-domain strictly convex
   % quadratic, QuaPar for a single bounded triangle, QuaParCPLQ for a general bounded domain), so
   % no caller could write `g = f.conj()` without dispatching on class first. Every one of those
   % is now a RatPar, so `conj`/`biconj`/`convEnv`/... are statically typed: they return a RatPar,
   % and `kind(g)` says what it actually is. See RETURN_TYPE.md.
   %
   % WHY THIS PARENT, AND WHY THESE CHILDREN. The name follows DESIGN.md II.3: RatPar = RATional
   % (quadratic/linear) on a PARabolic subdivision, the most general type in the family. The
   % specializations are:
   %      RatPar (rational on parabolic)                <- this class, abstract
   %        |-- RatPol   rational on POLYhedral         (adds `den`: per-face linear denominator)
   %        |-- QuaPar   QUAdratic on PARabolic         (adds `Ec` : per-edge conic)
   %        |-- QuaPoly  QUAdratic on POLYhedral        (adds neither)
   %        |-- QuaParCPLQ  a QuaPar still in cPLQ's own symbolic form (mesh not reconstructed)
   %
   % NOTE that QuaPoly is mathematically a special case of RatPol (denominator identically 1), and
   % RatPol/QuaPar are both special cases of RatPar -- but the MATLAB hierarchy deliberately makes
   % all four SIBLINGS under RatPar rather than chaining QuaPoly < RatPol. Chaining would make
   % `isa(aQuaPoly,'RatPol')` true, and code that (correctly) reads `p.den` after such a test would
   % then fail on a QuaPoly, which has no denominator to read -- e.g. conjPieceCPLQ's rational-piece
   % guard. Siblings keep every existing `isa` test meaning exactly what it meant before.
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
   % QuaPoly conjugate/biconjugate; the RatPols and QuaPars that actually arise from it are very
   % special cases (e.g. a parabolic edge only ever occurs surrounded by two parallel rays, and
   % hyperbolic edges never arise at all -- see [COAP]/[JOGO] and SUPPORT_MATRIX.md section 0).
   % RatPar exists here as a COMMON TYPE, not as a general rational-cubic-on-parabolic engine.

   properties
        % The mesh every child shares, declared here ONCE. These nine were previously repeated
        % verbatim (identical names, sizes and validators) in QuaPoly.m, RatPol.m and QuaPar.m;
        % they are reproduced unchanged, so behaviour is identical. Children add only what is
        % genuinely their own: RatPol adds `den`, QuaPar adds `Ec`, QuaPoly adds nothing.
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
   % Consequence for callers: `QuaPoly()`, `QuaPar()`, `RatPol()` with no arguments now return a
   % blank object instead of raising an error. That is also what MATLAB itself needs for
   % `ClassName.empty`, object arrays and `load`, so it is a gain rather than a concession.
   % ---------------------------------------------------------------------------------------

   methods
       function k = kind(obj)
       % objective: which concrete type this RatPar actually is -- the "type flag" callers switch
       %            on after receiving a statically-typed RatPar from conj/biconj/convEnv/...
       % [output] k : one of 'QuaPoly' | 'RatPol' | 'QuaPar' | 'QuaParCPLQ'
       %
       % Derived from the object's real class rather than stored, so it can never disagree with
       % what the object IS -- see this class's header. Prefer kind() over class() in user code:
       % class() also returns the name, but kind() is the documented, stable contract (class names
       % are an implementation detail that a future refactor may change).
            k = class(obj);
       end

       function tf = isMeshed(obj)
       % objective: whether this object carries an explicit V/E/F mesh (true for QuaPoly, RatPol
       %            and QuaPar) or is still held in cPLQ's symbolic form with the mesh not yet
       %            reconstructed (false for QuaParCPLQ, whose inherited mesh properties are empty
       %            precisely because reconstructing V/E/Ec/F/P from it is still open -- see
       %            QuaParCPLQ.m and DESIGN.md II.5.1).
       % [output] tf : logical
            tf = ~isempty(obj.f);
       end
   end
end
