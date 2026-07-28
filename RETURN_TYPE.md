# The `conj` return-type story

## The problem

`conj` currently returns a **different class depending on the shape of its input**:

| Input | `conj` returns | Where |
|---|---|---|
| Full-domain strictly convex quadratic | `QuaPoly` | `conjCPLQ.m` Case A |
| Single bounded triangle | `QuaPar` | Case B (`conjSingleTriangle`) |
| General bounded domain | `QuaParCPLQ` | Case C (symbolic wrapper) |

A caller cannot write `g = f.conj(); h = g.add(...)` without knowing which of three classes came
back. Worse, `QuaParCPLQ` is not a peer of the other two — it is a thin wrapper around the
vendored cPLQ `functionNDomain`, carrying a symbolic representation rather than a `V/E/Ec/f/F`
mesh. It was given the *operator surface* of a `QuaPar` (`conj`/`add`/`scalarMul`/`addQuadratic`/
`addScaledEnergy`/`eval`) precisely so that `infConv`/`moreau`/`proxAverage`/`QuaPar.biconj`
compose with it unchanged — but that is duck typing, not a type guarantee.

This is the top item blocking a general release: it is the one gap a downstream project hits on
its very first line of code.

---

## Constraint that removes two of the naive options

`QuaPoly` is a **special case of `RatPol`** (a `RatPol` with denominator `[0 0 1]`). So the
candidate return types are not three independent classes — really only **`RatPol`** and
**`QuaPar`** are distinct, and both are special cases of **`RatPar`** (rational cubic/linear on a
parabolic subdivision), the parent `DESIGN.md` II.3 already proposes but which was never built.

---

## Options

### Option A — status quo: return whichever class fits

- **Pro:** no work; each result is the tightest type that can hold it.
- **Con:** callers must dispatch on class. Operators must be implemented N×N across classes.
  `QuaParCPLQ` remains a duck-typed outlier. **Not viable for a public API.**

### Option B — always return `QuaPar`

Widen every result to `QuaPar`, converting `QuaPoly` results up (`toQuaPar.m` already does this).

- **Pro:** single return type; small change; `toQuaPar` exists.
- **Con:** **wrong for the biconjugate direction.** `f**` of a nonconvex `QuaPoly` is a `RatPol`
  (a genuinely rational convex envelope, [COAP] A.3/A.4) — `QuaPar` has no denominator and cannot
  represent it. Also forces a full-domain quadratic to carry an empty mesh. Fails on the type
  cycle `DESIGN.md` II.2 actually needs.

### Option C — build `RatPar` as a real parent; return `RatPar` with a type flag *(recommended)*

Make `RatPol` and `QuaPar` (and hence `QuaPoly`) inherit from `RatPar`. Every operator returns a
`RatPar` carrying a **type flag** recording what it actually is — `QuaPoly` | `RatPol` | `QuaPar` —
so callers get one static type while the tight structure stays available and cheaply checkable.

- **Pro:**
  - One return type; downstream code needs no class dispatch.
  - Matches the type cycle in `DESIGN.md` II.2 (`conj: QuaPoly → QuaPar`, `convEnv: QuaPoly →
    RatPol`) without lying about either.
  - Operators can be written once on `RatPar` and specialized where a subtype allows something
    faster, instead of N×N pairings.
  - Gives `QuaParCPLQ` a principled home: it becomes a `RatPar` subtype (or a flagged
    representation) rather than a duck-typed lookalike.
  - Preserves exactness — the flag records the true structure, so nothing is widened lossily.
- **Con:** the largest change of the three. Touches all three existing `classdef`s, their
  constructors, and every operator. Needs a decision on flag-vs-`isa` (below).

### Option D — return `RatPar` with no flag, inferring the type on demand

Same as C but derive "is this really a `QuaPoly`?" by inspecting the data (all denominators
`[0 0 1]`, all `Ec` rows zero, …).

- **Pro:** no flag to keep in sync — cannot go stale.
- **Con:** repeated O(n) scans; and inference cannot distinguish *"this is structurally a QuaPoly"*
  from *"this is a RatPol that happens to be polynomial right now"*, which matters for choosing an
  algorithm. A flag set at construction states intent; inference only observes the current values.

---

## Outcome (implemented 2026-07-27)

**Option C**, refined during design review into a **two-axis product lattice**. A piecewise
function is built from a function type (`Rat` ⊃ `Qua`) and a subdivision type (`Par` ⊃ `Pol`); the
four named classes are the cells of that grid and inheritance is the grid's partial order:

```
      RatPar < Rat & Par            (abstract)
      /    \
 RatPol    QuaPar                   < RatPar & Pol   /   < RatPar & Qua
      \    /
      QuaPol                        < RatPol & QuaPar
```

Four MATLAB-specific facts settled the mechanics, each verified by probe before committing:

| Probe | Result | Consequence |
|---|---|---|
| Diamond with disjoint properties | works, merged once | the lattice is expressible |
| Same property defined in two superclasses | **fatal**, `conflictingSuperClassProperty` | — |
| Child redefining it | **fatal**, `RedefinedProperty` | **unresolvable** — so pinning cannot be done by overriding |
| Same *method* in two superclasses, child overrides | works | method-level conflicts *are* resolvable |
| Property-less abstract markers under MI | works, `isa` true on both axes | traits are viable |
| Base constructor in a diamond | **runs twice, and clobbers** | forced the leaf-only constructor protocol |

So `den` and `Ec` live on `RatPar` alone, and each type's pinned value (`den ≡ 1` for a `Qua`,
`Ec ≡ 0` for a `Pol`) is enforced by `set.` validators that read the object's own **traits** — one
definition site covering the whole lattice, and no type can be made to lie about itself.

`kind()` remains a method over `class(obj)` rather than a stored field, so it cannot drift.

**Degree is deliberately not a type.** `Rat` is really *cubic numerator over linear denominator*, so
one might add `Cub` between `Rat` and `Qua` — but nothing in the toolbox dispatches on degree (every
consumer either rejects cubics or ignores the distinction; only `isConvex` accepts them). A type
should track what changes the *algorithm*, not what changes the data's shape. Promoting degree would
make the lattice a 2×2×2 product — eight combinations — for no dispatch benefit. Revisit only if an
operator ever *produces* cubics rather than merely tolerating them.

### Original recommendation, for the record

**Option C**, matching your reading: everything returns a `RatPar` carrying a type flag indicating
whether it is in fact a `RatPol`, `QuaPoly`, or `QuaPar`.

Two design points worth settling before implementation:

1. **Flag vs. `isa`.** If `RatPol`/`QuaPar`/`QuaPoly` become genuine subclasses, MATLAB's own
   `isa`/`class` already *is* the flag, and a separate field risks disagreeing with the actual
   class. Suggested rule: **use real subclasses as the source of truth**, and expose a
   `kind()` method returning `'QuaPoly'|'RatPol'|'QuaPar'` for callers who want a plain string to
   switch on. That gives the ergonomics of a flag with no possibility of drift.
   The alternative — one concrete `RatPar` class with a stored flag and no subclasses — is simpler
   to implement but gives up MATLAB's dispatch, which the current operator code relies on heavily
   (e.g. `maxQuaPar`/`infConv` compose with `QuaParCPLQ` *only* because MATLAB dispatches on the
   operand's actual class).
2. **Value-class semantics.** All three classes are currently value classes; `RatPar` must stay
   one, or every operator's copy semantics change.

**Sequencing note.** This is a breaking API change. It is much cheaper to do it *before* tagging
0.1 than after — a 0.1 that pins the current shape-dependent return type would have to be broken
by 0.2. If 0.1 is meant to be usable by other projects, this should land first.

---

## Scope guard

Building `RatPar` does **not** mean implementing every possible `RatPar`. Per `SUPPORT_MATRIX.md`
§0, CCA2's goal is `QuaPoly` conjugate/biconjugate; the `RatPol`s and `QuaPar`s that arise are very
special cases (e.g. a parabolic edge only ever occurs surrounded by two parallel rays, and
hyperbolic edges never arise at all). `RatPar` is needed here as a **common type**, not as a
general rational-cubic-on-parabolic engine.
