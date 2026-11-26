# Referee report

**Paper:** *Exact Conjugates and Convex Envelopes of Bivariate Piecewise Quadratic Functions: Algorithms, Proofs, and a Verified Implementation*  
**Venue:** Mathematical Programming Computation  
**Manuscript:** `manuscript/main.pdf`  
**Recommendation: Major revision**

## Summary of the contribution

The manuscript develops exact algorithms for the Fenchel conjugate and closed convex envelope of rational bivariate piecewise linear–quadratic functions. It characterizes possible algebraic forms, challenges a published parabolic-subdivision claim, and reports three implementations with computational checks. The structural results and software could be valuable to MPC readers working in computational optimization.

## Strengths

- The paper has a clear and interesting mathematical narrative: the answer type is richer than the previously claimed parabolic subdivision, and this affects the data structures an exact implementation must support (Sections 3.1–3.4).
- The distinction between rational face functions, conic boundaries, algebraic vertices, and algebraic envelope faces is concrete and useful. The degree-four vertex example and the elliptic-edge counterexample give the claims memorable witnesses.
- The authors are unusually candid about scope and limitations, including named refusals, thin conjugates, slow cases, uncompleted ranges, and the absence of a citable release (Sections 8.4 and 9).
- The implementation account distinguishes exact decisions from sampling and reports multiple independent checking approaches, cross-language parity, and reproducible benchmark inputs. The shared-machine caveat is stated rather than hidden.

## Major concerns / weaknesses

1. **The central mathematical claims need proofs in the paper.** Section 3 says “The full proofs are in the repository.” Theorem 2 receives only a brief explanation; Theorem 3 has a useful outline; Theorem 5 gives a mechanism but not a complete derivation; and the counterexample’s key algebraic and geometric assertions are asserted with computational certificates. For a paper whose title and abstract emphasize algorithms and proofs, a reader should not need repository materials to assess the main results. Include self-contained proofs or appendices for the degree bounds, the biconjugate face-type theorem and its hypotheses, and the counterexample (including the switching region and nondegeneracy argument). State precisely which steps are symbolic derivations and which are certified computations.

2. **The contribution and novelty need sharper positioning.** The manuscript presents three contributions—new algorithms, new face-type results, and a verification harness—but does not yet make their boundaries easy to assess against the two closest predecessor papers. In particular, Section 2.2 says the earlier algorithms remain valid on inputs where their claim holds, while Section 3 argues that their stated general face classification is false. Explain exactly what remains usable from those algorithms, what the new algorithms do differently, and which mathematical results are new independently of the implementation. A concise comparison table by operator, input assumptions, output representation, complexity claim, and proof status would help.

3. **The completeness and exactness statements need a consistent scope.** Exact arithmetic for decisions is discussed alongside claims such as “complete,” “exact,” and “no wrong value,” while the evidence is bounded by tested classes, finite corpora, refusals, timeouts, and sampled points. Separate theorem-level guarantees from implementation behavior and validation. State the supported input conditions and whether termination and completeness are proved; clarify “four oracles that share no code” given the stated absence of an independent external comparison.

4. **The computational evidence is extensive but its independence is overstated.** Section 9.2 marks “independent verification of results” as met, yet Section 8.4 says there is no independent external comparison and Section 10 lists the external oracle set as open. Distinguish internal cross-checking from external validation in the compliance table, and provide a citable, versioned artifact with reproduction instructions. The shared-machine timings are appropriately caveated; retain them as indicative observations until controlled measurements are available.

5. **The manuscript is too repository-centered for a standalone article.** Numerous claims point to tests and internal records rather than giving readers the mathematical evidence. Work through the three-piece counterexample from input to switching curve, and illustrate the biconjugate pipeline with a reproducible example. Move implementation history and low-level test details to an appendix or supplement.

## Minor concerns

- Figure 1’s face labels are difficult to read at the printed size; enlarge the detailed panels or provide a companion table/key. Figure 2’s legend is also very small. The captions explain the purpose but the figures should remain legible in print.
- The abstract is dense; a shorter statement of the problem, main result, and evidence would improve accessibility.
- The introduction first says continuity is assumed and then introduces lower-semicontinuous interpretation for discontinuous input. Reconcile these statements by stating the main theorem’s assumptions there and separating the extended implementation behavior.
- Section 8.4’s 25-test result is implementation validation, not a mathematical certificate; state the corresponding proposition or finite test domain.
- Section 9.2’s venue-compliance table should distinguish “code and data available” from “citable artifact”; currently the former is marked met while the lack of a citable release is separately listed as a submission-blocking gap.

## Assessment of fit for Mathematical Programming Computation

The topic fits MPC well: it combines exact optimization algorithms, mathematical structure, and substantial computational artifacts. The paper could make a worthwhile contribution to the journal if the proofs are brought into the article, the algorithmic guarantees are separated cleanly from empirical evidence, and the artifact and verification claims are made consistent. In its present form, the mathematical contribution is promising but not yet assessable from the manuscript alone, and the computational reproducibility case is described more strongly in the compliance table than the stated evidence supports.

## Score sheet

| Criterion | Rating (1-5) | Notes |
|---|---:|---|
| Originality / novelty | 4 | The face-type results and counterexamples appear potentially significant; comparison with prior work needs sharpening. |
| Technical soundness / correctness | 3 | The statements are plausible and some arguments are sketched, but several central proofs are deferred to the repository. |
| Significance / impact | 4 | Exact operators for nonconvex bivariate PLQ functions have a credible optimization use case. |
| Clarity & presentation | 3 | Strong organization and candid limitations, offset by density, repository-dependent exposition, and small figure labels. |
| Fit for MPC | 4 | Good subject fit, contingent on stronger mathematical exposition and a clearer computational-evidence account. |
| Reproducibility & artifact quality | 3 | Multiple implementations and data are described, but external validation and a citable release remain absent. |

## Overall

I recommend major revision. The paper has a compelling core result and a serious implementation effort, but the central theorems and algorithmic guarantees are not sufficiently established in the manuscript itself. A revised version should supply the proofs, make the novelty and scope precise, and align its reproducibility claims with the acknowledged absence of external validation and a citable artifact.

## Priority list for the authors

1. Put complete proofs or detailed appendices for the principal structural theorems and counterexamples in the submission.
2. State theorem-level guarantees separately from empirical results, and reconcile the “independent verification” compliance claim with the absence of an external comparison.
3. Position the algorithms against the predecessor work with a precise comparison of assumptions, outputs, and guarantees.
4. Make the three-piece example and one envelope computation self-contained, while moving repository-specific history and test bookkeeping out of the main narrative.
5. Improve figure legibility and provide a citable, versioned artifact before making a full reproducibility claim.
