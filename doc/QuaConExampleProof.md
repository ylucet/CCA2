# Why the \(I_3\,|\,I_5\) arc really is a hyperbola

_Written 2026-09-10. A self-contained proof for the reader of `doc/QuaConExample.md` who is not
convinced that a **non-degenerate** hyperbolic arc can occur in the conjugate of a bivariate PLQ
function, because every example she built came out as a **degenerate** one — a pair of straight
lines._

**The short version.** Her examples are not wrong and the document is not wrong. She has been
comparing pieces that share an edge, and for such a pair the conic is *forced* to degenerate into
two lines — that is Proposition 5 below, and it is the reason a two-piece example can never
produce a curved arc. The document's arc comes from two pieces that share **no** edge, and for
those nothing is forced. Whether a given pair gives a genuine hyperbola, a genuine ellipse, or two
lines is decided by **two determinants of the input data**, computable by hand before any conjugate
is calculated (Theorem 6). For the pair in question they are \(-265\) and \(144\), and that is the
whole proof.

Everything here is exact rational arithmetic. Section 9 lists what was recomputed, and with which
independent tool.

---

## 1. Setting and notation

Let \(f:\mathbb{R}^2\to\mathbb{R}\cup\{+\infty\}\) be piecewise linear-quadratic in the sense of
Rockafellar–Wets Definition 10.20: \(\operatorname{dom} f\) is a finite union of polyhedra
\(T_1,\dots,T_N\) and on each of them \(f\) agrees with one quadratic expression

\[
q_k(x) \;=\; \tfrac12\,\langle x, Q_k x\rangle \;+\; \langle \beta_k, x\rangle \;+\; \gamma_k ,
\]

with \(Q_k\) symmetric. Throughout we assume every \(Q_k\) is **positive definite**, which is the
case in `doc/QuaConExample.md` §4.1; \(f\) itself is *not* convex, because the gradient jumps
across the internal edges.

Since \(\operatorname{dom} f\) is the union of the \(T_k\),

\[
f^*(s) \;=\; \sup_{x\in \operatorname{dom} f}\bigl[\langle s,x\rangle - f(x)\bigr]
       \;=\; \max_{1\le k\le N} h_k(s),
\qquad
h_k(s) \;:=\; \sup_{x\in T_k}\bigl[\langle s,x\rangle - q_k(x)\bigr],
\tag{1.1}
\]

so \(f^*\) is a finite max of \(N\) functions, each of which is the conjugate of one quadratic
*restricted to one polyhedron*. Write \(g_k\) for the conjugate of the same quadratic taken over
all of \(\mathbb{R}^2\):

\[
g_k(s) \;:=\; \sup_{x\in\mathbb{R}^2}\bigl[\langle s,x\rangle - q_k(x)\bigr].
\]

**Lemma 1 (the interior cell).**
\(g_k(s)=\tfrac12\langle s-\beta_k,\;Q_k^{-1}(s-\beta_k)\rangle-\gamma_k\), and the unique
maximiser is \(x_k(s)=Q_k^{-1}(s-\beta_k)\).

*Proof.* The objective \(x\mapsto \langle s,x\rangle - q_k(x)\) is strictly concave (because
\(Q_k\succ 0\)), so it has a unique maximiser, characterised by first-order optimality
\(s-Q_kx-\beta_k=0\), i.e. \(x=x_k(s)\). Substituting, and writing \(x=x_k(s)\) for brevity,

\[
\begin{aligned}
\langle s,x\rangle - \tfrac12\langle x,Q_kx\rangle - \langle\beta_k,x\rangle - \gamma_k
 &= \langle s-\beta_k, x\rangle - \tfrac12\langle x,Q_kx\rangle - \gamma_k \\
 &= \langle Q_kx, x\rangle - \tfrac12\langle x,Q_kx\rangle - \gamma_k
    &&\text{(since } s-\beta_k=Q_kx\text{)}\\
 &= \tfrac12\langle x,Q_kx\rangle - \gamma_k \\
 &= \tfrac12\langle s-\beta_k,\;Q_k^{-1}(s-\beta_k)\rangle - \gamma_k .
\end{aligned}
\]
\(\blacksquare\)

Define the **interior cell** of piece \(k\),

\[
C_k \;:=\; \{\, s \;:\; x_k(s)\in\operatorname{int} T_k \,\},
\tag{1.2}
\]

an open convex polyhedron (it is the preimage of \(\operatorname{int} T_k\) under the affine
bijection \(s\mapsto Q_k^{-1}(s-\beta_k)\)). On \(C_k\) the unconstrained maximiser is feasible for
\(T_k\), hence

\[
s\in C_k \quad\Longrightarrow\quad h_k(s)=g_k(s).
\tag{1.3}
\]

These are the faces \(I_1,\dots,I_5\) of `doc/QuaConExample.md` §3.2. If, on some open set, the max
in (1.1) is attained by piece \(i\) on one part and by piece \(j\) on another, and both parts meet
the respective interior cells, then the boundary between them is contained in the level set

\[
Z_{ij} \;:=\; \{\, s \;:\; g_i(s)=g_j(s) \,\},
\tag{1.4}
\]

which is a conic, since \(g_i-g_j\) is a quadratic polynomial in \(s\). **The question "is the arc
curved?" is the question "is \(Z_{ij}\) a non-degenerate conic?"** — a degenerate conic is a union
of lines, so an arc of one is a straight segment.

---

## 2. The two invariants of a conic, and what they decide

Write a quadratic polynomial in \(s=(s_1,s_2)\) as

\[
P(s) \;=\; A\,s_1^2 + 2B\,s_1s_2 + C\,s_2^2 + 2p_1 s_1 + 2p_2 s_2 + R
\tag{2.1}
\]

and attach to it the symmetric \(3\times3\) matrix and its leading \(2\times2\) block

\[
M \;=\; \begin{pmatrix} A & B & p_1 \\ B & C & p_2 \\ p_1 & p_2 & R \end{pmatrix},
\qquad
A_2 \;=\; \begin{pmatrix} A & B \\ B & C \end{pmatrix},
\qquad
p \;=\; (p_1,p_2),
\]

so that \(P(s)=\bigl\langle (s,1),\, M\,(s,1)\bigr\rangle\). What follows are three elementary
facts, proved here rather than cited because the precise form of the second is what the rest of the
note leans on.

**Proposition 2 (centre and canonical form).** If \(\det A_2\neq 0\) then \(P\) has a unique
critical point (its *centre*)

\[
c \;=\; -A_2^{-1}p,
\qquad\text{and}\qquad
P(c+u) \;=\; \langle u, A_2 u\rangle + k,
\qquad
k \;=\; R - \langle p, A_2^{-1}p\rangle \;=\; \frac{\det M}{\det A_2}.
\tag{2.2}
\]

*Proof.* Expanding, \(P(s)=\langle s,A_2s\rangle+2\langle p,s\rangle+R\), so with \(s=c+u\),

\[
P(c+u) \;=\; \langle u,A_2u\rangle \;+\; 2\langle A_2c+p,\,u\rangle
            \;+\; \bigl[\langle c,A_2c\rangle+2\langle p,c\rangle+R\bigr].
\]

The linear term vanishes exactly when \(A_2c+p=0\), i.e. \(c=-A_2^{-1}p\), which is the unique
solution because \(A_2\) is invertible. The constant is then

\[
\langle A_2^{-1}p,\,p\rangle - 2\langle p,\,A_2^{-1}p\rangle + R
 \;=\; R - \langle p, A_2^{-1}p\rangle \;=\; k .
\]

For the second expression for \(k\), use the Schur-complement determinant formula on the block
matrix \(M=\begin{pmatrix}A_2 & p\\ p^{\mathsf T} & R\end{pmatrix}\):

\[
\det M \;=\; \det A_2 \cdot \bigl(R - p^{\mathsf T}A_2^{-1}p\bigr) \;=\; k\,\det A_2 .
\]
\(\blacksquare\)

The **type** of the conic is read off \(\det A_2\). Writing (2.1) in the more familiar
\(a s_1^2 + b s_1 s_2 + c s_2^2 + \dots\) form, \(a=A\), \(b=2B\), \(c=C\), so

\[
b^2-4ac \;=\; 4B^2-4AC \;=\; -4\det A_2 ,
\tag{2.3}
\]

and \(\det A_2<0\) — equivalently \(b^2-4ac>0\), equivalently \(A_2\) **indefinite** — is the
hyperbolic case. Diagonalising \(A_2 = O\operatorname{diag}(\lambda_1,\lambda_2)O^{\mathsf T}\) by
a rotation \(O\), Proposition 2 turns \(P=0\) into

\[
\lambda_1 u^2 + \lambda_2 v^2 + k \;=\; 0,
\qquad \lambda_1\lambda_2 \;=\; \det A_2 .
\tag{2.4}
\]

**Corollary 3 (in the hyperbolic case, \(\det M\) is exactly the "two lines or not" switch).**
Suppose \(\det A_2<0\), so \(\lambda_1<0<\lambda_2\). Then:

* if \(k\neq 0\) (equivalently \(\det M\neq 0\)) the solution set of (2.4) is a **genuine
  hyperbola**: two disjoint unbounded branches, centre \(c\), semi-axes
  \(\sqrt{|k/\lambda_1|}\) and \(\sqrt{|k/\lambda_2|}\), and two distinct real asymptotes through
  \(c\) whose directions are the two real roots \(m\) of the leading form \(A+2Bm+Cm^2=0\);
* if \(k=0\) (equivalently \(\det M=0\)) then (2.4) reads \(\lambda_1u^2=-\lambda_2v^2\), i.e.
  \(v=\pm\sqrt{-\lambda_1/\lambda_2}\,u\): exactly **two distinct real lines through \(c\)** — a
  *degenerate* hyperbola.

*Proof.* Immediate from (2.4): divide by \(-k\) in the first case, which after possibly swapping
the two axes gives \(u^2/a^2-v^2/b^2=1\) with \(a^2=|k/\lambda_1|\), \(b^2=|k/\lambda_2|\); and
factor the difference of two squares in the second. Finally \(k=0\iff\det M=0\) because
\(\det M=k\det A_2\) with \(\det A_2\neq 0\). \(\blacksquare\)

So for the pair we care about, *everything* hinges on the sign of \(\det A_2\) and the vanishing of
\(\det M\). It is worth recording, in general form, the one direction of the factorisation
criterion that we use later.

**Lemma 4 (a product of affine forms has \(\det M=0\)).** If \(P=\ell\,m\) with
\(\ell(s)=a_1s_1+a_2s_2+a_3\) and \(m(s)=b_1s_1+b_2s_2+b_3\) affine, then, writing
\(\alpha=(a_1,a_2,a_3)\) and \(\beta=(b_1,b_2,b_3)\) in \(\mathbb{R}^3\),

\[
M \;=\; \tfrac12\bigl(\alpha\beta^{\mathsf T} + \beta\alpha^{\mathsf T}\bigr),
\qquad\text{hence}\qquad \det M = 0,
\]

and moreover

\[
\det A_2 \;=\; a_1b_1\,a_2b_2 - \Bigl(\frac{a_1b_2+a_2b_1}{2}\Bigr)^{\!2}
         \;=\; -\Bigl(\frac{a_1b_2-a_2b_1}{2}\Bigr)^{\!2} \;\le\; 0 ,
\tag{2.5}
\]

with strict inequality iff \(\ell\) and \(m\) are non-parallel. In words: **a real pair of lines is
always of hyperbolic type — or, if the lines are parallel, parabolic — never elliptic.**

*Proof.* Homogenise: \(\ell(s)=\langle\alpha,(s,1)\rangle\) and \(m(s)=\langle\beta,(s,1)\rangle\),
so \(P(s)=\bigl\langle (s,1),\,\alpha\beta^{\mathsf T}(s,1)\bigr\rangle\); symmetrising the matrix,
which does not change the quadratic form, gives
\(M=\tfrac12(\alpha\beta^{\mathsf T}+\beta\alpha^{\mathsf T})\). That is a sum of two rank-one
matrices, so \(\operatorname{rank} M\le 2<3\) and \(\det M=0\). For (2.5), expand \(\ell m\) to
read off \(A=a_1b_1\), \(C=a_2b_2\), \(2B=a_1b_2+a_2b_1\), and apply the elementary identity
\(\bigl(\frac{p+q}{2}\bigr)^2-pq=\bigl(\frac{p-q}{2}\bigr)^2\) with \(p=a_1b_2\), \(q=a_2b_1\).
\(\blacksquare\)

> **Aside, and it matters for the reader's confusion.** By Lemma 4 the degenerate conics she has
> been getting — two crossing lines — *are* hyperbolas. They are the degenerate members of the
> family. So the question was never "hyperbola or something else"; it was always "degenerate or
> not", i.e. \(\det M=0\) or \(\det M\neq0\).

---

## 3. The conic \(Z_{ij}\) and its primal companion

Fix two pieces \(i\neq j\) and introduce

\[
D \;:=\; Q_i^{-1}-Q_j^{-1},
\qquad
w \;:=\; Q_i^{-1}\beta_i - Q_j^{-1}\beta_j,
\]
\[
\kappa \;:=\; \tfrac12\bigl(\langle\beta_i,Q_i^{-1}\beta_i\rangle
                           - \langle\beta_j,Q_j^{-1}\beta_j\rangle\bigr)
              \;-\; (\gamma_i-\gamma_j).
\]

Expanding the two expressions of Lemma 1,

\[
g_i(s)-g_j(s) \;=\; \tfrac12\langle s, D s\rangle - \langle w, s\rangle + \kappa ,
\]

so, doubling to clear the \(\tfrac12\) and comparing with (2.1),

\[
2\,(g_i-g_j)(s) \;=\; \langle s,Ds\rangle - 2\langle w,s\rangle + 2\kappa
\qquad\text{has matrix}\qquad
M_d \;=\; \begin{pmatrix} D & -w \\ -w^{\mathsf T} & 2\kappa \end{pmatrix}
\tag{3.1}
\]

— the **dual** matrix, i.e. \(A_2=D\), \(p=-w\), \(R=2\kappa\) in the notation of §2. Do the same on
the primal side. With \(\delta:=\beta_i-\beta_j\),

\[
2\,(q_i-q_j)(x) \;=\; \bigl\langle x,(Q_i-Q_j)x\bigr\rangle + 2\langle\delta,x\rangle
                      + 2(\gamma_i-\gamma_j)
\]
\[
\text{has matrix}\qquad
M_p \;=\; \begin{pmatrix} Q_i-Q_j & \delta \\
                          \delta^{\mathsf T} & 2(\gamma_i-\gamma_j) \end{pmatrix}
\tag{3.2}
\]

— the **primal** matrix.

Both matrices are built from the same twelve numbers. The next two lemmas say they carry the same
information.

---

## 4. Lemma A — the type transfers

**Lemma A.** \(\det D \cdot \det Q_i \cdot \det Q_j = \det(Q_i-Q_j)\), and consequently, since
\(\det Q_i\det Q_j>0\),

\[
\operatorname{sign}\bigl(\det D\bigr) \;=\; \operatorname{sign}\bigl(\det(Q_i-Q_j)\bigr).
\]

*Proof.* Factor the difference of inverses:

\[
D \;=\; Q_i^{-1}-Q_j^{-1} \;=\; Q_i^{-1}\,(Q_j-Q_i)\,Q_j^{-1},
\]

which is checked by multiplying out: \(Q_i^{-1}(Q_j-Q_i)Q_j^{-1}=Q_i^{-1}-Q_j^{-1}\). Taking
determinants, \(\det D=\det(Q_j-Q_i)/(\det Q_i\det Q_j)\). Finally, in dimension \(2\),
\(\det(-N)=(-1)^2\det N=\det N\), so \(\det(Q_j-Q_i)=\det(Q_i-Q_j)\). \(\blacksquare\)

Combining Lemma A with (2.3) applied to (3.1) gives a formula worth displaying, because it lets one
compute the discriminant of the arc from the two input Hessians and nothing else. For the
*unscaled* difference \(g_i-g_j\) — whose matrix is \(M_d/2\), so \(A_2=D/2\) —

\[
b^2-4ac \;=\; -4\det\!\Bigl(\frac{D}{2}\Bigr) \;=\; -\det D
        \;=\; -\,\frac{\det(Q_i-Q_j)}{\det Q_i\,\det Q_j}.
\tag{4.1}
\]

**Interpretation.** \(Q_i-Q_j\) indefinite means neither piece is uniformly more curved than the
other — neither \(Q_i\succeq Q_j\) nor \(Q_j\succeq Q_i\) in the Loewner order. *That* is the
geometric content of "the arc is hyperbolic rather than elliptic", and it is a statement purely
about the input.

---

## 5. Lemma B — the degeneracy transfers

**Lemma B.** \(\det M_d \cdot \det Q_i \cdot \det Q_j = -\det M_p\). In particular
\(\det M_d=0 \iff \det M_p=0\).

*Proof.* Both sides are polynomials in the twelve entries
\((Q_i,Q_j,\beta_i,\beta_j,\gamma_i,\gamma_j)\) — on the left, after clearing the inverses through
the factor \(\det Q_i\det Q_j\) — so it suffices to prove the identity on the dense open set where
in addition \(S:=Q_j-Q_i\) is invertible; it then extends to all data by continuity, a polynomial
vanishing on a dense set being identically zero.

Abbreviate \(A:=Q_i\), \(B:=Q_j\), \(a:=\beta_i\), \(b:=\beta_j\), \(c:=\gamma_i-\gamma_j\),
\(\delta:=a-b\), \(S:=B-A\). Apply Proposition 2's Schur formula to each matrix:

\[
\det M_d \;=\; \det D\cdot\bigl(2\kappa - \langle w, D^{-1}w\rangle\bigr),
\tag{5.1}
\]
\[
\det M_p \;=\; \det(A-B)\cdot\bigl(2c - \langle \delta, (A-B)^{-1}\delta\rangle\bigr).
\tag{5.2}
\]

By Lemma A, \(\det D=\det(A-B)/(\det A\det B)\). So the claim
\(\det M_d\,\det A\,\det B=-\det M_p\) reduces, after cancelling \(\det(A-B)\neq 0\), to the scalar
identity

\[
2\kappa - \langle w, D^{-1}w\rangle
 \;=\; -\bigl(2c - \langle\delta,(A-B)^{-1}\delta\rangle\bigr),
\]

i.e., since \(2\kappa=\langle a,A^{-1}a\rangle-\langle b,B^{-1}b\rangle-2c\) and
\((A-B)^{-1}=-S^{-1}\), to

\[
\langle a,A^{-1}a\rangle - \langle b,B^{-1}b\rangle - \langle w, D^{-1}w\rangle
 \;=\; -\,\langle \delta, S^{-1}\delta\rangle .
\tag{5.3}
\]

Now compute. From \(D=A^{-1}SB^{-1}\) (Lemma A's factorisation) we get \(D^{-1}=BS^{-1}A\). Put

\[
u := A^{-1}a, \qquad v := B^{-1}b, \qquad z := u-v,
\]

so that \(a=Au\), \(b=Bv\), \(w=u-v=z\) and \(\delta=Au-Bv\). In these variables (5.3) reads

\[
\langle u,Au\rangle - \langle v,Bv\rangle - \langle z, BS^{-1}Az\rangle
 + \bigl\langle Au-Bv,\; S^{-1}(Au-Bv)\bigr\rangle \;=\; 0 .
\tag{5.4}
\]

Handle the third and fourth terms separately, using \(A=B-S\) in the third and
\(Au-Bv=Bz-Su\) in the fourth. (Check the latter:
\(Bz-Su=Bu-Bv-Su=(B-S)u-Bv=Au-Bv\).)

**Third term.** \(BS^{-1}A=BS^{-1}(B-S)=BS^{-1}B-B\), hence

\[
\langle z,BS^{-1}Az\rangle \;=\; \langle z,BS^{-1}Bz\rangle-\langle z,Bz\rangle .
\]

**Fourth term.** With \(S^{-1}\) symmetric,

\[
\bigl\langle Bz-Su,\,S^{-1}(Bz-Su)\bigr\rangle
 \;=\; \langle z,BS^{-1}Bz\rangle - 2\langle Bz,\,S^{-1}Su\rangle
        + \langle Su,\,S^{-1}Su\rangle
 \;=\; \langle z,BS^{-1}Bz\rangle - 2\langle Bz,u\rangle + \langle Su,u\rangle .
\]

Substituting both into (5.4), the \(\langle z,BS^{-1}Bz\rangle\) terms cancel and we are left with

\[
\langle u,Au\rangle - \langle v,Bv\rangle + \langle z,Bz\rangle
 - 2\langle Bz,u\rangle + \langle u,Su\rangle .
\]

Use \(A+S=B\) on the first and last terms: \(\langle u,Au\rangle+\langle u,Su\rangle
=\langle u,Bu\rangle\). Expand
\(\langle z,Bz\rangle=\langle u,Bu\rangle-2\langle u,Bv\rangle+\langle v,Bv\rangle\) and
\(-2\langle Bz,u\rangle=-2\langle u,Bu\rangle+2\langle u,Bv\rangle\). The total is

\[
\langle u,Bu\rangle - \langle v,Bv\rangle
 + \langle u,Bu\rangle - 2\langle u,Bv\rangle + \langle v,Bv\rangle
 - 2\langle u,Bu\rangle + 2\langle u,Bv\rangle \;=\; 0 ,
\]

which is (5.4). \(\blacksquare\)

Lemmas A and B together say something clean enough to remember:

> **The Legendre transform preserves the projective type of the difference conic.** For two
> strictly convex quadratics, the primal conic \(\{q_i=q_j\}\) and the dual conic
> \(\{q_i^*=q_j^*\}\) have discriminants of the same sign and \(3\times3\) determinants that vanish
> together. Ellipse goes to ellipse, hyperbola to hyperbola, line pair to line pair.

Two convenient closed forms follow, obtained by combining Lemma B with the scaling rules
(multiplying a conic by \(t\) scales \(\det A_2\) by \(t^2\) and \(\det M\) by \(t^3\)); they are
what makes the whole question a hand computation. For the unscaled \(g_i-g_j\),

\[
b^2-4ac \;=\; -\,\frac{\det(Q_i-Q_j)}{\det Q_i\,\det Q_j},
\qquad\qquad
\det M(g_i-g_j) \;=\; -\,\frac{\det M_p}{8\,\det Q_i\,\det Q_j}.
\tag{5.5}
\]

---

## 6. Lemma C — why an adjacent pair can never give a curved arc

**Proposition 5 (adjacency forces degeneracy).** Suppose \(f\) is finite and continuous on
\(T_i\cup T_j\), and \(T_i\cap T_j\) contains a segment of positive length lying on the line
\(L=\{x:\ell(x)=0\}\) for some non-zero affine \(\ell\). Then

\[
q_i - q_j \;=\; \ell\cdot m \qquad\text{for some affine } m,
\]

hence \(\det M_p=0\) (Lemma 4), hence \(\det M_d=0\) (Lemma B): the conic \(Z_{ij}\) is
**degenerate**. Moreover by (2.5) it is either two distinct real lines (the generic case) or a
double line — never a genuine hyperbola and never a genuine ellipse.

*Proof.* Continuity of \(f\) on \(T_i\cup T_j\) means \(f\) is single-valued there, so
\(q_i(x)=q_j(x)\) at every \(x\in T_i\cap T_j\), in particular on a segment
\([x_0,x_1]\subseteq L\) with \(x_0\neq x_1\). Parametrise \(L\) by
\(t\mapsto x_0+t(x_1-x_0)\). The composition \(t\mapsto (q_i-q_j)\bigl(x_0+t(x_1-x_0)\bigr)\) is a
univariate polynomial of degree at most \(2\) vanishing for every \(t\in[0,1]\), hence infinitely
often, hence identically. So \(q_i-q_j\) vanishes on **all** of \(L\), not just on the shared
segment.

Now choose affine coordinates \((y_1,y_2)\) in which \(\ell(y)=y_1\), so \(L=\{y_1=0\}\). In these
coordinates \(P:=q_i-q_j\) is again a polynomial of degree \(\le 2\):

\[
P \;=\; c_{20}y_1^2 + c_{11}y_1y_2 + c_{02}y_2^2 + c_{10}y_1 + c_{01}y_2 + c_{00}.
\]

\(P\) vanishes identically on \(\{y_1=0\}\), i.e. \(c_{02}y_2^2+c_{01}y_2+c_{00}=0\) for all
\(y_2\), so \(c_{02}=c_{01}=c_{00}=0\). Therefore

\[
P \;=\; y_1\bigl(c_{20}y_1 + c_{11}y_2 + c_{10}\bigr) \;=\; \ell\cdot m,
\qquad m \text{ affine},
\]

as claimed. The rest is Lemma 4 and Lemma B. \(\blacksquare\)

Two consequences, which between them explain the reader's experience completely.

**Corollary 5a (two pieces are never enough).** Let \(f\) have two pieces whose intersection
contains a segment — the normal situation for a two-piece subdivision of a connected planar domain,
since the two polyhedra have to meet along their common boundary. Then the single pair available is
adjacent, so by Proposition 5 the only conic in sight is a pair of lines. **No such two-piece
example can exhibit a curved arc, whatever the two Hessians are.** This is exactly the family the
reader built, which is why the outcome was always the same.

**Corollary 5b (what is needed instead).** A non-degenerate arc requires a pair \((i,j)\) with
\(T_i\cap T_j\) containing no segment. Note that sharing a single *point* is harmless: it imposes
the one scalar equation \(q_i=q_j\) at that point, not the two-parameter family of equations that
forces a factorisation. The first configuration that admits such a pair is **three pieces with the
middle one separating the other two** — which is precisely `doc/QuaConExample.md` §2. (In §4.1's
five-piece fan, pieces \(3\) and \(5\) share only the apex \((0,0)\), where both quadratics
vanish.)

---

## 7. The theorem, and the whole example in one table

Collecting Lemmas A, B and Corollary 3:

**Theorem 6.** Let \(T_i,T_j\) be two pieces with positive definite Hessians \(Q_i,Q_j\), and let
\(M_p\) be the primal \(3\times3\) matrix (3.2). Then the conic \(Z_{ij}=\{g_i=g_j\}\) separating
the two interior faces is

* of **hyperbolic** type iff \(\det(Q_i-Q_j)<0\), of elliptic type iff \(>0\), parabolic iff
  \(=0\);
* **non-degenerate** iff \(\det M_p\neq 0\);

and therefore a **genuine hyperbola** iff \(\det(Q_i-Q_j)<0\) **and** \(\det M_p\neq 0\). Both
tests are integer determinants of the input data; no conjugate need be computed.

Applied to `doc/QuaConExample.md` §4.1, whose data are

\[
\begin{array}{llll}
Q_1=\begin{pmatrix}3&-1\\-1&5\end{pmatrix} & \det Q_1 = 14
  & \beta_1=(0,1) & \gamma_1=0\\
Q_2=\begin{pmatrix}3&-5\\-5&17\end{pmatrix} & \det Q_2 = 26
  & \beta_2=(6,-8) & \gamma_2=0\\
Q_3=\begin{pmatrix}15&-2\\-2&11\end{pmatrix} & \det Q_3 = 161
  & \beta_3=(2,-6) & \gamma_3=0\\
Q_4=\begin{pmatrix}7&2\\2&17\end{pmatrix} & \det Q_4 = 115
  & \beta_4=(8,-3) & \gamma_4=0\\
Q_5=\begin{pmatrix}11&11\\11&35\end{pmatrix} & \det Q_5 = 264
  & \beta_5=(10,0) & \gamma_5=0
\end{array}
\]

the theorem gives the following for all **ten** pairs of the five pieces. The last two columns are
the closed forms (5.5) and may be compared line by line with the tables in
`doc/QuaConExample.md` §3.3 and §4.2.

| pair | adjacent? | \(\det(Q_i-Q_j)\) | \(\det M_p\) | verdict | \(b^2-4ac\) | \(\det M\) of \(g_i-g_j\) |
|---|---|---|---|---|---|---|
| \(1\,\vert\,2\) | yes | \(-16\) | \(\mathbf{0}\) | two lines | \(4/91\) | \(0\) |
| \(2\,\vert\,3\) | yes | \(-81\) | \(\mathbf{0}\) | two lines | \(81/4186\) | \(0\) |
| \(3\,\vert\,4\) | yes | \(-64\) | \(\mathbf{0}\) | two lines | \(64/18515\) | \(0\) |
| \(4\,\vert\,5\) | yes | \(-9\) | \(\mathbf{0}\) | two lines | \(3/10120\) | \(0\) |
| \(1\,\vert\,3\) | no | \(+71\) | \(584\) | **ellipse** | \(-71/2254\) | \(-73/2254\) |
| \(1\,\vert\,4\) | no | \(+39\) | \(1024\) | **ellipse** | \(-39/1610\) | \(-64/805\) |
| \(1\,\vert\,5\) | no | \(+96\) | \(3248\) | **ellipse** | \(-2/77\) | \(-29/264\) |
| \(2\,\vert\,4\) | no | \(-49\) | \(-40\) | **hyperbola** | \(49/2990\) | \(1/598\) |
| \(2\,\vert\,5\) | no | \(-112\) | \(-224\) | **hyperbola** | \(7/429\) | \(7/1716\) |
| \(3\,\vert\,5\) | no | \(-265\) | \(144\) | **HYPERBOLA** | \(265/42504\) | \(-3/7084\) |

Read the table twice.

1. **Every one of the four adjacent pairs has \(\det M_p=0\) exactly, and hyperbolic type.** That
   is Proposition 5 in action: adjacent comparisons give *degenerate hyperbolas*, i.e. pairs of
   straight lines, and they can give nothing else. This is what the reader has been observing, and
   it is a theorem about her construction, not evidence about ours.
2. **Every one of the six non-adjacent pairs has \(\det M_p\neq 0\)** — non-degeneracy is the
   generic situation the moment adjacency stops forcing the issue — and their types are then
   decided by \(\det(Q_i-Q_j)\): three ellipses and three hyperbolas. For the three pairs that meet
   at the triple point \(p\), the rows reproduce exactly the six numbers printed in
   `doc/QuaConExample.md` (\(-71/2254\), \(-2/77\), \(265/42504\) for \(b^2-4ac\), and
   \(-73/2254\), \(-29/264\), \(-3/7084\) for the \(3\times3\) determinants), recomputed here from
   the primal data through (4.1) and (5.5) without touching a conjugate. That agreement is a strong
   consistency check on the document in both directions.
3. **A caveat on how much the table claims.** It classifies the conic \(Z_{ij}\) for every pair. It
   does *not* claim that every pair actually contributes an edge to \(f^*\): whether the two
   interior cells \(C_i,C_j\) meet, and whether some third piece takes the maximum first, is
   separate information, and `doc/QuaConExample.md` §3.1 records that six of the 23 formulas never
   win anywhere. Section 9 below is what establishes that the \(3\,|\,5\) conic is genuinely
   realised as an edge of \(f^*\). The table's role is the converse direction, which is the one
   under dispute: for an adjacent pair no amount of geometry could produce a curved edge, because
   there is no curve to produce.

---

## 8. The witness, spelled out

Pieces \(3\) and \(5\). They are separated by pieces \(4\) and \(2\) and meet only at the apex, so
Proposition 5 does not apply and nothing constrains \(q_3-q_5\).

**Type.** \(Q_3-Q_5=\begin{pmatrix}4&-13\\-13&-24\end{pmatrix}\). Its indefiniteness needs no
eigenvalues, just two vectors:

\[
\bigl\langle e_1,(Q_3-Q_5)e_1\bigr\rangle = 4 > 0,
\qquad
\bigl\langle e_2,(Q_3-Q_5)e_2\bigr\rangle = -24 < 0,
\]

so \(Q_3-Q_5\) is indefinite; concretely
\(\det(Q_3-Q_5)=4\cdot(-24)-(-13)^2=-96-169=-265<0\). Piece 3 is the stiffer one along \(x_1\),
piece 5 along \(x_2\) — neither dominates the other. Hyperbolic type. By (4.1) the arc's
discriminant is

\[
b^2-4ac \;=\; \frac{265}{161\cdot 264} \;=\; \frac{265}{42504} \;>\; 0 .
\]

**Non-degeneracy.** \(\beta_3-\beta_5=(-8,-6)\) and \(\gamma_3-\gamma_5=0\), so

\[
M_p \;=\; \begin{pmatrix} 4 & -13 & -8 \\ -13 & -24 & -6 \\ -8 & -6 & 0 \end{pmatrix},
\]
\[
\det M_p \;=\; 4(0-36) - (-13)(0-48) + (-8)(78-192) \;=\; -144 - 624 + 912 \;=\; 144 \;\neq\; 0 .
\]

By Lemma B, \(\det M_d=-144/42504\neq 0\). **By Theorem 6 the arc \(\{g_3=g_5\}\) is a
non-degenerate hyperbola.** \(\blacksquare\)

**The curve itself.** Written out from Lemma 1 and cleared of denominators by
\(2\cdot 42504 = 2\det Q_3\det Q_5\), \(Z_{35}\) is

\[
-2731\,s_1^2 + 4598\,s_1s_2 + 2189\,s_2^2 + 107420\,s_1 + 9988\,s_2 - 421996 \;=\; 0 ,
\]

which is the equation in `doc/QuaConExample.md` §3.3. Its invariants are

\[
\begin{aligned}
A_2 &= \begin{pmatrix} -2731 & 2299 \\ 2299 & 2189 \end{pmatrix},
 & \det A_2 &= -11\,263\,560 = -265\cdot 42504, \\
\det M &= -260\,148\,962\,304 = -144\cdot 42504^2,
 & c &= \Bigl(\tfrac{2496}{265},\, -\tfrac{3226}{265}\Bigr)
      = (9.41886\ldots,\, -12.17358\ldots), \\
k &= \frac{\det M}{\det A_2} = \frac{6\,120\,576}{265} = 23\,096.5132\ldots
\end{aligned}
\]

— the two Lemma A / Lemma B predictions confirmed on the nose (the powers of \(42504\) are the
scaling rules for \(\det A_2\) and \(\det M\)). Since \(k\neq 0\), Corollary 3 applies: the
eigenvalues of \(A_2\) are \(-271-\sqrt{11\,337\,001}\) and \(-271+\sqrt{11\,337\,001}\), i.e.
\(-3638.0463\ldots\) and \(+3096.0463\ldots\), and the canonical form is

\[
\frac{u^2}{a^2} - \frac{v^2}{b^2} \;=\; 1,
\qquad a = 2.51964\ldots,\quad b = 2.73130\ldots,\quad e = 1.47481\ldots
\]

with asymptote slopes the roots of \(2189\,m^2+4598\,m-2731=0\), namely \(0.4829249\ldots\) and
\(-2.5834274\ldots\) (the discriminant \(45\,054\,240\) is not a perfect square, so the slopes are
irrational conjugates). Had \(\det M\) been \(0\), \(k\) would have been \(0\) and those two slopes
would have *been* the curve.

**Why it looks straight, and why no plot can settle this.** The realised edge of \(f^*\) — verified
in §9 — runs from the triple point outwards, and every point of it sits far out on one branch: the
sampled endpoints are at distances \(15.287\) and \(55.970\) from the centre, i.e. \(6.1\) and
\(22.2\) semi-axes \(a\). A hyperbola is asymptotically straight, so over that stretch the arc is
nearly a line segment: with the two endpoints joined by a chord of length \(40.68\), the maximum
perpendicular deviation of the arc from that chord is \(0.0479\), about **1 part in 850**. Its
curvature at the midpoint is \(1.18\times 10^{-3}\). An eyeball check, or a plot at any reasonable
scale, will call this a straight line. **The determinant will not.**

---

## 9. What was recomputed, and with what

Independent recomputation on 2026-09-10 (MATLAB R2024b on the UBC VPN; Python 3 / SymPy 1.14;
SNAT at `../SNAT`; Lean 4.33 + Mathlib).

* **CCA2, its own code path.** `QuaPol` accepts §4.1's five rows and reports
  `isEdgeContinuous = 1` on all four cevians and `isConvex = 0` — the input is a continuous
  nonconvex PLQ, as claimed. `conjConvexOverPiece` (the KKT active-set construction, not the
  `conj('cplq')` pipeline, which still dies at its Step-3 gate — a separate finding) rebuilds all
  five interior faces, and they agree character for character with §3.2's \(I_1,\dots,I_5\). From
  the two relevant formulas: the conic above, with \(b^2-4ac = 45\,054\,240\) and
  \(\det M = -260\,148\,962\,304\) for the cleared form; and both Lemma A and Lemma B verified
  symbolically (\(-265\) versus \(-265\), \(-144\) versus \(-144\)). At ten points of the conic,
  \(f^*=\max_k g_k\) with \(g_3-g_5=0\) to machine zero and a margin of \(0.18\) to \(0.97\) over
  pieces \(1,2,4\) — so the arc really is an edge of \(f^*\) and both faces attain the max there.
* **An oracle sharing no code with any of the above** (40-digit `mpmath`; per piece, the sup over
  the closed triangle by enumeration — interior stationary point when feasible, each closed edge as
  a clamped 1-D concave quadratic, plus the vertices). The \(I_3|I_5\) boundary of \(f^*\) was
  located by **bisecting on the argmax**, never using the conic; the located points satisfy the
  conic equation to \(10^{-35}\). Three of them span a triangle of area \(0.0085\)
  (\(2\times\text{area}=0.017\neq 0\)): the edge contains **three non-collinear points**, which a
  pair of lines could not do. Boundary points were found for
  \(s_2 = 2,\,2.5,\,\dots,\,6,\,10,\,20,\,40\) and none at \(s_2=80\).
* **SNAT** (its own convexity-agnostic 1-D conjugate engine, no CCA2 code): slicing \(f\) at fixed
  height \(y\) and conjugating the 1-D *nonconvex* PLQ slice, then maximising over \(y\). At
  \(s=(17/5,\,5/2)\) the outer maximisers sit at \(y=1303/1610\) (piece 3) and \(y=91/240\)
  (piece 5) and SNAT returns \(115291/32200=g_3(s)\) and \(3557/960=g_5(s)\) **exactly, as
  rationals**; at \(s=(18/5,\,5/2)\) it returns \(116631/32200=g_3(s)\) and
  \(37087/10560=g_5(s)\). The winner swaps between the two points, which is the arc being crossed.
  Both face formulas therefore hold in a second, independent engine.
* **Lean 4 + Mathlib**, compiled clean with no `sorry` (axioms only `propext`, `Classical.choice`,
  `Quot.sound`): Lemma 4 (`det3_of_line_pair`, and `det2_of_line_pair` for the identity
  \(\det A_2=-((a_1b_2-a_2b_1)/2)^2\)), Lemmas A and B in division-free polynomial form
  (`transfer_det2`, `transfer_det3`), the derivation of the printed integer conic from the primal
  entries (`dual_coeffs_are_the_printed_conic`), and the two non-factorisation statements
  (`primal_not_a_line_pair`, `dual_not_a_line_pair`).

---

## 10. If you want to build your own

Theorem 6 is a recipe. Pick two positive definite \(2\times2\) integer matrices with

\[
\det(Q_i-Q_j) < 0
\qquad\text{(equivalently: neither } Q_i\succeq Q_j \text{ nor } Q_j\succeq Q_i\text{)}
\]

— easiest by making one stiffer in \(x_1\) and the other stiffer in \(x_2\), as \(Q_3,Q_5\) are —
then pick \(\beta_i,\beta_j,\gamma_i,\gamma_j\) and check \(\det M_p\neq 0\), which is generic.
Finally arrange the geometry so the two pieces **share no segment**: a fan of three or more
triangles from a common apex, with the pair you care about separated by at least one other piece,
is the smallest such arrangement. Then verify, as §9 does, that the two interior cells really do
meet — non-degeneracy of the conic is a statement about the curve, and one still has to check that
\(f^*\) actually switches between those two faces somewhere rather than handing the maximum to a
third piece first.

A warning worth repeating: if you compare the two Hessians and get \(\det(Q_i-Q_j)>0\) you will
find a perfectly genuine **ellipse** arc, as pairs \(1|3\), \(1|4\) and \(1|5\) do. That is also a
refutation of a parabolic subdivision, and it is in fact what `doc/QuaConExample.md` §2's
three-piece minimal example produces. The hyperbola needs the extra condition that the two Hessians
be Loewner-incomparable.

---

_Companion documents: `doc/QuaConExample.md` (the computed conjugate, the 17 faces and the 41
edges), `CONJ_FIELD_PROOF.md` §4 (the input and the \(S_4\) vertex) and its Theorem 3 (the adjacent
case, of which Proposition 5 here is a self-contained proof)._
