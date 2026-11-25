# Hardin and Barker's strip-coloring formulas: an extension of Britnell and Wildon's enumeration

The note proves the generating function for the normalized six-column arrays
of OEIS A233254, together with its fourth-order recurrence for heights at least
five. It also proves Barker's neighboring formulas in A233251, A233252 and
A233253, the column recurrences of A233256 for widths one through six, and the
first two row recurrences of that table. The extended results prove the
seventh-column recurrence of A233255, Barker's generating functions and closed
form in A233257 and A233258, the recurrence of A233259, and every empirical
row and column recurrence listed in A233256 through dimension seven.
Britnell and Wildon's 2015 marked-partition enumeration, Proposition 8.1,
already establishes the one-column case A233162. Its consequences for
Barker's 2018 formulas are credited explicitly; the note extends that
enumeration to graph-restricted partitions and the other strip formulas.

## Files

- `note.tex`: self-contained LaTeX source, including its bibliography.
- `note.pdf`: typeset note, supplied separately from the source archive.
- `verify.py`: standard-library Python checker, requiring Python 3.9 or later.
- `finite-checks.json`: exact expected finite-check output.
- `metadata.json`: title, creators, license, source locators and proof coverage.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
python3 verify.py --output regenerated-checks.json
cmp finite-checks.json regenerated-checks.json
```

Alternatively, run `pdflatex note.tex` twice. The bibliography is embedded;
BibTeX and additional bibliography files are unnecessary. The LaTeX packages
are `fontenc`, `lmodern`, `geometry`, `amsmath`, `amssymb`, `amsthm`, `booktabs`,
`array`, and `hyperref`.

## Mathematical coverage

The written proof identifies normalized six-colorings with partitions into
independent sets carrying binary choices. For any nonempty ordered graph
with N vertices it gives

`A(G) = 2^N (P_G(3) + 3 P_G(2) + 15 P_G(1)) / 48`.

For the six-column connected bipartite graph this reduces to
`a(n) = 64^n (P_G(3) + 6) / 48`. The displayed eight-state integer matrix and
its annihilated residual prove the formula for every positive height.

The finite checker reconstructs transitions by enumeration, checks the
matrix identities and coefficients through degree 30, independently counts
the literal normalized six-color arrays at heights 1 through 6, and tests
the graph normalization formula on all 75 simple graphs with 1 through 4
vertices. All arithmetic is exact. Finite comparisons supplement the written
all-height proof.

For fixed height h, the column states are all canonical words of length h
on at most three labels. The matrix counts next columns subject to the
horizontal and antidiagonal inequalities. Exact vector annihilation
identities in dimensions 2, 5, 14, 41, 122, and 365 imply the stated
recurrences at every subsequent width. In particular, the orders and first
widths for heights 5, 6, 7 are (7,12), (10,17), (28,35).
The finite checker verifies every coordinate of these identities,
the 16-state seventh-column identity, and coefficients through width 40
or height 30. The matrix identities are exact finite computations; the
transfer construction and infinite-sequence implications are written proofs.

In the [pinned Lean file](https://github.com/the-omega-institute/trureturing/blob/81ea1b50693b21a1d6370aa3e18815afe5f96af2/D5/S3/Combinatorics/TraceAvoidingArrays/HardinA233254TriangleStrip.lean),
namespace `D5.S3.Combinatorics.TraceAvoidingArrays.HardinA233254TriangleStrip`,
`result` proves the denominator-multiplied generating-function identity for
the literal filtered-array count `a`, with `a(0)=0`. Its
`finite_certificate` states the exact one-row class counts and weighted
moments. The graph normalization and neighboring-width arguments are written
proofs, separate from the scope of these declarations.
The same pinned source's
[three-column file](https://github.com/the-omega-institute/trureturing/blob/81ea1b50693b21a1d6370aa3e18815afe5f96af2/D5/S3/Combinatorics/TraceAvoidingArrays/HardinA233251TriangleStrip.lean)
defines `widthNormal` and `widthLegal`, the literal array predicates, and its
[pair-strip file](https://github.com/the-omega-institute/trureturing/blob/81ea1b50693b21a1d6370aa3e18815afe5f96af2/D5/S3/Combinatorics/TraceAvoidingArrays/HardinA233175PairStrip.lean)
defines the first-occurrence flags through `allowed`, `good`, and `finish`.

## Sources

The OEIS locators are A233254 revision 8, A233251 revision 9, A233252 revision
8, A233253 revision 7, and A233256 revision 6. Barker's formula is the signed
FORMULA line dated 11 October 2018 in A233254. Hardin supplied the array
definition on 6 December 2013. The note cites Salas and Sokal's 2001 article
for the general transfer-matrix method, not as a proof of these normalized
array formulas. No priority assertion is made.
The additional locators are A233255 revision 9 (2 September 2024), A233257
revision 8 and A233258 revision 8 (11 October 2018), and A233259 revision 6
(2 March 2023). A233256's unsigned empirical row and column statements are
credited to Hardin; Barker's signed formulas are credited separately.
The first abstract sentence gives the exact original six-column posing locator.
The prior one-column source is J. R. Britnell and M. Wildon,
arXiv:1507.04803v1, Proposition 8.1, pp. 23-24 and Section 6.3 (2015),
later published in Journal of Combinatorial Theory, Series A 148 (2017),
116-144, DOI 10.1016/j.jcta.2016.12.003. The other two shifted-sequence
identities use already recorded formulas in A026244 and A007582.

## Creators and license

1. Haobo Ma — ChronoAI Pte Ltd, Singapore.
2. Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

Copyright 2026 Haobo Ma and Wenlin Zhang. Licensed under
[Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/)
(CC BY 4.0; `cc-by-4.0`).
