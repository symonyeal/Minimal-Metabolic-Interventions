# Where MMCS spends its search: a census on metabolic cut-set instances

Version 0.1.0, 2026-10-27. Frozen measurement; the numbers below are reproduced by the commands in
[EVIDENCE.md](EVIDENCE.md) and stored in [`results/`](results/).

## 1. Question

MMCS (Murakami and Uno 2014) enumerates the inclusion-minimal hitting sets `Tr(H)` of a hypergraph
`H` by depth-first search. This study asked whether search states that are equivalent up to vertex
relabelling recur often enough, after the established reductions, to justify reusing their
completions. The predeclared decision rule: a reuse method is worth building only if the subtrees
it could skip hold at least 1/6 of the search nodes, the smallest share compatible with a 1.2x
speedup before any detection cost (`1/(1 − f) ≥ 1.2`). The census measures that share without
pruning anything.

## 2. Formulation

Notation: `ev[e]` is the vertex set of edge `e`; a search node carries the selected set `S`, the
candidate set `A`, the uncovered edges `U` and, for each `s ∈ S`, its critical (private) edges
`P_s`. MMCS branches on an uncovered edge minimizing `|ev[e] ∩ A|`, adds one candidate per child and
discards a child in which some `P_s` becomes empty (violator pruning).

**Completion family.** The subtree of node `(S, A)` emits exactly the minimal hitting sets `M` with
`S ⊆ M ⊆ S ∪ A`. Writing `M = S ∪ T`, `T` is a valid completion if and only if

1. `T ∈ Tr(U')`, where `U' = min {ev[e] ∩ A : e ∈ U}`, and
2. for every `s ∈ S`, `T` misses some member of `O_s = min {ev[e] ∩ A : e ∈ P_s}`.

Condition 1 combines coverage of `H` with a private edge for each `t ∈ T`: an edge private to `t`
must be uncovered by `S`, and only its part inside `A` can meet `T`. Condition 2 keeps each selected
vertex irredundant. Taking inclusion-minimal members (`min`) changes neither condition. An owner
`s` with an empty member of `O_s` can never lose its witness and is dropped.

The completion family is therefore determined by the **restricted state**
`(b, U', {O_s})`: the remaining cardinality budget `b`, the family `U'` and the multiset of owner
families. Owner identities, `|S|` and vertices outside `A` affect only the labelled output.

**Census quantities.** For each node the census computes a 128-bit hash of the restricted state
(exact key) and a colour-refinement invariant of the coloured incidence structure of `U'` and the
owner families (WL key). A node is a repeat root when its key equals that of an earlier node whose
subtree has finished; nested repeats are counted once. Equal WL keys are necessary for
isomorphism, so the WL share bounds isomorphic reuse from above. A node is dead when its subtree
emits nothing.

| Quantity | Definition |
|---|---|
| Dead | Share of nodes whose subtree emits no solution |
| Exact | Share of nodes inside subtrees rooted at an exact repeat |
| WL bound | Share of nodes inside subtrees rooted at a WL repeat (upper bound for isomorphic repeats) |
| Gap | WL bound − exact: the most isomorphism can add over exact memoization |

## 3. Data

| Lineage | Source | Instances |
|---|---|---|
| A | METATOOL elementary modes of an E. coli central-metabolism network on five substrates (Haus, Klamt and Stephen 2008), as published in `VeraLiconaResearchGroup/MHSGenerationAlgorithms` @ `f043cc3` | Growth instances: modes containing biomass reaction `mue`, with `mue` removed; all other reactions eligible. 30, 103, 274, 879 and 2,891 edges. |
| B | BiGG `e_coli_core` (Orth, Fleming and Palsson 2010); 100,274 elementary modes under the model's bound signs, computed here by `efm` | Ethanol instance: the 48,770 modes using `EX_etoh_e`; exchange, demand, sink, biomass and ATP maintenance reactions and reactions whose bounds exclude zero are not eligible (69 eligible). |

Lineage B is an independent reconstruction of the same organism. No tractable mode family of a
different organism was found: the smallest non-E. coli BiGG model (`iAB_RBC_283`) has a flux-cone
kernel of dimension 130, against 24 for `e_coli_core`. No claim across organisms is made.

## 4. Validation

| Check | Result |
|---|---|
| Brute force, 406 hypergraphs, all algorithm settings, census on, `k = 0..3` | Exact families |
| Python Berge oracle (Hypergraph-Dualization), small example, Fano plane, A acetate-growth | Exact families (5, 7, 617) |
| SHD (Uno), all five A growth instances and four unfiltered A instances | Exact families (617 to 73,760) |
| SHD, B ethanol | Exact family (34,389) |
| LP oracle on the model, independent of the modes, B growth and ethanol, cuts of size ≤ 2 | Exact families (128, 67) |
| `efm` mode count against the published count for `e_coli_core` | 100,274 = 100,274 |
| Every B mode: nullity 1 and sign-feasible null vector (500 also in exact arithmetic) | 100,274 pass |

The census does not alter the search: the families with the census on equal those with it off.

## 5. Results

Comparator configuration: static twins merged, candidates in increasing unhit degree
(Hörmann and Schirneck 2026). All runs complete.

| Instance | Nodes | Dead | Exact | WL bound | Gap |
|---|---:|---:|---:|---:|---:|
| A acetate-growth | 77 | 26.0% | 1.3% | 5.2% | 3.9% |
| A succinate-growth | 476 | 60.3% | 15.8% | 22.3% | 6.5% |
| A glycerol-growth | 1,162 | 62.8% | 18.2% | 25.3% | 7.1% |
| A glucose-growth | 3,901 | 76.6% | 21.0% | 28.0% | 6.9% |
| A combined-growth | 38,391 | 79.3% | 15.7% | 19.8% | 4.1% |
| B ethanol | 41,948 | 77.8% | 35.9% | 42.7% | 6.8% |

Without twins the gap is 13.7% to 27.0% (all 24 configurations in
[`results/census-lineage-a.md`](results/census-lineage-a.md) and
[`results/census-lineage-b.md`](results/census-lineage-b.md)). Static twins remove most of the
symmetry that an isomorphism test would otherwise find.

Findings:

1. After twins and unhit ordering, isomorphic repeats add at most 7.1 percentage points over exact
   repeats on every instance, below the 1/6 threshold. Isomorphic continuation reuse was stopped
   at this gate.
2. Dead subtrees hold 26% to 79% of the nodes; on all instances but acetate-growth, more than 60%.
3. Exact repeats hold 1.3% to 35.9% of the nodes. On B ethanol, 13,171 of 41,948 nodes (31%) lie in
   exact repeats of dead states, the share an exact cache of infeasible states could skip at zero
   cost.
4. Computing restricted keys at every node costs about 100 times a plain search node in this
   implementation, so the census shares are upper bounds on savings, not measured speedups.

## 6. Decision record

| Decision | Reason |
|---|---|
| Stop isomorphic continuation reuse | Gap at most 7.1 points after twins; the rule required 1/6 |
| Publish the tool and the census | The census is a measurement independent of whether a method succeeds |
| Exact dead-state caching not pursued here | Exact memoization of subproblems has prior art (operation caches in Toda 2013; component caching in model counters); a claim needs a focused comparison and cheap incremental keys |

## 7. Limitations

- Two lineages, one organism. Lineage A substrates share one network and are not independent samples.
- The B growth census (196,481 nodes) did not finish within 240 s and is not reported.
- WL equality is an upper bound for isomorphism; exact isomorphic shares were not computed.
- Times are from one Windows 11 laptop and include process start; no benchmark-machine runs.
- Census shares count nodes, not time; per-node cost varies with the number of uncovered edges.

## References

- Haus, Klamt and Stephen (2008). Computing knock-out strategies in metabolic networks. *J. Comput. Biol.* 15(3). arXiv:0801.0082.
- Murakami and Uno (2014). Efficient algorithms for dualizing large-scale hypergraphs. *Discrete Appl. Math.* 170. arXiv:1102.3813.
- Hörmann and Schirneck (2026). Minimal-to-maximal conversion search is not output-polynomial. arXiv:2608.02159.
- Gainer-Dewar and Vera-Licona (2017). The minimal hitting set generation problem: algorithms and computation. *SIAM J. Discrete Math.* 31(1). arXiv:1601.02939.
- Toda (2013). Hypergraph transversal computation with binary decision diagrams. SEA 2013, LNCS 7933.
- Orth, Fleming and Palsson (2010). Reconstruction and use of microbial metabolic networks: the core *Escherichia coli* metabolic model as an educational guide. *EcoSal Plus*.
- Clement, Baalhuis, Teusink, Bruggeman, Planqué and de Groot (2020). Elementary conversion modes enable exploration of all metabolic capabilities of any organism. bioRxiv 10.1101/2020.06.06.137554.
