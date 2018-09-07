This repository is a counterpart to our paper, ARXIV LINK, where we give additional details of the experiment described in Section 6. This includes:
1. A description of the random CP-net generator used, as well as the code for this generator
2. An explanation of the different leaf prioritisations tested and which prioritisations were used in the presented results in ARXIV LINK.
3. The code for our different dominance testing functions.
4. The raw results of our experiments (using all possible prioritisation options).
5. Plots of these results, showing average perfromance of the different dominance testing functions. These are similar to those given in the paper, only now presenting the results for all prioritisation options.

# Experiment
As detailed in our paper, ARXIV LINK, we conducted an experiment to compare the performance of several dominance testing functions for CP-Nets.
Each dominance testing function answers dominance queries by building the associated search tree.
To improve efficiency, a pruning method is used to prune certain branches of this tree as it is constructed.
In order to fully specify how a dominance testing function works, one must also specify the method of leaf prioritisation used (when building up the tree).
The pruning method and prioritisation technique combinations we tested are listed in the Dominance Tessting Functions section below.
The performance of a dominance testing function on a given query is measured both by outcomes traversed and time elapsed.
For full details of this terminology, we refer you to section 6 of the paper.

Given *n* (number of variables in the CP-net), 100 CP-nets with *n* variables were randomly generated. Then, for each CP-net, 10 dominance queries were randomly generated. Each dominance testing function was then applied to this set of 1000 dominance queries and performance of the functions was recorded in each case.

This was done twice. First, for the case of binary CP-nets, for *n*=3-10.
Second, for the case of multivalued CP-nets, for *n*=3-8. By multivalued CP-nets, we mean that the variables could have domain size between 2-5.

In this repository, we provide the CP-net generator code we used. We also provide the code for each of the dominance testing functions.
These are all provided as R scripts and allow you to repeat these experiments with different values of *n*.
You may also vary the values of the maximum possible size of variable domain and the maximum number of parents a variable can have (we left this as *(n-1)* in our experiments, which does not impose any condition upon the generated CP-nets).

# CP-Net Generator
The CP-net generator we used for these experiments takes some inspiration from the CP-net generator described by Allen et al. (2016). In particular, the usage of the dagcode representation of DAGs, and the test for conditional preference table (CPT) degeneracy. However, our generator does not guarantee a specific distribution over the generated CP-nets. This is because the probability vectors required by Allen et al. (2016), to make their generation uniform, required more precision than was possible for our computational resources. Our generator also differs from theirs in that it allows variables to have different domain sizes (whereas Allen et al. (2016) requires all variables to have the same domain size).

Our CP-net generator function, `rand.cpn`, is given in the R script `CPNGenerator.R`, along with all of its dependent functions. Note that these functions require the R libraries `primes` and `Rmpfr`. In this script, we also give an illustrative example of what the output CP-nets look like as R objects. We also show how to generate an associated outcome. These are the required inputs for the dominance testing functions.

The CP-net generator function works as follows:
1. Generate a valid dagcode.
2. Convert this dagcode to the adjacency matrix, *A*, of the associated DAG.
3. `for` each variable: Generate a random CPT. Let *CPT* denote this list of *n* CPTs.
4. `if` any of the CPTs in *CPT* are degenerate `then` go back to step 3
5. Return the CP-net as the pair (*A*, *CPT*).

For the theory of dagcodes and obtaining DAGs from dagcodes, we refer you to Allen et al. (2016). 

A CPT, say the CPT of variable *X* (CPT(*X*)), is degenerate (invalid) if it implies that one of the parents of *X* is not valid. A parent, *Y*, of *X* is valid if you can change the user's preference over *X* by changing the value taken by *Y only*. Thus, a CPT of *X* is valid (non-degenerate) if, for every parent *Y*, there are two rows where the parental assignment differs only on the value taken by *Y* and the preference over *X* is different in these two rows.

The function `rand.cpn(n,p,d)` outputs a random CP-net with *n* variables, where each variable has no more than *p* parents, and each variable has a domain size in 2-*d*.

# Dominance Testing Functions
## Pruning Methods
We tested three different pruning methods:
* Rank Pruning - ARXIV LINK
* Penalty Pruning - (Li et al., 2011)
* Suffix Fixing - (Boutilier et al., 2004)

We tested each of these methods used individually, all possible pairwise combinations, and all three methods used together.
This gives us 7 pruning schema options.
Again, we refer you to section 6 of ARXIV LINK for a detailed description of how these three pruning methods work, as well as how they can be combined.

## Leaf Prioritisation methods
There have been several methods of leaf prioritisation suggested in the existing literature.
However, there has been no experimental analysis of how effective the different methods are.
We decided, for the purposes of keeping our dominance testing functions efficient, to only implement leaf prioritisation methods that do not require additional calculations.
Thus, we consider the following prioritisation heuristics:
* Minimal Depth - Simply selects a leaf at minimal depth in the current search tree.
* Penalty Prioritisation (Li et al., 2011) - Selects a leaf which has minimal *f* (evaluation function for penalty pruning) value.
* Rank Prioritisation - Selects a leaf with maximal rank, *r*, value.
* Rank + Difference Prioritisation - Selects a leaf, *o*, with maximal *r(o) + L_D(o, o_1)* value (for the query '*o_1>o_2*?').

The latter two are our suggested heuristics, based upon our rank values introduced in ARXIV LINK.
Search directions with higher rank (or rank + diff.) values are more likely to either be successful in reaching o_1 or to terminate quickly.
Thus, it makes sense to prioritise these directions.

## Combinations
As we mentioned above, only prioritisation methods that do not require additional calculations were tested.
Below we show, for each pruning schema, which leaf prioritisation methods were tested.
Note that, as minimal depth is not an intelligent heuristic, it was tested only for the pruning method that did not permit the others (suffix fixing alone).

(R- rank pruning, P- penalty pruning, S- suffix fixing)
* *Pruning Method - Tested Prioritisation Hueristics*
* R - Rank Prioritisation, Rank + Diff. Prioritisation
* P - Penalty Prioritisation
* S - Min. Depth
* R & P - Rank Prioritisation, Rank + Diff. Prioritisation, Penalty Prioritisation
* R & S - Rank Prioritisation, Rank + Diff. Prioritisation
* P & S - Penalty Prioritisation
* R & P & S - Rank Prioritisation, Rank + Diff. Prioritisation, Penalty Prioritisation

This makes the total number of tested functions 13.

In ARXIV LINK, for each pruning option, only one function's performance was presented. For those pruning methods with multiple prioritisation options, it is the function that uses rank pruning that is presented in the paper. The reasoning for this is given below in the Results Plots section.

## R Functions
The R functions to answer dominance queries using these different methods are given in the R Script `DominanceTestingFunctions.R`.
Suppose we wish to answer the dominance query 'o1>o2 ?' (o1, o2 outcomes) for the CP-net N.
Let `A` be the adjacency matrix of the structure of N and let `CPT` be a list of the conditional preference tables (CPTs) of N.
To answer this dominance query using the different pruning and prioritisation combinations above, use the following R functions, found in `DominanceTestingFunctions.R`:

* **Rank Pruning, Rank Prioritisation**\
`DQ.Rank(o1,o2,A,CPT,priority="rank",suffix=FALSE,dig)`
* **Rank Pruning, Rank + Difference Prioritisation**\
`DQ.Rank(o1,o2,A,CPT,priority="rank.diff",suffix=FALSE,dig)`
* **Penalty Pruning, Penalty Prioritisation**\
`DQ.Penalty(o1,o2,A,CPT,suffix=FALSE,dig)`
* **Suffix Fixing, Minimal Depth Prioritisation**\
`DQ.SF(o1,o2,A,CPT,dig)`
* **Rank Pruning + Penalty Pruning, Rank Prioritisation**\
`DQ.PR(o1,o2,A,CPT,priority="rank",suffix=FALSE,dig)`
* **Rank Pruning + Penalty Pruning, Rank + Difference Prioritisation**\
`DQ.PR(o1,o2,A,CPT,priority="rank.diff",suffix=FALSE,dig)`
* **Rank Pruning + Penalty Pruning, Penalty Prioritisation**\
`DQ.PR(o1,o2,A,CPT,priority="penalty",suffix=FALSE,dig)`
* **Rank Pruning + Suffix Fixing, Rank Prioritisation**\
`DQ.Rank(o1,o2,A,CPT,priority="rank",suffix=TRUE,dig)`
* **Rank Pruning + Suffix Fixing, Rank + Difference Prioritisation**\
`DQ.Rank(o1,o2,A,CPT,priority="rank.diff",suffix=TRUE,dig)`
* **Penalty Pruning + Suffix Fixing, Penalty Prioritisation**\
`DQ.Penalty(o1,o2,A,CPT,suffix=TRUE,dig)`
* **Rank Pruning + Penalty Pruning + Suffix Fixing, Rank Prioritisation**\
`DQ.PR(o1,o2,A,CPT,priority="rank",suffix=TRUE,dig)`
* **Rank Pruning + Penalty Pruning + Suffix Fixing, Rank + Difference Prioritisation**\
`DQ.PR(o1,o2,A,CPT,priority="rank.diff",suffix=TRUE,dig)`
* **Rank Pruning + Penalty Pruning + Suffix Fixing, Penalty Prioritisation**\
`DQ.PR(o1,o2,A,CPT,priority="penalty",suffix=TRUE,dig)`

Note that `dig` is a specified level of precision. How this is calculated for the CP-net of interest is given in `DominanceTestingFunctions.R`. Each of the above functions relies on a set of minor functions (`DP`,`ancestor`,`n.val`,`Rank`,`Penalty`,and `Pen.Rank`) that are given at the start of `DominanceTestingFunctions.R` and must be loaded before the dominance testing functions may be used. These functions also require two R libraries, `primes` and `Rmpfr`.

The above dominance testing functions output three results: the outcome of the dominance query, the number of outcomes traversed, and the time taken to answer the query.

# Results
The results of the experiments described above are given in the file `Raw Performance Results`. Recall that we performed the following experiment twice (once in the case of binary CP-nets, `d=2`, with `n=3-10` and once in the case of multivalued CP-nets, `d=5`, with `n=3-8`). First, we generated 100 CP-nets. The *kth* CP-net is given in the file `CPN.npd.k.RData` (`p=n-1`, `k` in 1-100). The format of these CP-nets is as discussed in `CPNGenerator.R`.

For each of these CP-nets, 10 dominance queries were generated. The 10 dominance queries for the `CPN.npd.k.RData` CP-net can be found in the file `DQ.npd.k.RData`. A dominance query asks if '*o1>o2*' is entailed. The 10 *o1* outcomes are given in order in the vector `O1` (which has length *10 x n*). Similarly, the *o2* outcomes are given in vector `O2`.

Each of these 1000 queries was then answered by the 13 different dominance testing functions described above. The query outcomes and performance of the function for the 10 queries in `DQ.npd.k.RData` (for the CP-net `CPN.npd.k.RData`) can be found in the following files:

* **Rank Pruning, Rank Prioritisation** `Rank.R.npd.k.RData`
* **Rank Pruning, Rank + Difference Prioritisation** `Rank.RDiff.npd.k.RData`
* **Penalty Pruning, Penalty Prioritisation** `Pen.npd.k.RData`
* **Suffix Fixing, Minimal Depth Prioritisation** `SF.npd.k.RData`
* **Rank Pruning + Penalty Pruning, Rank Prioritisation** `Pen.Rank.R.npd.k.RData`
* **Rank Pruning + Penalty Pruning, Rank + Difference Prioritisation** `Pen.Rank.RDiff.npd.k.RData`
* **Rank Pruning + Penalty Pruning, Penalty Prioritisation** `Pen.Rank.P.npd.k.RData`
* **Rank Pruning + Suffix Fixing, Rank Prioritisation** `Rank.SF.R.npd.k.RData`
* **Rank Pruning + Suffix Fixing, Rank + Difference Prioritisation** `Rank.SF.RDiff.npd.k.RData`
* **Penalty Pruning + Suffix Fixing, Penalty Prioritisation** `Pen.SF.npd.k.RData`
* **Rank Pruning + Penalty Pruning + Suffix Fixing, Rank Prioritisation** `Pen.Rank.SF.R.npd.k.RData`
* **Rank Pruning + Penalty Pruning + Suffix Fixing, Rank + Difference Prioritisation** `Pen.Rank.SF.RDiff.npd.k.RData`
* **Rank Pruning + Penalty Pruning + Suffix Fixing, Penalty Prioritisation** `Pen.Rank.SF.P.npd.k.RData`

Each of these files contains a list and two vectors, all of length 10. The list, `Result`, gives the outcomes of the 10 dominance queries, as found by this function (note: the `Result` list for all 13 functions should be identical). This will either be "False, N does not entail o1 > o2" or "True, N does entail o1 > o2" for each entry of `Result`. The two vectors should end in `Count` and `Time`, the begining is determined by the function used. For example, `Pen.SF.npd.k.RData` countains {`Result`, `Pen.F.SUFF.Count`, `Pen.F.SUFF.Time`}. This prefix is indicative of the function used, though some of the prefixes contain naming conventions from previous versions of this experiment, e.g. the F in the given example. The `Count` vector gives the number of outcomes considered before an answer cound be determined, for each of the 10 queries (when using the given dominance testing function). The `Time` vector gives the time elapsed before an answer cound be determined, for each of the 10 queries (when using the given dominance testing function).

Note that the data given is not complete. For the non-binary case (`d=5`), the experiments have not completed in the following cases. Suffix fixing for `n=7,8`. Penalty pruning for `n=8`. Penalty pruning + suffix fixing for `n=8`. Thus, these results are not included in the file. Once they have completed, the results files will be added. As you will see from the data/plots, these three functions perform significantly worse than the others. That is, they take considerably longer on average to answer dominance queries, which is why these experiments in particular have not yet completed.

# Result Plots
For each `n` and `d` combination, each of the 13 dominance testing functions are used to answer a set of 1000 dominance queries and the outcomes considered and time elapsed is recorded in each instance. We have averaged these results over the 1000 queries, giving us 4 data sets: Binary CP-nets - Outcomes Considered, Binary CP-nets - Time Elapsed, Multivalued CP-nets - Outcomes Considered, Multivalued CP-nets - Time Elapsed. Each data set has an average value for each dominance testing function, for each value of `n`. These are plotted on the graphs in `Result Plots`.

Each data set is presented over three plots. The prefix of the PDF identifies which data set it is for (`Binary-Outcomes`,`Binary-Time`, `Multivalued-Outcomes`,`Multivalued-Time`). Each data set has three plots, with the title suffixes `-Single`, `-Pairs`, and `-Triples`. The `Single` plot shows the performance of the single pruning methods (Rank, Penalty, Suffix Fixing), with all possible leaf prioritisations. This plot is given on a log scale to preserve legibility. The `Pairs` plot shows the performance of all pairwise combinations (Rank + Penalty, Penalty + Suffix Fixing, Rank + Suffix Fixing), with all possible leaf prioritisations. The `Triples` plot shows the performance of the dominance testing functions that use all three pruning measures. Some functions appear in multiple plots so that performance can be compared between different plots.

The shaded areas in the plot represent +/-standard error intervals. For a discussion of the variability represented by this interval, we refer you to ARXIV LINK. In the `Single` plot, it is the standard error interval of the 'Rank Pruning with Rank + Difference Prioritisation' function. In the `Pairs` plot, it is the standard error interval of the 'Rank Pruning + Suffix Fixing with Rank Prioritisation' function. In the `Triples` plot, it is the standard error interval of the 'Penalty Pruning + Rank Pruning +Suffix Fixing with Penalty Prioritisation' function.

For a discussion of these results, please consult ARXIV LINK. The additional results presented here (as opposed to in the paper) provide a comparison of the same pruning method when tested with different prioritisation techniques. For the functions where prioritisation was allowed to vary, rank prioritsation performs better than penalty prioritisation which in turn performs better than rank + difference prioritisation in the majority of cases. This is why, for all such functions, the rank prioritisation case is presented in ARXIV LINK as the optimal prioritisation choice. However, this relative performance of the different prioritisations is not always the case. Further, the difference made to performance by switching prioritisation method is generally not significant. Thus, from these results, we do not feel able to conclude that one prioritisation choice is better than another in general.

**ACKNOWLEDGEMENTS**\
The work of K. Laing was supported by a University of Leeds Research Scholarship. This work was undertaken on ARC3, part of the High Performance Computing facilities at the University of Leeds, UK.

**REFERENCES**\
Allen, T.E., Goldsmith, G., Justice, H.E., Mattei, N., and Raines, K. (2016). Generating CP-nets Uniformly at Random. *Proc. of 30th Annual Converence of the Association for the Advancement of Artificial Intelligence*, pages 872-878, Arizona, USA.

Boutilier, C., Brafman, R. I., Domshlak, C., Hoos, H. H., and Poole, D. (2004). CP-nets: A tool for representing and reasoning with conditional *Ceteris Paribus* preference statements. *Journal of Artifcial Intelligence Research*, 21:135-191.

Li, M., Vo, Q. B., and Kowalczyk, R. (2011). Efficient heuristic approach to dominance testing in CP-nets. In Tumer, Yolum, Sonenberg, and Stone, editors, *Proc. of 10th International Conference on Autonomous Agents and Multiagent Systems*, pages 353-360, Taipei, Taiwan.
