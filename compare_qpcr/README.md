README
================

Data
----

qPCR dataset used in papers such as [Ramsuram et al, HMG, 2015](https://doi.org/10.1093/hmg/ddv158), [Ramsuram et al, JI, 2017](http://www.jimmunol.org/content/198/6/2320), and [Kulkarni et al, PNAS, 2013](https://doi.org/10.1073/pnas.1312237110), were obtained from Mary Carrington by personal communication.

Comparison qPCR vs RNAseq (HLApers)
-----------------------------------

![](./FigSA.png)

Statistical analysis
--------------------

We selected lineages with ≥10 individuals in both datasets. Then we performed a pairwise Mann-Whitney U test in qPCR and RNAseq separately to identify pairs of alleles with significantly different distribuitions of expression values (FDR ≤ 5%).

We then selected all significant pairs in the qPCR dataset, we checked their status in the RNAseq dataset.

In the table below,

`a1` and `a2`: a pair of lineages;

`p_qpcr` and `p_hlapers`: p-value\* from a pairwise Mann-Whitney U test for each lineage pair in qPCR and RNA-seq data, respectively;

`direction_qpcr` and `direction_hlapers`: (1) if expression of `a1` &gt; expression of `a2`, or (-1) if expression of `a1` &lt; expression of `a2`.

| a1    | a2    |  p\_qpcr|  p\_hlapers|  direction\_qpcr|  direction\_hlapers|
|:------|:------|--------:|-----------:|----------------:|-------------------:|
| A\*01 | A\*02 |  0.00000|     0.00000|                1|                  -1|
| A\*01 | A\*03 |  0.00000|     0.00037|                1|                  -1|
| A\*01 | A\*25 |  0.00095|     0.00088|                1|                  -1|
| A\*01 | A\*31 |  0.00009|     0.04162|                1|                  -1|
| A\*02 | A\*03 |  0.00000|     0.00000|                1|                   1|
| A\*02 | A\*24 |  0.00000|     0.00000|               -1|                   1|
| A\*02 | A\*32 |  0.04117|     0.00023|                1|                   1|
| A\*03 | A\*11 |  0.00000|     0.00476|               -1|                   1|
| A\*03 | A\*30 |  0.00406|     0.01958|               -1|                   1|
| A\*11 | A\*31 |  0.04939|     0.03705|                1|                  -1|
| A\*24 | A\*25 |  0.00095|     0.01958|                1|                  -1|
| A\*24 | A\*68 |  0.02467|     0.02241|                1|                  -1|
| A\*30 | A\*68 |  0.04117|     0.00124|               -1|                  -1|
| A\*32 | A\*68 |  0.00256|     0.02241|               -1|                  -1|
| B\*07 | B\*35 |  0.02257|     0.03713|               -1|                  -1|
| B\*08 | B\*35 |  0.02257|     0.02422|               -1|                   1|
| B\*15 | B\*35 |  0.04862|     0.02441|               -1|                  -1|
| B\*35 | B\*40 |  0.04605|     0.01228|                1|                   1|
| C\*01 | C\*03 |  0.00083|     0.00122|                1|                   1|
| C\*02 | C\*03 |  0.00296|     0.00000|                1|                   1|
| C\*03 | C\*04 |  0.00000|     0.00000|               -1|                  -1|
| C\*03 | C\*05 |  0.00003|     0.00002|               -1|                  -1|
| C\*03 | C\*06 |  0.00000|     0.00000|               -1|                  -1|
| C\*03 | C\*07 |  0.01673|     0.00000|               -1|                  -1|
| C\*03 | C\*08 |  0.00002|     0.00849|               -1|                  -1|
| C\*03 | C\*12 |  0.00229|     0.01160|               -1|                  -1|
| C\*03 | C\*15 |  0.00271|     0.00034|               -1|                  -1|
| C\*04 | C\*07 |  0.00000|     0.00000|                1|                   1|
| C\*04 | C\*12 |  0.01816|     0.00000|                1|                   1|
| C\*06 | C\*07 |  0.00000|     0.00000|                1|                   1|
| C\*06 | C\*12 |  0.01430|     0.00006|                1|                   1|

Finally, we calculated the number of pairs with the same or different ordering in qPCR and RNAseq.

| locus |  same\_dir|  n\_significant|
|:------|----------:|---------------:|
| A     |          4|              14|
| B     |          3|               4|
| C     |         13|              13|
