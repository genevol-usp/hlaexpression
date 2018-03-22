Report
================

**All analyses were carried out using European individuals only (N = 358)**

Typing accuracies
=================

\*Concordance: the proportion of the called alleles that are concordant with the Gourraud et al (2014) typings

| locus |   star|  kallisto|
|:------|------:|---------:|
| A     |  97.88|     97.06|
| B     |  98.20|     98.04|
| C     |  97.22|     95.92|
| DQB1  |  98.37|     98.53|
| DRB1  |  99.51|     99.51|

Expression estimates
====================

<img src="./expression/plots/expression_boxplot.png" width="1600" />

HLA diversity vs reference transcriptome
----------------------------------------

<img src="./expression/plots/star_imgt_vs_pri_TPM.png" width="3000" />

STAR-Salmon vs kallisto
-----------------------

<img src="./expression/plots/star_vs_kallisto_TPM.png" width="2000" />

Distribution of TPM values
--------------------------

<img src="./expression/plots/tpm_distributions.png" width="2000" />

ASE
---

<img src="./expression/plots/ase.png" width="1600" />

### ASE distribution

<img src="./expression/plots/ase_histogram.png" width="1600" />

Correlation of expression
-------------------------

<img src="./expression/plots/correlation_decrease.png" width="2000" />

### Among the HLA genes

<img src="./expression/plots/correlations.png" width="3000" />

### CRDs

\*CRDs are defined in <http://dx.doi.org/10.1101/171694*>

<img src="./expression/plots/crd.png" width="1800" />

### Between Class II genes and CIITA

<img src="./expression/plots/trans_activ_corrs.png" width="2000" />

eQTLs
=====

PCA of genotypes
----------------

<img src="./qtls/plots/genotype_pca.png" width="3000" />

Number of eGenes according to index
-----------------------------------

<img src="./qtls/plots/n_of_egenes.png" width="1800" />

Distribution of eQTLs around the TSS
------------------------------------

<img src="./qtls/plots/qtls_landscape.png" width="3000" />

Spread of rank 0 eQTLs around the TSS
-------------------------------------

<img src="./qtls/plots/qtls_density_geneStart.png" width="1500" />

Investigating bias
------------------

<img src="./investigate_bias/plots/expression_by_refallele.png" width="2000" />

The last column is the partial correlation cor(expression ~ divergence | controlling for genotype). The idea here is that once we control for the genotype effect, the divergence of the HLA allele from the reference would not explain much of the variation in expression.

| index             | gene     |  expression~divergence|  divergence~genotype|  expression~divergence\_genotype|
|:------------------|:---------|----------------------:|--------------------:|--------------------------------:|
| HLA\_personalized | HLA-A    |                   0.42|                 0.59|                             0.08|
| Reference         | HLA-A    |                  -0.05|                -0.10|                            -0.13|
| HLA\_personalized | HLA-B    |                   0.21|                 0.37|                             0.04|
| Reference         | HLA-B    |                  -0.01|                 0.39|                            -0.20|
| HLA\_personalized | HLA-C    |                   0.05|                 0.34|                            -0.20|
| Reference         | HLA-C    |                  -0.12|                -0.09|                            -0.09|
| HLA\_personalized | HLA-DPB1 |                   0.60|                 0.94|                             0.14|
| Reference         | HLA-DPB1 |                   0.55|                 0.94|                             0.09|
| HLA\_personalized | HLA-DQA1 |                   0.66|                 0.70|                             0.45|
| Reference         | HLA-DQA1 |                  -0.64|                 0.11|                            -0.68|
| HLA\_personalized | HLA-DQB1 |                  -0.35|                 0.50|                            -0.24|
| Reference         | HLA-DQB1 |                  -0.52|                 0.65|                            -0.19|
| HLA\_personalized | HLA-DRB1 |                   0.60|                 0.66|                             0.31|
| Reference         | HLA-DRB1 |                   0.03|                 0.22|                            -0.15|

Intersect with Encode Elements
------------------------------

### HLA-personalized index

| locus    |  rank| rsid        |     dist| regulome\_score | tfbs                        | dhs | histone\_marks                                            |
|:---------|-----:|:------------|--------:|:----------------|:----------------------------|:----|:----------------------------------------------------------|
| HLA-A    |     0| rs28780070  |     8.68| 5               | NA                          | NA  | NA                                                        |
| HLA-A    |     1| rs1611615   |   -52.60| 2b              | NA                          | NA  | H3k4me2/H3k4me3/H3k79me2/H2az/H3k9ac/H3k27ac              |
| HLA-B    |     0| rs1265094   |   214.76| 4               | NA                          | NA  | NA                                                        |
| HLA-B    |     1| rs9380240   |    52.82| 4               | NA                          | NA  | H3k27me3                                                  |
| HLA-B    |     2| rs2523616   |     1.09| 3a              | POL2-4/MTA3/STAT5/BCL3/POL2 | NA  | H4k20me1/H3k4me1/H3k36me3/H3k4me2/H3k79me2                |
| HLA-C    |     0| rs146911342 |     0.00| 3a              | POL2-4                      | NA  | H4k20me1/H3k36me3/H3k4me3/H3k4me2/H3k79me2/H3k4me1/H3k9ac |
| HLA-C    |     1| rs9263875   |    65.61| 1f              | IKZF1/TCF12                 | NA  | H3k4me1/H3k4me3/H3k27ac/H3k4me2/H3k9ac/H3k79me2/H2az      |
| HLA-C    |     2| rs1793890   |    14.47| 6               | NA                          | NA  | NA                                                        |
| HLA-DPB1 |     0| rs9277449   |     0.00| 7               | NA                          | NA  | H3k36me3                                                  |
| HLA-DPB1 |     1| rs9296068   |   -55.01| 6               | NA                          | NA  | H3k4me1/H2az                                              |
| HLA-DPB1 |     2| rs688209    |   645.40| 4               | NA                          | NA  | NA                                                        |
| HLA-DQA1 |     0| rs9272316   |     0.00| 6               | NA                          | NA  | H3k27ac/H3k4me2/H3k9ac/H3k4me1/H3k4me3/H3k79me2           |
| HLA-DQA1 |     1| rs9271375   |    -8.89| 7               | NA                          | NA  | H3k9me3                                                   |
| HLA-DQA1 |     2| rs148834340 |  -106.13| 4               | NA                          | NA  | H3k79me2                                                  |
| HLA-DQB1 |     0| rs9274688   |    -0.74| 5               | POL2-4                      | NA  | H3k79me2/H3k4me2/H3k4me3/H3k27ac/H3k9ac/H3k4me1           |
| HLA-DQB1 |     1| rs3134978   |   -15.32| 7               | NA                          | NA  | NA                                                        |
| HLA-DRB1 |     0| rs3129759   |   -27.54| 6               | NA                          | NA  | NA                                                        |
| HLA-DRB1 |     1| rs9273375   |   -69.16| 4               | POL2-4                      | NA  | H3k36me3/H3k4me1/H3k27ac/H3k79me2                         |
| HLA-DRB1 |     2| rs34693360  |     0.00| 4               | PU1                         | DHS | H3k79me2/H3k9ac/H3k4me2/H3k4me3/H3k27ac/H3k36me3          |

### Reference transcriptome index

| locus    |  rank| rsid        |     dist| regulome\_score | tfbs                                                                                                              | dhs | histone\_marks                                       |
|:---------|-----:|:------------|--------:|:----------------|:------------------------------------------------------------------------------------------------------------------|:----|:-----------------------------------------------------|
| HLA-A    |     0| rs2523764   |   -91.42| 2b              | NA                                                                                                                | NA  | NA                                                   |
| HLA-A    |     1| rs2975046   |     1.49| 1f              | POL2-4/NFATC1/BCL3                                                                                                | NA  | H4k20me1/H3k36me3                                    |
| HLA-A    |     2| rs9259825   |   -14.19| 2b              | POL2-4/SIN3A/IRF4/POU2F2                                                                                          | DHS | H3k4me2/H3k4me3/H3k79me2/H2az/H3k9ac/H3k27ac         |
| HLA-B    |     0| rs2853926   |    58.60| 6               | NA                                                                                                                | NA  | NA                                                   |
| HLA-B    |     1| rs1265109   |   202.06| 6               | NA                                                                                                                | NA  | NA                                                   |
| HLA-B    |     2| rs36057735  |     1.73| 4               | POL2-4/NFATC1/MTA3/POL2/NFIC                                                                                      | NA  | H4k20me1/H3k4me1/H3k36me3                            |
| HLA-C    |     0| rs4947308   |   -62.38| 6               | NA                                                                                                                | NA  | H3k9me3                                              |
| HLA-C    |     1| rs2074491   |     0.00| 4               | POL2-4/NRSF/BCL11A/TCF3/TCF12/MXI1/EBF1/POU2F2/TBP/TAF1/POL2/NF-YA/CHD2/BHLHE40/CFOS/RFX5/NFKB/NF-YB/IRF4/PU1/SP1 | DHS | H3k4me3/H3k4me2/H3k79me2/H3k9ac/H3k27ac/H2az/H3k4me1 |
| HLA-C    |     2| rs1064627   |   537.99| 1f              | STAT3                                                                                                             | NA  | H3k36me3/H3k4me2                                     |
| HLA-C    |     3| rs9501587   |  -107.03| 6               | NA                                                                                                                | NA  | H3k27me3                                             |
| HLA-DPB1 |     0| rs9277538   |     0.07| 4               | NA                                                                                                                | NA  | NA                                                   |
| HLA-DPB1 |     1| rs9296068   |   -55.01| 6               | NA                                                                                                                | NA  | H3k4me1/H2az                                         |
| HLA-DPB1 |     2| rs688209    |   645.40| 4               | NA                                                                                                                | NA  | NA                                                   |
| HLA-DQA1 |     0| rs9270521   |   -36.43| 6               | NA                                                                                                                | NA  | NA                                                   |
| HLA-DQA1 |     1| rs3129758   |   -11.33| 5               | NA                                                                                                                | NA  | H3k9me3                                              |
| HLA-DQB1 |     0| rs1770      |     0.00| 1f              | POL2-4/CHD2/TBP/TAF1/EBF1/ELF1/TCF12/POL2                                                                         | DHS | H3k4me1/H3k27ac/H3k4me2/H3k36me3/H3k79me2            |
| HLA-DQB1 |     1| rs9274622   |     0.00| 1f              | POL2-4/IKZF1/AFT2                                                                                                 | NA  | H3k79me2/H3k4me2/H3k4me3/H3k27ac/H3k9ac/H2az/H3k4me1 |
| HLA-DQB1 |     2| rs28891461  |    -5.86| 7               | NA                                                                                                                | NA  | NA                                                   |
| HLA-DQB1 |     3| rs116102092 |    88.57| 6               | NA                                                                                                                | NA  | NA                                                   |
| HLA-DQB1 |     4| rs139547197 |  -408.00| 4               | POL2-4/MTA3                                                                                                       | NA  | H3k79me2/H3k4me3/H3k4me2/H3k9ac/H3k27ac/H2az         |
| HLA-DRB1 |     0| rs73729140  |     0.00| 6               | NA                                                                                                                | NA  | NA                                                   |
| HLA-DRB1 |     1| rs28383307  |   -29.13| 7               | NA                                                                                                                | NA  | H3k9me3                                              |

RTC between IMGT and Ref Transcriptome eQTLs
--------------------------------------------

Variants with RTC &gt; 0.95 likely mark the same biological signal.

| gene\_imgt | variant\_imgt |  rank\_imgt| gene\_pri | variant\_pri |  rank\_pri|  d\_prime|   rtc|
|:-----------|:--------------|-----------:|:----------|:-------------|----------:|---------:|-----:|
| HLA-A      | rs28780070    |           0| HLA-A     | rs2523764    |          0|      0.97|  0.64|
| HLA-A      | rs1611615     |           1| HLA-A     | rs2523764    |          0|      0.74|  0.98|
| HLA-C      | rs146911342   |           0| HLA-C     | rs4947308    |          0|      0.97|  0.95|
| HLA-C      | rs9263875     |           1| HLA-C     | rs2074491    |          1|      0.94|  0.97|
| HLA-C      | rs1793890     |           2| HLA-C     | rs2074491    |          1|      1.00|  0.87|
| HLA-B      | rs1265094     |           0| HLA-B     | rs1265109    |          1|      0.75|  0.99|
| HLA-B      | rs9380240     |           1| HLA-B     | rs2853926    |          0|      1.00|  0.96|
| HLA-B      | rs2523616     |           2| HLA-B     | rs2853926    |          0|      0.53|  0.47|
| HLA-DRB1   | rs3129759     |           0| HLA-DRB1  | rs73729140   |          0|      0.23|  0.24|
| HLA-DRB1   | rs9273375     |           1| HLA-DRB1  | rs28383307   |          1|      1.00|  0.29|
| HLA-DRB1   | rs34693360    |           2| HLA-DRB1  | rs28383307   |          1|      0.58|  0.98|
| HLA-DQA1   | rs9272316     |           0| HLA-DQA1  | rs3129758    |          1|      0.99|  0.77|
| HLA-DQA1   | rs9271375     |           1| HLA-DQA1  | rs3129758    |          1|      0.64|  0.79|
| HLA-DQA1   | rs148834340   |           2| HLA-DQA1  | rs3129758    |          1|      0.23|  0.60|
| HLA-DQB1   | rs9274688     |           0| HLA-DQB1  | rs1770       |          0|      1.00|  1.00|
| HLA-DQB1   | rs3134978     |           1| HLA-DQB1  | rs9274622    |          1|      1.00|  0.90|
| HLA-DPB1   | rs9277449     |           0| HLA-DPB1  | rs9277538    |          0|      0.99|  1.00|
| HLA-DPB1   | rs9296068     |           1| HLA-DPB1  | rs9296068    |          1|      1.00|  1.00|
| HLA-DPB1   | rs688209      |           2| HLA-DPB1  | rs688209     |          2|      1.00|  1.00|

RTC with previous eQTLs
-----------------------

| gene     |  rank| rsid        | qtl\_previous |  d\_prime|   rtc| study\_pval                                                 |
|:---------|-----:|:------------|:--------------|---------:|-----:|:------------------------------------------------------------|
| HLA-A    |     0| rs3823342   | rs3823342     |      1.00|  1.00| geuvadis\_exon (82.9)                                       |
| HLA-A    |     1| rs1655924   | rs1655924     |      1.00|  1.00| geuvadis\_gene (16.1)/geuvadis\_exon (38.5)                 |
| HLA-C    |     0| rs146911342 | rs146911342   |      1.00|  1.00| geuvadis\_gene (24.1)/geuvadis\_exon (26.2)                 |
| HLA-C    |     0| rs146911342 | rs41561715    |      1.00|  1.00| geuvadis\_gene (31.9)/geuvadis\_exon (37.2)                 |
| HLA-C    |     1| rs12199223  | rs12199223    |      1.00|  1.00| geuvadis\_exon (6.8)                                        |
| HLA-C    |     2| rs9266301   | rs9266301     |      1.00|  1.00| geuvadis\_exon (7.1)                                        |
| HLA-B    |     0| rs1265094   | rs1265094     |      1.00|  1.00| geuvadis\_exon (15.4)                                       |
| HLA-B    |     1| rs9264803   | rs9264803     |      1.00|  1.00| geuvadis\_exon (11.4)                                       |
| HLA-B    |     2| rs2308655   | rs2308655     |      1.00|  1.00| geuvadis\_exon (34.6)                                       |
| HLA-DRB1 |     0| rs3104412   | rs3104412     |      1.00|  1.00| geuvadis\_exon (17.7)                                       |
| HLA-DRB1 |     1| rs28483671  | rs28483671    |      1.00|  1.00| geuvadis\_exon (36.8)                                       |
| HLA-DQA1 |     0| rs75170544  | rs75170544    |      1.00|  1.00| geuvadis\_gene (28.4)/geuvadis\_exon (48.3)/gtex\_v7 (17.3) |
| HLA-DQA1 |     1| rs9271375   | rs9271375     |      1.00|  1.00| geuvadis\_exon (13.4)                                       |
| HLA-DQB1 |     0| rs1770      | rs1770        |      1.00|  1.00| geuvadis\_gene (43.4)/geuvadis\_exon (68.3)/gtex\_v7 (21.2) |
| HLA-DQB1 |     1| rs9272209   | rs9272209     |      1.00|  1.00| geuvadis\_gene (7)                                          |
| HLA-DPB1 |     0| rs9277449   | rs9277538     |      0.99|  1.00| geuvadis\_exon (57.9)                                       |
| HLA-DPB1 |     1| rs9296068   | rs34885310    |      0.96|  0.98| delaneau (3.7)                                              |
| HLA-DPB1 |     2| rs688209    | rs184379497   |      0.68|  0.95| NA (NA)                                                     |

Association with GWAS traits
----------------------------

| gene     |  rank| variant     | trait (GWAS variant)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
|:---------|-----:|:------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| HLA-A    |     0| rs28780070  | age at menarche (rs16896742); drug-induced liver injury (rs2523822); serum IgE measurement (rs2571391)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| HLA-A    |     1| rs1611615   | NA (NA)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| HLA-B    |     0| rs1265094   | membranous glomerulonephritis (rs1265159, rs3096697, rs3134945); neoplasm of mature B-cells (rs6457327); response to reverse transcriptase inhibitor (rs1265112)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-B    |     1| rs9380240   | Crohn's disease (rs9264942); HIV-1 infection (rs9264942); inflammatory bowel disease (rs9264942); Vitiligo (rs9468925)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| HLA-B    |     2| rs2523616   | serum IgE measurement (rs3130941)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| HLA-C    |     0| rs146911342 | body height (rs1265097); leukocyte count (rs2853946)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| HLA-C    |     1| rs9263875   | chronic obstructive pulmonary disease, surfactant protein D measurement (rs1265093, rs2074488); coronary heart disease (rs3869109)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| HLA-C    |     2| rs1793890   | chronic hepatitis B infection (rs2853953); cutaneous psoriasis measurement, psoriasis (rs10484554); membranous glomerulonephritis (rs1265159); monocyte count (rs3095254); nasopharyngeal neoplasm (rs2894207); neoplasm of mature B-cells (rs6457327); psoriasis (rs1265181, rs12191877, rs10484554); psoriasis vulgaris (rs4406273); psoriatic arthritis (rs13191343, rs12191877)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| HLA-DPB1 |     0| rs9277449   | NA (NA)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| HLA-DPB1 |     1| rs9296068   | NA (NA)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| HLA-DPB1 |     2| rs688209    | NA (NA)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| HLA-DQA1 |     0| rs9272316   | antibody measurement, Epstein-Barr virus infection (rs477515); antinuclear antibody measurement (rs2395185); Cystic fibrosis, lung disease severity measurement (rs9268947); Hodgkins lymphoma (rs2395185); inflammatory bowel disease (rs477515); lung carcinoma (rs2395185); lymphoma (rs9268853); NA (rs9268923); response to vaccine (rs477515); ulcerative colitis (rs9268853, rs9268923, rs2395185)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| HLA-DQA1 |     1| rs9271375   | asthma (rs7775228); chronic lymphocytic leukemia (rs9273363); Crohn's disease (rs9273363); Graves disease (rs2273017); inflammatory bowel disease (rs9273363); multiple sclerosis (rs3135338, rs3129871); multiple sclerosis, oligoclonal band measurement (rs3129871, rs3828840); Parkinson's disease (rs9275326); rheumatoid arthritis (rs6910071); seasonal allergic rhinitis (rs7775228); systemic scleroderma (rs3129763); type I diabetes mellitus (rs9268645); ulcerative colitis (rs9273363)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| HLA-DQA1 |     2| rs148834340 | alopecia areata (rs9275572); body height (rs6457620); chronic hepatitis B infection (rs2856718); Chronic Hepatitis C infection (rs4273729, rs9275572); cryoglobulinemia, Chronic Hepatitis C infection (rs9461776); Graves disease (rs6457617); hepatitis B infection (rs2856718); hepatitis C induced liver cirrhosis (rs910049); hepatocellular carcinoma (rs9275572); HPV seropositivity (rs9357152); multiple sclerosis (rs3129871); multiple sclerosis, IgG index (rs6457617); multiple sclerosis, oligoclonal band measurement (rs3129871, rs9275563); rheumatoid arthritis (rs6457617, rs6457620, rs9275572); systemic scleroderma (rs6457617); Vogt-Koyanagi-Harada disease (rs3021304)                                                                                                                                                                                                                                                                                                                                                                                                                           |
| HLA-DQB1 |     0| rs9274688   | lymphoma (rs2647045)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| HLA-DQB1 |     1| rs3134978   | antibody measurement, Epstein-Barr virus infection (rs2854275); asthma (rs3129943); autoimmune hepatits type 1 (rs2187668); autoimmune thyroid disease, type I diabetes mellitus (rs1270942, rs1980493); cancer, response to pazopanib, serum alanine aminotransferase measurement (rs1800625); celiac disease (rs2187668); chronic lymphocytic leukemia (rs674313, rs9273363); complement C4 measurement (rs2071278); Crohn's disease (rs9273363); cutaneous lupus erythematosus (rs9267531, rs2187668); hepatitis B infection (rs652888); inflammatory bowel disease (rs9273363); leprosy (rs602875); membranous glomerulonephritis (rs3130618, rs652888, rs389884, rs7775397, rs3129939, rs1980493, rs2187668, rs1480380); multiple sclerosis (rs2040406); protein measurement (rs2187668); rheumatoid arthritis (rs558702); schizophrenia (rs622076, rs3131296, rs9274623, rs73396800); systemic lupus erythematosus (rs3131379, rs1270942, rs1150757, rs1150754, rs1150753, rs2187668, rs3129716, rs3957147); systemic scleroderma (rs3129763); type I diabetes mellitus (rs2647044); ulcerative colitis (rs9273363) |
| HLA-DRB1 |     0| rs3129759   | Alzheimers disease (rs9271192); Crohn's disease (rs9271060); leprosy (rs9271100); seasonal allergic rhinitis, asthma (rs9273373); systemic lupus erythematosus (rs9271100); type I diabetes mellitus (rs9272346); ulcerative colitis (rs9271100, rs9271209)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| HLA-DRB1 |     1| rs9273375   | Crohn's disease (rs9271366); drug-induced liver injury (rs3129900); hepatitis C induced liver cirrhosis (rs3129860); multiple sclerosis (rs3129934, rs3135388, rs3129889, rs9271366); protein measurement (rs9271366); systemic lupus erythematosus (rs9270984, rs9273076); ulcerative colitis (rs9271366); ulcerative colitis, Crohn's disease (rs9271366)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| HLA-DRB1 |     2| rs34693360  | antinuclear antibody measurement (rs2395185); Hodgkins lymphoma (rs2395185); lung carcinoma (rs2395185); lymphoma (rs9268853); NA (rs9268923); ulcerative colitis (rs9268853, rs9268923, rs2395185)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |

HLA lineages
------------

<img src="./qtls/plots/lineage_and_effects.png" width="2100" />

### F-test: is there a difference between lineages?

#### traditional ANOVA

| locus    |   df|        F| p.value  |
|:---------|----:|--------:|:---------|
| HLA-A    |   12|   11.699| 1.15e-21 |
| HLA-B    |   14|    8.185| 2.77e-16 |
| HLA-C    |   11|   46.858| 8.18e-76 |
| HLA-DPB1 |    8|   12.557| 8.11e-17 |
| HLA-DQA1 |    4|  153.725| 1.10e-94 |
| HLA-DQB1 |    4|   43.112| 2.13e-32 |
| HLA-DRB1 |    9|   31.370| 4.85e-46 |

#### Welch ANOVA

| locus    |  num.df|  denom.df|        F| p.value  |
|:---------|-------:|---------:|--------:|:---------|
| HLA-A    |      12|   108.711|   11.465| 1.49e-14 |
| HLA-B    |      14|   147.301|    8.756| 1.21e-13 |
| HLA-C    |      11|   134.520|   32.954| 6.34e-33 |
| HLA-DPB1 |       8|    64.651|    8.369| 1.02e-07 |
| HLA-DQA1 |       4|   213.233|  118.651| 4.41e-53 |
| HLA-DQB1 |       4|   256.272|   48.426| 2.64e-30 |
| HLA-DRB1 |       9|   167.952|   44.434| 5.70e-40 |
