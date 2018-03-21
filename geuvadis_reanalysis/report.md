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

### TPM

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

Numbers in parentheses are the CaVEMaN probabilities of being causal.

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
| HLA\_personalized | HLA-A    |                   0.41|                 0.75|                            -0.09|
| Reference         | HLA-A    |                   0.32|                -0.11|                             0.32|
| HLA\_personalized | HLA-B    |                   0.21|                 0.37|                             0.05|
| Reference         | HLA-B    |                  -0.03|                 0.03|                            -0.05|
| HLA\_personalized | HLA-C    |                   0.17|                 0.34|                            -0.07|
| Reference         | HLA-C    |                   0.15|                 0.34|                            -0.10|
| HLA\_personalized | HLA-DPB1 |                   0.62|                 0.94|                             0.14|
| Reference         | HLA-DPB1 |                   0.59|                 0.94|                             0.08|
| HLA\_personalized | HLA-DQA1 |                   0.59|                 0.66|                             0.32|
| Reference         | HLA-DQA1 |                  -0.57|                 0.62|                            -0.34|
| HLA\_personalized | HLA-DQB1 |                  -0.38|                 0.65|                            -0.02|
| Reference         | HLA-DQB1 |                  -0.46|                 0.37|                            -0.40|
| HLA\_personalized | HLA-DRB1 |                   0.58|                 0.41|                             0.45|
| Reference         | HLA-DRB1 |                   0.26|                 0.30|                             0.09|

Intersect with Encode Elements
------------------------------

### HLA-personalized index

| locus    |  rank| rsid        |       pos| tf                                        | dhs | chrom\_state | histone\_marks                                                            |
|:---------|-----:|:------------|---------:|:------------------------------------------|:----|:-------------|:--------------------------------------------------------------------------|
| HLA-A    |     0| rs3823342   |  29945290| POL2-4/POL2                               | NA  | TSS          | H3k27ac/H3k9ac/H3k4me2/H3k4me3/H3k79me2/H4k20me1/H3k36me3                 |
| HLA-A    |     1| rs1655924   |  29937342| NA                                        | NA  | R            | NA                                                                        |
| HLA-B    |     0| rs1265094   |  31139116| NA                                        | NA  | R            | NA                                                                        |
| HLA-B    |     1| rs9264803   |  31302674| NA                                        | NA  | R            | NA                                                                        |
| HLA-B    |     2| rs2308655   |  31354526| POL2-4                                    | DHS | TSS          | H4k20me1/H3k36me3/H3k4me2/H3k79me2/H3k4me3/H3k9ac/H3k27ac/H3k4me1/H3k9me3 |
| HLA-C    |     0| rs146911342 |  31270002| POL2-4                                    | NA  | TSS          | H4k20me1/H3k36me3/H3k4me3/H3k4me2/H3k79me2/H3k4me1/H3k9ac                 |
| HLA-C    |     1| rs12199223  |  31274954| NA                                        | NA  | R            | NA                                                                        |
| HLA-C    |     2| rs9266301   |  31361729| NA                                        | NA  | R            | NA                                                                        |
| HLA-DPB1 |     0| rs9277449   |  33084838| NA                                        | NA  | T            | H3k36me3                                                                  |
| HLA-DPB1 |     1| rs9296068   |  33020918| NA                                        | NA  | R            | H3k4me1/H2az                                                              |
| HLA-DPB1 |     2| rs688209    |  33732599| NA                                        | NA  | R            | NA                                                                        |
| HLA-DQA1 |     0| rs75170544  |  32545508| NA                                        | NA  | R            | NA                                                                        |
| HLA-DQA1 |     1| rs9271375   |  32619290| NA                                        | NA  | R            | H3k9me3                                                                   |
| HLA-DQB1 |     0| rs1770      |  32660056| POL2-4/CHD2/TBP/TAF1/EBF1/ELF1/TCF12/POL2 | DHS | E            | H3k4me1/H3k27ac/H3k4me2/H3k36me3/H3k79me2                                 |
| HLA-DQB1 |     1| rs9272209   |  32634200| RUNX3/NFIC/FOXM1                          | NA  | E            | H3k27ac/H3k4me2/H3k9ac/H3k4me1/H3k4me3                                    |
| HLA-DRB1 |     0| rs3104412   |  32618190| NA                                        | NA  | R            | H3k9me3                                                                   |
| HLA-DRB1 |     1| rs28483671  |  32559077| NA                                        | NA  | NA           | H3k79me2                                                                  |

### Reference transcriptome index

| locus    |  rank| rsid        |       pos| tf                      | dhs | chrom\_state | histone\_marks                                            |
|:---------|-----:|:------------|---------:|:------------------------|:----|:-------------|:----------------------------------------------------------|
| HLA-A    |     0| rs200093949 |  29932006| NA                      | NA  | T            | NA                                                        |
| HLA-A    |     1| rs2523764   |  29849840| NA                      | NA  | R            | NA                                                        |
| HLA-A    |     2| rs2734919   |  29918937| NA                      | NA  | R            | NA                                                        |
| HLA-B    |     0| rs3130949   |  31224090| NA                      | NA  | T            | NA                                                        |
| HLA-B    |     1| rs72443738  |  31305881| NA                      | NA  | R            | NA                                                        |
| HLA-B    |     2| rs2844623   |  31264766| NA                      | NA  | R            | NA                                                        |
| HLA-B    |     3| rs28380929  |  31351826| POL2-4/POL2/NFATC1      | NA  | T            | H4k20me1/H3k4me1/H3k36me3                                 |
| HLA-B    |     4| rs67565791  |  31093236| NA                      | NA  | R            | NA                                                        |
| HLA-C    |     0| rs146911342 |  31270002| POL2-4                  | NA  | TSS          | H4k20me1/H3k36me3/H3k4me3/H3k4me2/H3k79me2/H3k4me1/H3k9ac |
| HLA-C    |     1| rs1058067   |  31353904| POL2-4                  | NA  | TSS          | H4k20me1/H3k4me1/H3k36me3/H3k4me2/H3k79me2/H3k4me3/H3k9ac |
| HLA-C    |     2| rs9264185   |  31252205| NA                      | NA  | R            | NA                                                        |
| HLA-DPB1 |     0| rs9277449   |  33084838| NA                      | NA  | T            | H3k36me3                                                  |
| HLA-DPB1 |     1| rs9296068   |  33020918| NA                      | NA  | R            | H3k4me1/H2az                                              |
| HLA-DPB1 |     2| rs688209    |  33732599| NA                      | NA  | R            | NA                                                        |
| HLA-DQA1 |     0| rs28724008  |  32582773| NA                      | DHS | PF           | H3k36me3/H3k79me2                                         |
| HLA-DQA1 |     1| rs9275553   |  32709211| NA                      | NA  | R            | NA                                                        |
| HLA-DQB1 |     0| rs3830059   |  32666230| STAT5/SIN3A/POL2-4/WHIP | NA  | TSS          | H3k79me2/H3k4me2/H3k4me3/H3k27ac/H3k9ac/H2az              |
| HLA-DQB1 |     1| rs9273595   |  32661314| POL2-4                  | DHS | NA           | H3k4me1/H3k27ac/H3k4me2/H3k36me3/H3k79me2                 |
| HLA-DQB1 |     2| rs114969562 |  32626522| NA                      | NA  | NA           | H3k27ac/H3k4me1/H3k79me2                                  |
| HLA-DQB1 |     3| rs3134988   |  32672995| NA                      | NA  | T            | NA                                                        |
| HLA-DRB1 |     0| rs9269749   |  32580359| NA                      | NA  | T            | H3k36me3                                                  |
| HLA-DRB1 |     1| rs9271365   |  32619017| NA                      | NA  | R            | H3k9me3                                                   |
| HLA-DRB1 |     2| rs182670197 |  32534983| NA                      | NA  | R            | H3k9me3                                                   |

RTC between IMGT and Ref Transcriptome eQTLs
--------------------------------------------

Variants with RTC &gt; 0.95 likely mark the same biological signal.

| gene\_imgt | variant\_imgt |  rank\_imgt| gene\_pri | variant\_pri |  rank\_pri|  d\_prime|   rtc|
|:-----------|:--------------|-----------:|:----------|:-------------|----------:|---------:|-----:|
| HLA-A      | rs3823342     |           0| HLA-A     | rs200093949  |          0|      0.99|  0.92|
| HLA-A      | rs1655924     |           1| HLA-A     | rs2523764    |          1|      0.68|  0.94|
| HLA-C      | rs146911342   |           0| HLA-C     | rs146911342  |          0|      1.00|  1.00|
| HLA-C      | rs12199223    |           1| HLA-C     | rs9264185    |          2|      1.00|  0.97|
| HLA-C      | rs9266301     |           2| HLA-C     | rs1058067    |          1|      0.37|  0.98|
| HLA-B      | rs1265094     |           0| HLA-B     | rs3130949    |          0|      0.93|  0.98|
| HLA-B      | rs9264803     |           1| HLA-B     | rs2844623    |          2|      0.95|  0.97|
| HLA-B      | rs2308655     |           2| HLA-B     | rs3130949    |          0|      0.47|  0.91|
| HLA-DRB1   | rs3104412     |           0| HLA-DRB1  | rs9271365    |          1|      0.97|  1.00|
| HLA-DRB1   | rs28483671    |           1| HLA-DRB1  | rs9269749    |          0|      0.45|  0.71|
| HLA-DQA1   | rs75170544    |           0| HLA-DQA1  | rs28724008   |          0|      0.58|  0.97|
| HLA-DQA1   | rs9271375     |           1| HLA-DQA1  | rs9275553    |          1|      0.11|  0.22|
| HLA-DQB1   | rs1770        |           0| HLA-DQB1  | rs3830059    |          0|      0.95|  0.99|
| HLA-DQB1   | rs9272209     |           1| HLA-DQB1  | rs114969562  |          2|      0.76|  0.98|
| HLA-DPB1   | rs9277449     |           0| HLA-DPB1  | rs9277449    |          0|      1.00|  1.00|
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

| gene     | variant     |  rank| info                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
|:---------|:------------|-----:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| HLA-A    | rs3823342   |     0| age at menarche (rs16896742); nasopharyngeal neoplasm (rs2517713); multiple sclerosis (rs2523393); IGA glomerulonephritis (rs2523946); beta-2 microglobulin measurement (rs9260489); multiple sclerosis (rs9260489)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| HLA-A    | rs1655924   |     1| NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| HLA-B    | rs1265094   |     0| response to reverse transcriptase inhibitor (rs1265112); membranous glomerulonephritis (rs1265159); psoriasis (rs1265181); membranous glomerulonephritis (rs3096697); cutaneous lupus erythematosus (rs3130564); membranous glomerulonephritis (rs3130564); membranous glomerulonephritis (rs3134945); neoplasm of mature B-cells (rs6457327)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| HLA-B    | rs9264803   |     1| cutaneous psoriasis measurement, psoriasis (rs10484554); psoriasis (rs10484554); psoriasis (rs12191877); psoriatic arthritis (rs12191877); psoriatic arthritis (rs13191343); lymphocyte count (rs2524079); psoriasis vulgaris (rs4406273); HIV-1 infection (rs4418214); Crohn's disease (rs9264942); HIV-1 infection (rs9264942); inflammatory bowel disease (rs9264942)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| HLA-B    | rs2308655   |     2| membranous glomerulonephritis (rs1265159); psoriasis (rs1265181); body height (rs2247056); triglyceride measurement (rs2247056); body height (rs2256183); cervical carcinoma (rs2516448); diffuse large B-cell lymphoma (rs2523607); CD4:CD8 lymphocyte ratio (rs2524054); schizophrenia (rs2596500); rheumatoid arthritis (rs2596565); marginal zone B-cell lymphoma (rs2922994); myositis (rs3094013); immune system disease (rs3130544); membranous glomerulonephritis (rs3130544); cutaneous lupus erythematosus (rs3130564); membranous glomerulonephritis (rs3130564); cutaneous lupus erythematosus (rs3131060); myositis (rs3131619); membranous glomerulonephritis (rs3134792); psoriasis (rs3134792); body height (rs6457374); membranous glomerulonephritis (rs7750641); autoimmune thyroid disease, type I diabetes mellitus (rs886424)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| HLA-C    | rs146911342 |     0| body height (rs1265097)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| HLA-C    | rs12199223  |     1| cutaneous psoriasis measurement, psoriasis (rs10484554); psoriasis (rs10484554); psoriasis (rs12191877); psoriatic arthritis (rs12191877); response to reverse transcriptase inhibitor (rs1265112); membranous glomerulonephritis (rs1265159); psoriasis (rs1265181); psoriatic arthritis (rs13191343); chronic hepatitis B infection (rs2853953); nasopharyngeal neoplasm (rs2894207); psoriasis vulgaris (rs4406273)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| HLA-C    | rs9266301   |     2| chronic obstructive pulmonary disease, surfactant protein D measurement (rs2074488); chronic obstructive pulmonary disease, surfactant protein D measurement (rs9266629)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| HLA-DPB1 | rs9277449   |     0| NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| HLA-DPB1 | rs9296068   |     1| lymphoma (rs4530903)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| HLA-DPB1 | rs688209    |     2| NA                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| HLA-DQA1 | rs75170544  |     0| antinuclear antibody measurement (rs2395185); Hodgkins lymphoma (rs2395185); lung carcinoma (rs2395185); ulcerative colitis (rs2395185); antibody measurement, Epstein-Barr virus infection (rs477515); inflammatory bowel disease (rs477515); response to vaccine (rs477515); Vogt-Koyanagi-Harada disease (rs9268838); lymphoma (rs9268853); ulcerative colitis (rs9268853); NA (rs9268923); ulcerative colitis (rs9268923); Cystic fibrosis, lung disease severity measurement (rs9268947)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| HLA-DQA1 | rs9271375   |     1| Graves disease (rs2273017); systemic scleroderma (rs3129763); multiple sclerosis (rs3135338); multiple sclerosis, oligoclonal band measurement (rs3828840); type I diabetes mellitus (rs9268645)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| HLA-DQB1 | rs1770      |     0| lymphoma (rs2647045); ulcerative colitis (rs6927022); ulcerative colitis (rs9268877); type I diabetes mellitus (rs9272346); chronic lymphocytic leukemia (rs9273363); Crohn's disease (rs9273363); inflammatory bowel disease (rs9273363); ulcerative colitis (rs9273363); seasonal allergic rhinitis, asthma (rs9273373)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| HLA-DQB1 | rs9272209   |     1| systemic lupus erythematosus (rs1150753); systemic lupus erythematosus (rs1150757); autoimmune thyroid disease, type I diabetes mellitus (rs1270942); systemic lupus erythematosus (rs1270942); membranous glomerulonephritis (rs1480380); autoimmune thyroid disease, type I diabetes mellitus (rs1980493); membranous glomerulonephritis (rs1980493); multiple sclerosis (rs2040406); autoimmune hepatits type 1 (rs2187668); celiac disease (rs2187668); cutaneous lupus erythematosus (rs2187668); membranous glomerulonephritis (rs2187668); protein measurement (rs2187668); systemic lupus erythematosus (rs2187668); type I diabetes mellitus (rs2647044); antibody measurement, Epstein-Barr virus infection (rs2854275); systemic lupus erythematosus (rs3129716); systemic scleroderma (rs3129763); membranous glomerulonephritis (rs3130618); schizophrenia (rs3131296); systemic lupus erythematosus (rs3131379); hepatitis C induced liver cirrhosis (rs3135363); response to vaccine (rs3135363); membranous glomerulonephritis (rs389884); systemic lupus erythematosus (rs3957147); rheumatoid arthritis (rs558702); leprosy (rs602875); schizophrenia (rs622076); chronic lymphocytic leukemia (rs674313); schizophrenia (rs73396800); membranous glomerulonephritis (rs7775397); cutaneous lupus erythematosus (rs9267531); schizophrenia (rs9274623) |
| HLA-DRB1 | rs3104412   |     0| type I diabetes mellitus (rs9272346); seasonal allergic rhinitis, asthma (rs9273373)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     |
| HLA-DRB1 | rs28483671  |     1| hepatitis C induced liver cirrhosis (rs3129860); multiple sclerosis (rs3129889); multiple sclerosis (rs3135388); systemic lupus erythematosus (rs9270984); Crohn's disease (rs9271366); multiple sclerosis (rs9271366); protein measurement (rs9271366); ulcerative colitis (rs9271366); ulcerative colitis, Crohn's disease (rs9271366)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |

HLA lineages
------------

<img src="./qtls/plots/lineage_and_effects.png" width="2100" />

### F-test: is there a difference between lineages?

#### traditional ANOVA

| locus    |   df|        F| p.value  |
|:---------|----:|--------:|:---------|
| HLA-A    |   12|   12.138| 1.52e-22 |
| HLA-B    |   14|    8.381| 9.59e-17 |
| HLA-C    |   11|   56.132| 3.66e-88 |
| HLA-DPB1 |    8|   13.568| 3.05e-18 |
| HLA-DQA1 |    4|  140.118| 3.07e-88 |
| HLA-DQB1 |    4|   47.233| 3.30e-35 |
| HLA-DRB1 |    9|   34.090| 1.16e-49 |

#### Welch ANOVA

| locus    |  num.df|  denom.df|        F| p.value  |
|:---------|-------:|---------:|--------:|:---------|
| HLA-A    |      12|   111.205|   12.112| 2.34e-15 |
| HLA-B    |      14|   147.755|    9.538| 8.23e-15 |
| HLA-C    |      11|   139.582|   35.195| 6.45e-35 |
| HLA-DPB1 |       8|    63.895|    9.484| 1.55e-08 |
| HLA-DQA1 |       4|   215.519|  107.624| 3.07e-50 |
| HLA-DQB1 |       4|   256.090|   53.511| 9.80e-33 |
| HLA-DRB1 |       9|   167.656|   48.290| 4.59e-42 |
