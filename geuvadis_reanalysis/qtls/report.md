Report
================

eQTLs
=====

**All analyses were carried out using European individuals only**

PCA of genotypes
----------------

![](./plots/genotype_pca.png)

Number of eGenes according to method of phenotype correction
------------------------------------------------------------

![](./plots/pca_vs_peer.png)

Distribution of eQTLs around the TSS
------------------------------------

![](./plots/qtls_landscape.png)

HLA lineages
------------

![](./plots/lineage_and_effects.png)

RTC
===

Geuvadis
--------

Some of our eQTLs seem to mark the same biological signal as the eQTLs found by Geuvadis (RTC &gt; 0.9)

| gene     | variant     |  rank| geuvadis\_gene    | geuvadis\_variant |    rtc|
|:---------|:------------|-----:|:------------------|:------------------|------:|
| HLA-A    | rs3823342   |     0| HLA-A             | rs114565353       |  0.358|
| HLA-A    | rs1655924   |     1| HLA-A             | rs114565353       |  0.998|
| HLA-B    | rs1265081   |     0| HLA-C             | rs115899777       |  0.982|
| HLA-B    | rs9265885   |     1| HLA-B             | rs137939159       |  0.298|
| HLA-C    | rs113169753 |     0| HLA-C             | rs115899777       |  0.977|
| HLA-C    | rs3134776   |     1| HLA-B             | rs137939159       |  0.915|
| HLA-C    | rs17408553  |     2| HLA-B             | rs137939159       |  0.351|
| HLA-DQA1 | rs9272316   |     0| HLA-DRB1          | rs116405062       |  0.970|
| HLA-DQA1 | rs9271375   |     1| HLA-DQA1/HLA-DQB1 | rs9274660         |  0.822|
| HLA-DQA1 | rs114968045 |     2| HLA-DQA1/HLA-DQB1 | rs9274660         |  0.771|
| HLA-DQB1 | rs1770      |     0| HLA-DQA1/HLA-DQB1 | rs9274660         |  0.996|
| HLA-DQB1 | rs4713570   |     1| HLA-DRB1          | rs116405062       |  0.708|
| HLA-DRB1 | rs9270698   |     0| HLA-DQA1/HLA-DQB1 | rs9274660         |  0.920|
| HLA-DRB1 | rs1048372   |     1| HLA-DQA1/HLA-DQB1 | rs9274660         |  0.723|
| HLA-DRB1 | rs9274063   |     2| HLA-DQA1/HLA-DQB1 | rs9274660         |  0.952|

Association with GWAS traits
----------------------------

| gene     |  rank| variant     | gwas\_variant |    rtc| trait                                            | studies                                                                                                                                                                                                                       |
|:---------|-----:|:------------|:--------------|------:|:-------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| HLA-A    |     0| rs3823342   | rs2517713     |  0.984| nasopharyngeal neoplasm                          | www.ncbi.nlm.nih.gov/pubmed/19664746                                                                                                                                                                                          |
| HLA-A    |     1| rs1655924   | rs2571391     |  0.943| serum IgE measurement                            | www.ncbi.nlm.nih.gov/pubmed/22075330                                                                                                                                                                                          |
| HLA-B    |     0| rs1265081   | rs1265112     |  0.997| response to reverse transcriptase inhibitor      | www.ncbi.nlm.nih.gov/pubmed/21810746                                                                                                                                                                                          |
| HLA-B    |     1| rs9265885   | rs2442719     |  0.999| HIV-1 infection, HIV viral set point measurement | www.ncbi.nlm.nih.gov/pubmed/26039976                                                                                                                                                                                          |
| HLA-C    |     0| rs113169753 | rs9263963     |  0.971| serum IgG glycosylation measurement              | www.ncbi.nlm.nih.gov/pubmed/23382691                                                                                                                                                                                          |
| HLA-C    |     1| rs3134776   | rs1265159     |  0.990| membranous glomerulonephritis                    | www.ncbi.nlm.nih.gov/pubmed/21323541                                                                                                                                                                                          |
| HLA-C    |     2| rs17408553  | rs2524079     |  0.998| lymphocyte count                                 | www.ncbi.nlm.nih.gov/pubmed/21738480                                                                                                                                                                                          |
| HLA-DQA1 |     0| rs9272316   | rs2395185     |  0.981| antinuclear antibody measurement                 | www.ncbi.nlm.nih.gov/pubmed/19122664 www.ncbi.nlm.nih.gov/pubmed/22286212 www.ncbi.nlm.nih.gov/pubmed/19915573 www.ncbi.nlm.nih.gov/pubmed/23143601 www.ncbi.nlm.nih.gov/pubmed/25186300 www.ncbi.nlm.nih.gov/pubmed/20228799 |
| HLA-DQA1 |     0| rs9272316   | rs2395185     |  0.981| Hodgkins lymphoma                                | www.ncbi.nlm.nih.gov/pubmed/19122664 www.ncbi.nlm.nih.gov/pubmed/22286212 www.ncbi.nlm.nih.gov/pubmed/19915573 www.ncbi.nlm.nih.gov/pubmed/23143601 www.ncbi.nlm.nih.gov/pubmed/25186300 www.ncbi.nlm.nih.gov/pubmed/20228799 |
| HLA-DQA1 |     0| rs9272316   | rs2395185     |  0.981| lung carcinoma                                   | www.ncbi.nlm.nih.gov/pubmed/19122664 www.ncbi.nlm.nih.gov/pubmed/22286212 www.ncbi.nlm.nih.gov/pubmed/19915573 www.ncbi.nlm.nih.gov/pubmed/23143601 www.ncbi.nlm.nih.gov/pubmed/25186300 www.ncbi.nlm.nih.gov/pubmed/20228799 |
| HLA-DQA1 |     0| rs9272316   | rs2395185     |  0.981| ulcerative colitis                               | www.ncbi.nlm.nih.gov/pubmed/19122664 www.ncbi.nlm.nih.gov/pubmed/22286212 www.ncbi.nlm.nih.gov/pubmed/19915573 www.ncbi.nlm.nih.gov/pubmed/23143601 www.ncbi.nlm.nih.gov/pubmed/25186300 www.ncbi.nlm.nih.gov/pubmed/20228799 |
| HLA-DQA1 |     0| rs9272316   | rs9268853     |  0.981| lymphoma                                         | www.ncbi.nlm.nih.gov/pubmed/21297633 www.ncbi.nlm.nih.gov/pubmed/23349640 www.ncbi.nlm.nih.gov/pubmed/23511034                                                                                                                |
| HLA-DQA1 |     0| rs9272316   | rs9268853     |  0.981| ulcerative colitis                               | www.ncbi.nlm.nih.gov/pubmed/21297633 www.ncbi.nlm.nih.gov/pubmed/23349640 www.ncbi.nlm.nih.gov/pubmed/23511034                                                                                                                |
| HLA-DQA1 |     0| rs9272316   | rs9268905     |  0.981| Cystic fibrosis                                  | www.ncbi.nlm.nih.gov/pubmed/21602797                                                                                                                                                                                          |
| HLA-DQA1 |     0| rs9272316   | rs9268923     |  0.981| ulcerative colitis                               | www.ncbi.nlm.nih.gov/pubmed/20228798 www.ncbi.nlm.nih.gov/pubmed/26819262                                                                                                                                                     |
| HLA-DQA1 |     0| rs9272316   | rs9268923     |  0.981| NA                                               | www.ncbi.nlm.nih.gov/pubmed/20228798 www.ncbi.nlm.nih.gov/pubmed/26819262                                                                                                                                                     |
| HLA-DQA1 |     1| rs9271375   | rs3135338     |  0.986| multiple sclerosis                               | www.ncbi.nlm.nih.gov/pubmed/20159113                                                                                                                                                                                          |
| HLA-DQA1 |     2| rs114968045 | rs9461776     |  0.997| cryoglobulinemia, Chronic Hepatitis C infection  | www.ncbi.nlm.nih.gov/pubmed/25030430                                                                                                                                                                                          |
| HLA-DQB1 |     0| rs1770      | rs6927022     |  0.998| ulcerative colitis                               | www.ncbi.nlm.nih.gov/pubmed/23128233                                                                                                                                                                                          |
| HLA-DQB1 |     1| rs4713570   | rs10484561    |  0.995| neoplasm of mature B-cells                       | www.ncbi.nlm.nih.gov/pubmed/20639881                                                                                                                                                                                          |
| HLA-DQB1 |     2| rs9275596   | rs9275596     |  1.000| IGA glomerulonephritis                           | www.ncbi.nlm.nih.gov/pubmed/21399633 www.ncbi.nlm.nih.gov/pubmed/25305756 www.ncbi.nlm.nih.gov/pubmed/25710614                                                                                                                |
| HLA-DQB1 |     2| rs9275596   | rs9275596     |  1.000| kidney disease                                   | www.ncbi.nlm.nih.gov/pubmed/21399633 www.ncbi.nlm.nih.gov/pubmed/25305756 www.ncbi.nlm.nih.gov/pubmed/25710614                                                                                                                |
| HLA-DQB1 |     2| rs9275596   | rs9275596     |  1.000| peanut allergy measurement                       | www.ncbi.nlm.nih.gov/pubmed/21399633 www.ncbi.nlm.nih.gov/pubmed/25305756 www.ncbi.nlm.nih.gov/pubmed/25710614                                                                                                                |
| HLA-DRB1 |     0| rs9270698   | rs9271060     |  0.992| Crohn's disease                                  | www.ncbi.nlm.nih.gov/pubmed/26891255                                                                                                                                                                                          |
| HLA-DRB1 |     0| rs9270698   | rs9271100     |  0.992| leprosy                                          | www.ncbi.nlm.nih.gov/pubmed/19838193 www.ncbi.nlm.nih.gov/pubmed/25642632 www.ncbi.nlm.nih.gov/pubmed/26192919                                                                                                                |
| HLA-DRB1 |     0| rs9270698   | rs9271100     |  0.992| systemic lupus erythematosus                     | www.ncbi.nlm.nih.gov/pubmed/19838193 www.ncbi.nlm.nih.gov/pubmed/25642632 www.ncbi.nlm.nih.gov/pubmed/26192919                                                                                                                |
| HLA-DRB1 |     0| rs9270698   | rs9271100     |  0.992| ulcerative colitis                               | www.ncbi.nlm.nih.gov/pubmed/19838193 www.ncbi.nlm.nih.gov/pubmed/25642632 www.ncbi.nlm.nih.gov/pubmed/26192919                                                                                                                |
| HLA-DRB1 |     0| rs9270698   | rs9271192     |  0.992| Alzheimers disease                               | www.ncbi.nlm.nih.gov/pubmed/24162737                                                                                                                                                                                          |
| HLA-DRB1 |     1| rs1048372   | rs9271348     |  0.999| rheumatoid arthritis                             | www.ncbi.nlm.nih.gov/pubmed/24532677                                                                                                                                                                                          |
| HLA-DRB1 |     2| rs9274063   | rs35597309    |  0.998| esophageal squamous cell carcinoma               | www.ncbi.nlm.nih.gov/pubmed/25129146                                                                                                                                                                                          |

Trans-eQTLs
-----------

-   Approximate pass as described on QTLtools website

| gene\_name | variant\_id |  variant\_chr|  variant\_pos|  log10\_nom\_pval|
|:-----------|:------------|-------------:|-------------:|-----------------:|
| HLA-DQA1   | rs3997880   |             2|     178431654|          22.02544|
| HLA-DQA1   | rs145339215 |             2|     178431521|          21.92733|
| HLA-DQA1   | rs3997878   |             2|     178431899|          21.20140|
| HLA-DQA1   | rs141354030 |             2|     178451007|          15.41356|
