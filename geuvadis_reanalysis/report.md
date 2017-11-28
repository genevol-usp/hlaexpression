Report
================

Typing accuracies
=================

\*Concordance: the proportion of the called alleles that are concordant with the Gourraud et al (2014) typings

| locus |   star|  kallisto|
|:------|------:|---------:|
| A     |  97.55|     97.06|
| B     |  98.04|     97.55|
| C     |  96.90|     95.92|
| DQB1  |  98.37|     98.53|
| DRB1  |  99.51|     99.51|

Expression estimates
====================

![](./expression/plots/expression_boxplot.png)

kallisto vs STAR-Salmon
-----------------------

### TPM

![](./expression/plots/star_vs_kallisto_TPM.png)

### PCA-corrected

![](./expression/plots/star_vs_kallisto_PCA.png)

HLA diversity vs reference chromosomes only
-------------------------------------------

### TPM

#### STAR

![](./expression/plots/star_imgt_vs_pri_TPM.png)

#### kallisto

![](./expression/plots/kallisto_imgt_vs_pri_TPM.png)

### PCA-corrected

#### STAR

![](./expression/plots/star_imgt_vs_pri_PCA.png)

#### kallisto

![](./expression/plots/kallisto_imgt_vs_pri_PCA.png)

Distribution of TPM values
--------------------------

![](./expression/plots/tpm_distributions.png)

ASE
---

### ASE by number of genotyping errors

\*Each point represents a heterozygous genotype in the intersect with Gourraud data.

![](./expression/plots/ase.png)

### ASE distribution

![](./expression/plots/ase_histogram.png)

Correlation of expression
-------------------------

![](./expression/plots/correlation_decrease.png)

### Among the HLA genes

![](./expression/plots/hlacorrelations.png)

### Between Class II genes and CIITA

![](./expression/plots/trans_activ_corrs.png)

### Between pairs of HLA genes on the same vs on different haplotypes

#### HLA-A vs HLA-B

![](./expression/plots/a_vs_b.png)

#### HLA-A vs HLA-C

![](./expression/plots/a_vs_c.png)

#### HLA-B vs HLA-C

![](./expression/plots/b_vs_c.png)

#### HLA-DQA1 vs HLA-DQB1

![](./expression/plots/dqa_vs_dqb.png)

#### HLA-DQA1 vs HLA-DRB1

![](./expression/plots/dqa_vs_drb.png)

eQTLs
=====

**All analyses were carried out using European individuals only**

PCA of genotypes
----------------

![](./qtls/plots/genotype_pca.png)

Number of eGenes according to index
-----------------------------------

![](./qtls/plots/n_of_egenes.png)

Distribution of eQTLs around the TSS
------------------------------------

### IMGT index

![](./qtls/plots/qtls_landscape_imgt.png)

### Reference transcriptome

![](./qtls/plots/qtls_landscape_pri.png)

RTC between IMGT and Ref Transcriptome eQTLs
--------------------------------------------

Variants with RTC &gt; 0.9 likely mark the same biological signal.

| gene\_imgt | variant\_imgt |  rank\_imgt| gene\_pri | variant\_pri |  rank\_pri|  d\_prime|   rtc|
|:-----------|:--------------|-----------:|:----------|:-------------|----------:|---------:|-----:|
| HLA-A      | rs7383537     |           2| HLA-A     | rs28724990   |          0|      0.71|  0.92|
| HLA-B      | rs1265094     |           0| HLA-B     | rs1265159    |          0|      0.99|  0.99|
| HLA-B      | rs1265094     |           0| HLA-C     | rs9468965    |          0|      0.91|  0.90|
| HLA-B      | rs1265094     |           0| HLA-C     | rs2394941    |          2|      0.91|  0.96|
| HLA-B      | rs9264803     |           1| HLA-B     | rs2844623    |          1|      0.95|  0.98|
| HLA-C      | rs41561715    |           0| HLA-C     | rs9468965    |          0|      0.96|  0.99|
| HLA-C      | rs2394941     |           1| HLA-B     | rs1265159    |          0|      0.99|  0.99|
| HLA-C      | rs2394941     |           1| HLA-C     | rs2394941    |          2|      1.00|  1.00|
| HLA-C      | rs2249741     |           2| HLA-B     | rs2844623    |          1|      1.00|  0.94|
| HLA-C      | rs2249741     |           2| HLA-C     | rs2074491    |          1|      1.00|  0.98|
| HLA-DPB1   | rs9277538     |           0| HLA-DPB1  | rs9277538    |          0|      1.00|  1.00|
| HLA-DPB1   | rs9296068     |           1| HLA-DPB1  | rs9296068    |          1|      1.00|  1.00|
| HLA-DQA1   | rs9272316     |           0| HLA-DQA1  | rs144150162  |          0|      0.96|  0.98|
| HLA-DQA1   | rs9272316     |           0| HLA-DQB1  | rs28746846   |          0|      0.82|  0.90|
| HLA-DQA1   | rs9271375     |           1| HLA-DRB1  | rs9271365    |          1|      0.92|  1.00|
| HLA-DQA1   | rs9271375     |           1| HLA-DQB1  | rs114969562  |          2|      0.72|  0.93|
| HLA-DQB1   | rs1770        |           0| HLA-DQB1  | rs28746846   |          0|      0.95|  0.99|
| HLA-DQB1   | rs9272209     |           1| HLA-DQB1  | rs9273595    |          1|      0.92|  0.97|
| HLA-DQB1   | rs9272209     |           1| HLA-DQB1  | rs114969562  |          2|      0.76|  0.98|
| HLA-DQB1   | rs9269386     |           2| HLA-DRB1  | rs182670197  |          2|      0.80|  0.96|
| HLA-DRB1   | rs3104412     |           0| HLA-DRB1  | rs9271365    |          1|      0.97|  1.00|
| HLA-DRB1   | rs3104412     |           0| HLA-DQB1  | rs3134988    |          3|      0.77|  0.94|

Comparison with previous eQTLs
------------------------------

### eQTLs in regulomeDB

| gene     |  rank| rsid       | score | qtl                         |
|:---------|-----:|:-----------|:------|:----------------------------|
| HLA-A    |     0| rs16896724 | 5     | NA                          |
| HLA-A    |     1| rs2517827  | 7     | NA                          |
| HLA-A    |     2| rs7383537  | 6     | NA                          |
| HLA-B    |     0| rs1265094  | 4     | NA                          |
| HLA-B    |     1| rs9264803  | 6     | NA                          |
| HLA-C    |     0| rs41561715 | 4     | NA                          |
| HLA-C    |     1| rs2394941  | 5     | NA                          |
| HLA-C    |     2| rs2249741  | 1f    | HLA-C\_eQTL\_Lymphoblastoid |
| HLA-DPB1 |     0| rs9277538  | 4     | NA                          |
| HLA-DPB1 |     1| rs9296068  | 6     | HLA-DPB1\_eQTL\_Monocytes   |
| HLA-DPB1 |     2| rs688209   | 4     | NA                          |
| HLA-DQA1 |     0| rs9272316  | 6     | NA                          |
| HLA-DQA1 |     1| rs9271375  | 7     | NA                          |
| HLA-DQB1 |     0| rs1770     | 1f    | dsQTL\_Lymphoblastoid       |
| HLA-DQB1 |     1| rs9272209  | 3a    | NA                          |
| HLA-DQB1 |     2| rs9269386  | 7     | NA                          |
| HLA-DRB1 |     0| rs3104412  | 6     | NA                          |
| HLA-DRB1 |     1| rs9270700  | 5     | NA                          |

### eQTLs in HaploReg

| gene     |  rank| rsID       | study           | tissue                                    |      pvalue|
|:---------|-----:|:-----------|:----------------|:------------------------------------------|-----------:|
| HLA-A    |     0| rs16896724 | GTEx2015\_v6    | Testis                                    |    6.970055|
| HLA-A    |     0| rs16896724 | GTEx2015\_v6    | Artery\_Tibial                            |    6.721690|
| HLA-A    |     0| rs16896724 | GTEx2015\_v6    | Esophagus\_Muscularis                     |    5.074220|
| HLA-A    |     1| rs2517827  | GTEx2015\_v6    | Artery\_Tibial                            |    7.596338|
| HLA-A    |     1| rs2517827  | GTEx2015\_v6    | Artery\_Aorta                             |    7.230838|
| HLA-A    |     1| rs2517827  | GTEx2015\_v6    | Cells\_Transformed\_fibroblasts           |    6.614386|
| HLA-A    |     1| rs2517827  | GTEx2015\_v6    | Nerve\_Tibial                             |    5.194803|
| HLA-A    |     2| rs7383537  | NA              | NA                                        |          NA|
| HLA-B    |     0| rs1265094  | Lappalainen2013 | Lymphoblastoid\_EUR\_exonlevel            |   15.356583|
| HLA-B    |     1| rs9264803  | NA              | NA                                        |          NA|
| HLA-C    |     0| rs41561715 | NA              | NA                                        |          NA|
| HLA-C    |     1| rs2394941  | GTEx2015\_v6    | Whole\_Blood                              |    6.510757|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Muscle\_Skeletal                          |   23.816378|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Cells\_Transformed\_fibroblasts           |   22.516603|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Artery\_Tibial                            |   20.557065|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Skin\_Sun\_Exposed\_Lower\_leg            |   19.716441|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Esophagus\_Mucosa                         |   14.186417|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Adipose\_Subcutaneous                     |   14.126082|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Esophagus\_Muscularis                     |   12.767794|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Nerve\_Tibial                             |   12.107584|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Lung                                      |   11.929606|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Heart\_Left\_Ventricle                    |   11.508504|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Thyroid                                   |   10.864279|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Adipose\_Visceral\_Omentum                |   10.352884|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Skin\_Not\_Sun\_Exposed\_Suprapubic       |   10.163933|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Breast\_Mammary\_Tissue                   |    9.649241|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Artery\_Aorta                             |    9.008754|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Brain\_Cerebellum                         |    8.752610|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Whole\_Blood                              |    8.719715|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Colon\_Transverse                         |    8.463598|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Stomach                                   |    8.358130|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Brain\_Hypothalamus                       |    8.044904|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Liver                                     |    7.497827|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Testis                                    |    7.097510|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Pituitary                                 |    6.961919|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Colon\_Sigmoid                            |    6.703395|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Pancreas                                  |    6.524238|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Brain\_Caudate\_basal\_ganglia            |    6.490060|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Cells\_EBV-transformed\_lymphocytes       |    6.148427|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Esophagus\_Gastroesophageal\_Junction     |    5.960101|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Adrenal\_Gland                            |    5.415203|
| HLA-C    |     2| rs2249741  | GTEx2015\_v6    | Brain\_Cerebellar\_Hemisphere             |    5.338115|
| HLA-DPB1 |     0| rs9277538  | Westra2013      | Whole\_Blood                              |  197.008179|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Skin\_Sun\_Exposed\_Lower\_leg            |   60.914401|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Lung                                      |   60.055495|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Thyroid                                   |   59.088221|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Artery\_Tibial                            |   56.734184|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Nerve\_Tibial                             |   50.606605|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Muscle\_Skeletal                          |   48.583396|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Esophagus\_Muscularis                     |   43.947374|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Whole\_Blood                              |   43.113227|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Artery\_Aorta                             |   40.465541|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Adipose\_Subcutaneous                     |   37.181776|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Skin\_Not\_Sun\_Exposed\_Suprapubic       |   32.334648|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Brain\_Cerebellum                         |   30.608913|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Esophagus\_Mucosa                         |   27.653941|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Breast\_Mammary\_Tissue                   |   26.428902|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Heart\_Atrial\_Appendage                  |   23.713445|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Brain\_Cerebellar\_Hemisphere             |   23.584761|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Colon\_Transverse                         |   23.552799|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Heart\_Left\_Ventricle                    |   22.574776|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Esophagus\_Gastroesophageal\_Junction     |   20.831067|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Testis                                    |   20.614541|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Artery\_Coronary                          |   20.006289|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Pancreas                                  |   19.875123|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Adipose\_Visceral\_Omentum                |   19.161209|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Stomach                                   |   18.064263|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Spleen                                    |   17.271497|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Cells\_Transformed\_fibroblasts           |   15.814070|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Liver                                     |   11.316593|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Brain\_Hypothalamus                       |   11.295993|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Prostate                                  |   11.170276|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Colon\_Sigmoid                            |   10.902804|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Cells\_EBV-transformed\_lymphocytes       |   10.395900|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Brain\_Hippocampus                        |   10.386991|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Small\_Intestine\_Terminal\_Ileum         |    9.910858|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Brain\_Cortex                             |    8.792221|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Brain\_Nucleus\_accumbens\_basal\_ganglia |    8.443727|
| HLA-DPB1 |     0| rs9277538  | GTEx2015\_v6    | Brain\_Putamen\_basal\_ganglia            |    6.384340|
| HLA-DPB1 |     1| rs9296068  | Westra2013      | Whole\_Blood                              |  105.585393|
| HLA-DPB1 |     2| rs688209   | NA              | NA                                        |          NA|
| HLA-DQA1 |     0| rs9272316  | NA              | NA                                        |          NA|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Muscle\_Skeletal                          |   22.685239|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Skin\_Sun\_Exposed\_Lower\_leg            |   21.250407|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Whole\_Blood                              |   19.123646|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Adipose\_Subcutaneous                     |   16.239829|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Thyroid                                   |   13.239817|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Artery\_Tibial                            |   12.612294|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Adipose\_Visceral\_Omentum                |   11.676662|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Esophagus\_Mucosa                         |   11.184226|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Skin\_Not\_Sun\_Exposed\_Suprapubic       |   11.056696|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Breast\_Mammary\_Tissue                   |   10.589586|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Lung                                      |    9.999028|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Testis                                    |    9.669677|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Esophagus\_Muscularis                     |    9.531619|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Stomach                                   |    9.509276|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Nerve\_Tibial                             |    9.084250|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Esophagus\_Gastroesophageal\_Junction     |    8.904301|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Heart\_Atrial\_Appendage                  |    8.895583|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Artery\_Aorta                             |    8.814401|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Colon\_Transverse                         |    8.038032|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Pituitary                                 |    7.233891|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Colon\_Sigmoid                            |    6.512684|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Heart\_Left\_Ventricle                    |    6.444339|
| HLA-DQA1 |     1| rs9271375  | GTEx2015\_v6    | Adrenal\_Gland                            |    5.820496|
| HLA-DQB1 |     0| rs1770     | Lappalainen2013 | Lymphoblastoid\_EUR\_exonlevel            |   68.288868|
| HLA-DQB1 |     0| rs1770     | Lappalainen2013 | Lymphoblastoid\_EUR\_genelevel            |   43.432284|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Skin\_Sun\_Exposed\_Lower\_leg            |   34.498779|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Whole\_Blood                              |   31.304094|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Muscle\_Skeletal                          |   28.954606|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Nerve\_Tibial                             |   28.028577|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Adipose\_Subcutaneous                     |   25.596807|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Lung                                      |   25.364010|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Thyroid                                   |   25.003063|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Esophagus\_Mucosa                         |   24.650149|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Skin\_Not\_Sun\_Exposed\_Suprapubic       |   22.520603|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Artery\_Tibial                            |   22.468819|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Heart\_Left\_Ventricle                    |   19.023305|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Stomach                                   |   18.008471|
| HLA-DQB1 |     0| rs1770     | Lappalainen2013 | Lymphoblastoid\_YRI\_exonlevel            |   16.284500|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Breast\_Mammary\_Tissue                   |   16.187618|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Esophagus\_Muscularis                     |   15.281161|
| HLA-DQB1 |     0| rs1770     | Lappalainen2013 | Lymphoblastoid\_YRI\_genelevel            |   15.001735|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Adipose\_Visceral\_Omentum                |   14.436470|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Colon\_Transverse                         |   13.656355|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Artery\_Aorta                             |   13.526991|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Heart\_Atrial\_Appendage                  |   12.665867|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Cells\_EBV-transformed\_lymphocytes       |   12.537449|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Pancreas                                  |   11.924593|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Esophagus\_Gastroesophageal\_Junction     |   11.575731|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Liver                                     |   10.861057|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Brain\_Frontal\_Cortex\_BA9               |   10.810076|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Artery\_Coronary                          |   10.801265|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Spleen                                    |    9.632579|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Brain\_Caudate\_basal\_ganglia            |    9.274431|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Colon\_Sigmoid                            |    9.260475|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Brain\_Anterior\_cingulate\_cortex\_BA24  |    9.235356|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Vagina                                    |    9.228642|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Testis                                    |    9.157207|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Small\_Intestine\_Terminal\_Ileum         |    8.479718|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Brain\_Hypothalamus                       |    7.969335|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Brain\_Cerebellar\_Hemisphere             |    7.931159|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Adrenal\_Gland                            |    7.767512|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Uterus                                    |    7.651149|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Pituitary                                 |    7.455304|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Ovary                                     |    7.169072|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Brain\_Putamen\_basal\_ganglia            |    6.676998|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Brain\_Cortex                             |    6.582086|
| HLA-DQB1 |     0| rs1770     | GTEx2015\_v6    | Brain\_Nucleus\_accumbens\_basal\_ganglia |    6.021652|
| HLA-DQB1 |     1| rs9272209  | GTEx2015\_v6    | Esophagus\_Mucosa                         |    5.467033|
| HLA-DQB1 |     1| rs9272209  | GTEx2015\_v6    | Muscle\_Skeletal                          |    5.375633|
| HLA-DQB1 |     2| rs9269386  | NA              | NA                                        |          NA|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Muscle\_Skeletal                          |   15.182436|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Skin\_Sun\_Exposed\_Lower\_leg            |   14.477008|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Adipose\_Subcutaneous                     |   13.044064|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Whole\_Blood                              |   12.411939|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Nerve\_Tibial                             |   10.647275|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Esophagus\_Mucosa                         |   10.404950|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Thyroid                                   |    8.838950|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Skin\_Not\_Sun\_Exposed\_Suprapubic       |    8.485531|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Lung                                      |    8.090744|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Adipose\_Visceral\_Omentum                |    7.694177|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Artery\_Tibial                            |    7.605068|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Artery\_Coronary                          |    7.140416|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Artery\_Aorta                             |    6.600229|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Esophagus\_Gastroesophageal\_Junction     |    6.265296|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Esophagus\_Muscularis                     |    6.092468|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Colon\_Transverse                         |    6.006194|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Testis                                    |    5.830163|
| HLA-DRB1 |     0| rs3104412  | GTEx2015\_v6    | Breast\_Mammary\_Tissue                   |    5.522108|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Whole\_Blood                              |   73.379518|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Adipose\_Subcutaneous                     |   54.972231|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Muscle\_Skeletal                          |   54.083624|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Artery\_Tibial                            |   50.214338|
| HLA-DRB1 |     1| rs9270700  | Lappalainen2013 | Lymphoblastoid\_EUR\_exonlevel            |   49.846491|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Skin\_Sun\_Exposed\_Lower\_leg            |   49.625861|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Lung                                      |   47.215927|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Thyroid                                   |   46.692452|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Nerve\_Tibial                             |   44.859765|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Esophagus\_Mucosa                         |   40.132813|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Artery\_Aorta                             |   30.486852|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Adipose\_Visceral\_Omentum                |   30.433806|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Heart\_Left\_Ventricle                    |   28.857675|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Skin\_Not\_Sun\_Exposed\_Suprapubic       |   27.486600|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Esophagus\_Muscularis                     |   27.414735|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Colon\_Transverse                         |   26.540620|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Testis                                    |   26.201386|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Breast\_Mammary\_Tissue                   |   25.767346|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Pancreas                                  |   25.388503|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Heart\_Atrial\_Appendage                  |   25.256971|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Stomach                                   |   20.235734|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Colon\_Sigmoid                            |   20.127755|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Esophagus\_Gastroesophageal\_Junction     |   18.118585|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Brain\_Cerebellum                         |   16.138904|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Cells\_Transformed\_fibroblasts           |   15.946338|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Brain\_Frontal\_Cortex\_BA9               |   15.801406|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Cells\_EBV-transformed\_lymphocytes       |   15.729568|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Artery\_Coronary                          |   15.502442|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Brain\_Hippocampus                        |   14.784203|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Brain\_Caudate\_basal\_ganglia            |   14.079291|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Brain\_Cortex                             |   14.039807|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Uterus                                    |   13.247042|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Brain\_Nucleus\_accumbens\_basal\_ganglia |   13.168902|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Liver                                     |   12.579252|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Brain\_Anterior\_cingulate\_cortex\_BA24  |   12.303637|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Brain\_Hypothalamus                       |   12.206855|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Prostate                                  |   12.187087|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Brain\_Putamen\_basal\_ganglia            |   11.391908|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Pituitary                                 |   11.379027|
| HLA-DRB1 |     1| rs9270700  | Lappalainen2013 | Lymphoblastoid\_YRI\_exonlevel            |   11.298293|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Brain\_Cerebellar\_Hemisphere             |   11.097188|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Vagina                                    |   10.523146|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Small\_Intestine\_Terminal\_Ileum         |    9.759650|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Ovary                                     |    9.667626|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Adrenal\_Gland                            |    9.041808|
| HLA-DRB1 |     1| rs9270700  | GTEx2015\_v6    | Spleen                                    |    8.960945|

### eQTLs in GRASP (via HaploReg)

| gene     |  rank| rsid      | phenotype                                                                                  |     pvalue| study                                |
|:---------|-----:|:----------|:-------------------------------------------------------------------------------------------|----------:|:-------------------------------------|
| HLA-B    |     0| rs1265094 | Gene expression of HCG22 \[probe 1560767\_at\] in lymphoblastoid cell lines                |  16.823909| www.ncbi.nlm.nih.gov/pubmed/17873877 |
| HLA-B    |     0| rs1265094 | Gene expression of HCG22 in lymphoblastoid cell lines                                      |  15.792999| www.ncbi.nlm.nih.gov/pubmed/17873877 |
| HLA-B    |     0| rs1265094 | Gene expression of TCF19 \[probe 223274\_at\] in lymphoblastoid cell lines                 |  10.522879| www.ncbi.nlm.nih.gov/pubmed/17873877 |
| HLA-B    |     0| rs1265094 | Gene expression of HCG27 \[probe 1559050\_at\] in lymphoblastoid cell lines                |  10.161151| www.ncbi.nlm.nih.gov/pubmed/17873877 |
| HLA-B    |     0| rs1265094 | Gene expression of TCF19 in lymphoblastoid cell lines                                      |   9.591999| www.ncbi.nlm.nih.gov/pubmed/17873877 |
| HLA-B    |     0| rs1265094 | Gene expression of HCG27 in lymphoblastoid cell lines                                      |   9.241000| www.ncbi.nlm.nih.gov/pubmed/17873877 |
| HLA-B    |     0| rs1265094 | Gene expression of HLA-C \[probe 216526\_x\_at\] in lymphoblastoid cell lines              |   9.086186| www.ncbi.nlm.nih.gov/pubmed/17873877 |
| HLA-B    |     0| rs1265094 | Gene expression of CCHCR1 in normal prepouch ileum                                         |   4.823275| www.ncbi.nlm.nih.gov/pubmed/23474282 |
| HLA-C    |     2| rs2249741 | Gene expression of HLA-C in CEU-CHB-JPT-YRI lymphoblastoid cell lines                      |  17.588901| www.ncbi.nlm.nih.gov/pubmed/17873874 |
| HLA-C    |     2| rs2249741 | Gene expression of hmm31752 in CEU-CHB-JPT lymphoblastoid cell lines                       |   5.941399| www.ncbi.nlm.nih.gov/pubmed/17873874 |
| HLA-DPB1 |     0| rs9277538 | Methylation levels at chr6:33151502-33151552 \[hg18 coord probe cg19921353\] in Cerebellum |  10.759451| www.ncbi.nlm.nih.gov/pubmed/20485568 |
| HLA-DPB1 |     0| rs9277538 | Differential splicing of HLA-DPA1 \[probeset 2950333\] in lymphoblastoid cell lines        |   8.096910| www.ncbi.nlm.nih.gov/pubmed/19052777 |
| HLA-DPB1 |     1| rs9296068 | Gene expression of HLA-DPB1 in blood                                                       |  39.207608| www.ncbi.nlm.nih.gov/pubmed/21829388 |
| HLA-DPB1 |     1| rs9296068 | Gene expression of HLA-DPB1 in peripheral blood monocytes                                  |  14.341989| www.ncbi.nlm.nih.gov/pubmed/20502693 |
| HLA-DPB1 |     1| rs9296068 | Gene expression of HLA-DOA in blood                                                        |   6.619789| www.ncbi.nlm.nih.gov/pubmed/21829388 |
| HLA-DPB1 |     1| rs9296068 | Gene expression of HLA-DPA1 in blood                                                       |   6.214670| www.ncbi.nlm.nih.gov/pubmed/21829388 |
| HLA-DPB1 |     1| rs9296068 | Gene expression of PSMB9 in blood                                                          |   3.585027| www.ncbi.nlm.nih.gov/pubmed/21829388 |
| HLA-DPB1 |     1| rs9296068 | Gene expression of HLA-DOB in blood                                                        |   3.000000| www.ncbi.nlm.nih.gov/pubmed/21829388 |
| HLA-DPB1 |     2| rs688209  | Gene expression of ITPR3 in peripheral blood monocytes                                     |   8.454693| www.ncbi.nlm.nih.gov/pubmed/20502693 |
| HLA-DPB1 |     2| rs688209  | Gene expression of DEFB125 in peripheral blood monocytes                                   |   5.892790| www.ncbi.nlm.nih.gov/pubmed/20502693 |
| HLA-DPB1 |     2| rs688209  | Gene expression of TAC3 in peripheral blood monocytes                                      |   5.754487| www.ncbi.nlm.nih.gov/pubmed/20502693 |
| HLA-DPB1 |     2| rs688209  | Gene expression of OR3A4 in peripheral blood monocytes                                     |   5.669586| www.ncbi.nlm.nih.gov/pubmed/20502693 |

### RTC

Here we can see that, when our eQTLs were not previously described in the databases, they tag some eQTL in the database.

Most of them are tagging an eQTL from regulomeDB. I believe that's because, for all sources not regulomeDB, I selected only the top variant for each tissue and gene. For regulomeDB that's no possible beucase the p-value is not directly available, so there are more variants per gene, which increases the chance that some of these variants will have high RTC with one of our eQTLs.

| gene     | rsid       |  rank| qtl\_previous |  d\_prime|   rtc| study        | tissue                                                                                                                                    |
|:---------|:-----------|-----:|:--------------|---------:|-----:|:-------------|:------------------------------------------------------------------------------------------------------------------------------------------|
| HLA-A    | rs16896724 |     0| rs2735097     |      0.99|  0.99| regulomeDB   | Lymphoblastoid                                                                                                                            |
| HLA-A    | rs2517827  |     1| rs2523759     |      1.00|  0.99| regulomeDB   | Lymphoblastoid                                                                                                                            |
| HLA-A    | rs2517827  |     1| rs2734980     |      1.00|  0.99| regulomeDB   | Lymphoblastoid                                                                                                                            |
| HLA-A    | rs2517827  |     1| rs2734984     |      1.00|  0.99| regulomeDB   | Lymphoblastoid                                                                                                                            |
| HLA-A    | rs2517827  |     1| rs5013088     |      1.00|  0.99| regulomeDB   | Lymphoblastoid,Monocytes                                                                                                                  |
| HLA-A    | rs7383537  |     2| rs5013310     |      0.84|  0.99| regulomeDB   | Lymphoblastoid                                                                                                                            |
| HLA-B    | rs1265094  |     0| rs1265181     |      1.00|  0.98| regulomeDB   | Monocytes                                                                                                                                 |
| HLA-B    | rs9264803  |     1| rs1050118     |      0.80|  1.00| Nedelec2016  | NonInfected\_Macrophages                                                                                                                  |
| HLA-C    | rs41561715 |     0| rs2073724     |      0.96|  0.98| regulomeDB   | Monocytes                                                                                                                                 |
| HLA-C    | rs2394941  |     1| rs2394894     |      1.00|  1.00| regulomeDB   | Monocytes                                                                                                                                 |
| HLA-C    | rs2249741  |     2| rs2249741     |      1.00|  1.00| regulomeDB   | Lymphoblastoid                                                                                                                            |
| HLA-DPB1 | rs9277538  |     0| rs9277538     |      1.00|  1.00| Westra2013   | Whole\_Blood                                                                                                                              |
| HLA-DPB1 | rs9296068  |     1| rs34885310    |      0.96|  0.98| Delaneau2017 | LCL                                                                                                                                       |
| HLA-DQA1 | rs9272316  |     0| rs9271488     |      0.96|  0.99| regulomeDB   | Lymphoblastoid                                                                                                                            |
| HLA-DQA1 | rs9271375  |     1| rs9268528     |      0.67|  0.99| regulomeDB   | Lymphoblastoid                                                                                                                            |
| HLA-DQB1 | rs1770     |     0| rs2647025     |      1.00|  1.00| regulomeDB   | Lymphoblastoid                                                                                                                            |
| HLA-DQB1 | rs9272209  |     1| rs2187668     |      0.98|  1.00| regulomeDB   | Monocytes                                                                                                                                 |
| HLA-DQB1 | rs9269386  |     2| rs3129872     |      0.56|  0.99| regulomeDB   | Monocytes                                                                                                                                 |
| HLA-DRB1 | rs3104412  |     0| rs1063355     |      0.99|  1.00| regulomeDB   | Lymphoblastoid                                                                                                                            |
| HLA-DRB1 | rs9270700  |     1| rs9270980     |      0.99|  1.00| GTEx2015\_v6 | Adipose\_Subcutaneous,Brain\_Anterior\_cingulate\_cortex\_BA24,Brain\_Cerebellum,Brain\_Frontal\_Cortex\_BA9,Breast\_Mammary\_Tissue,Lung |
| HLA-DRB1 | rs9270700  |     1| rs9271055     |      0.99|  1.00| regulomeDB   | Lymphoblastoid                                                                                                                            |

Association with GWAS traits
----------------------------

| gene     |  rank| variant    | gwas\_variant |   rtc| trait                                                                                                                                                  | studies                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
|:---------|-----:|:-----------|:--------------|-----:|:-------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| HLA-A    |     0| rs16896724 | rs2571391     |  0.99| serum IgE measurement                                                                                                                                  | www.ncbi.nlm.nih.gov/pubmed/22075330                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-A    |     1| rs2517827  | rs3095267     |  0.95| migraine disorder                                                                                                                                      | www.ncbi.nlm.nih.gov/pubmed/23793025                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-B    |     0| rs1265094  | rs3130573     |  0.92| systemic scleroderma                                                                                                                                   | www.ncbi.nlm.nih.gov/pubmed/21750679                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-B    |     1| rs9264803  | rs9264942     |  1.00| Crohn's disease/HIV-1 infection/inflammatory bowel disease                                                                                             | www.ncbi.nlm.nih.gov/pubmed/21051598 www.ncbi.nlm.nih.gov/pubmed/20041166 www.ncbi.nlm.nih.gov/pubmed/23128233 www.ncbi.nlm.nih.gov/pubmed/26192919 www.ncbi.nlm.nih.gov/pubmed/21051598 www.ncbi.nlm.nih.gov/pubmed/20041166 www.ncbi.nlm.nih.gov/pubmed/23128233 www.ncbi.nlm.nih.gov/pubmed/26192919 www.ncbi.nlm.nih.gov/pubmed/21051598 www.ncbi.nlm.nih.gov/pubmed/20041166 www.ncbi.nlm.nih.gov/pubmed/23128233 www.ncbi.nlm.nih.gov/pubmed/26192919                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| HLA-C    |     0| rs41561715 | rs9263963     |  0.97| serum IgG glycosylation measurement                                                                                                                    | www.ncbi.nlm.nih.gov/pubmed/23382691                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-C    |     1| rs2394941  | rs1265159     |  0.99| membranous glomerulonephritis                                                                                                                          | www.ncbi.nlm.nih.gov/pubmed/21323541                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-C    |     2| rs2249741  | rs9264942     |  1.00| Crohn's disease/HIV-1 infection/inflammatory bowel disease                                                                                             | www.ncbi.nlm.nih.gov/pubmed/21051598 www.ncbi.nlm.nih.gov/pubmed/20041166 www.ncbi.nlm.nih.gov/pubmed/23128233 www.ncbi.nlm.nih.gov/pubmed/26192919 www.ncbi.nlm.nih.gov/pubmed/21051598 www.ncbi.nlm.nih.gov/pubmed/20041166 www.ncbi.nlm.nih.gov/pubmed/23128233 www.ncbi.nlm.nih.gov/pubmed/26192919 www.ncbi.nlm.nih.gov/pubmed/21051598 www.ncbi.nlm.nih.gov/pubmed/20041166 www.ncbi.nlm.nih.gov/pubmed/23128233 www.ncbi.nlm.nih.gov/pubmed/26192919                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| HLA-DPB1 |     2| rs688209   | rs4711350     |  0.95| schizophrenia                                                                                                                                          | www.ncbi.nlm.nih.gov/pubmed/26198764                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DQA1 |     0| rs9272316  | rs2395185     |  0.98| antinuclear antibody measurement/Hodgkins lymphoma/lung carcinoma/ulcerative colitis                                                                   | www.ncbi.nlm.nih.gov/pubmed/19122664 www.ncbi.nlm.nih.gov/pubmed/22286212 www.ncbi.nlm.nih.gov/pubmed/19915573 www.ncbi.nlm.nih.gov/pubmed/23143601 www.ncbi.nlm.nih.gov/pubmed/25186300 www.ncbi.nlm.nih.gov/pubmed/20228799 www.ncbi.nlm.nih.gov/pubmed/19122664 www.ncbi.nlm.nih.gov/pubmed/22286212 www.ncbi.nlm.nih.gov/pubmed/19915573 www.ncbi.nlm.nih.gov/pubmed/23143601 www.ncbi.nlm.nih.gov/pubmed/25186300 www.ncbi.nlm.nih.gov/pubmed/20228799 www.ncbi.nlm.nih.gov/pubmed/19122664 www.ncbi.nlm.nih.gov/pubmed/22286212 www.ncbi.nlm.nih.gov/pubmed/19915573 www.ncbi.nlm.nih.gov/pubmed/23143601 www.ncbi.nlm.nih.gov/pubmed/25186300 www.ncbi.nlm.nih.gov/pubmed/20228799 www.ncbi.nlm.nih.gov/pubmed/19122664 www.ncbi.nlm.nih.gov/pubmed/22286212 www.ncbi.nlm.nih.gov/pubmed/19915573 www.ncbi.nlm.nih.gov/pubmed/23143601 www.ncbi.nlm.nih.gov/pubmed/25186300 www.ncbi.nlm.nih.gov/pubmed/20228799                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| HLA-DQA1 |     0| rs9272316  | rs9268853     |  0.98| lymphoma/ulcerative colitis                                                                                                                            | www.ncbi.nlm.nih.gov/pubmed/21297633 www.ncbi.nlm.nih.gov/pubmed/23349640 www.ncbi.nlm.nih.gov/pubmed/23511034 www.ncbi.nlm.nih.gov/pubmed/21297633 www.ncbi.nlm.nih.gov/pubmed/23349640 www.ncbi.nlm.nih.gov/pubmed/23511034                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 |
| HLA-DQA1 |     0| rs9272316  | rs9268905     |  0.98| Cystic fibrosis                                                                                                                                        | www.ncbi.nlm.nih.gov/pubmed/21602797                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DQA1 |     0| rs9272316  | rs9268923     |  0.98| ulcerative colitis/NA                                                                                                                                  | www.ncbi.nlm.nih.gov/pubmed/20228798 www.ncbi.nlm.nih.gov/pubmed/26819262 www.ncbi.nlm.nih.gov/pubmed/20228798 www.ncbi.nlm.nih.gov/pubmed/26819262                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| HLA-DQA1 |     1| rs9271375  | rs4321864     |  0.99| Staphylococcus aureus infection                                                                                                                        | www.ncbi.nlm.nih.gov/pubmed/26450422                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DQB1 |     0| rs1770     | rs6927022     |  1.00| ulcerative colitis                                                                                                                                     | www.ncbi.nlm.nih.gov/pubmed/23128233                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DQB1 |     1| rs9272209  | rs2187668     |  1.00| autoimmune hepatits type 1/celiac disease/cutaneous lupus erythematosus/membranous glomerulonephritis/protein measurement/systemic lupus erythematosus | www.ncbi.nlm.nih.gov/pubmed/17558408 www.ncbi.nlm.nih.gov/pubmed/18204098 www.ncbi.nlm.nih.gov/pubmed/20694011 www.ncbi.nlm.nih.gov/pubmed/21408207 www.ncbi.nlm.nih.gov/pubmed/20190752 www.ncbi.nlm.nih.gov/pubmed/24768677 www.ncbi.nlm.nih.gov/pubmed/21323541 www.ncbi.nlm.nih.gov/pubmed/25827949 www.ncbi.nlm.nih.gov/pubmed/26316170 www.ncbi.nlm.nih.gov/pubmed/17558408 www.ncbi.nlm.nih.gov/pubmed/18204098 www.ncbi.nlm.nih.gov/pubmed/20694011 www.ncbi.nlm.nih.gov/pubmed/21408207 www.ncbi.nlm.nih.gov/pubmed/20190752 www.ncbi.nlm.nih.gov/pubmed/24768677 www.ncbi.nlm.nih.gov/pubmed/21323541 www.ncbi.nlm.nih.gov/pubmed/25827949 www.ncbi.nlm.nih.gov/pubmed/26316170 www.ncbi.nlm.nih.gov/pubmed/17558408 www.ncbi.nlm.nih.gov/pubmed/18204098 www.ncbi.nlm.nih.gov/pubmed/20694011 www.ncbi.nlm.nih.gov/pubmed/21408207 www.ncbi.nlm.nih.gov/pubmed/20190752 www.ncbi.nlm.nih.gov/pubmed/24768677 www.ncbi.nlm.nih.gov/pubmed/21323541 www.ncbi.nlm.nih.gov/pubmed/25827949 www.ncbi.nlm.nih.gov/pubmed/26316170 www.ncbi.nlm.nih.gov/pubmed/17558408 www.ncbi.nlm.nih.gov/pubmed/18204098 www.ncbi.nlm.nih.gov/pubmed/20694011 www.ncbi.nlm.nih.gov/pubmed/21408207 www.ncbi.nlm.nih.gov/pubmed/20190752 www.ncbi.nlm.nih.gov/pubmed/24768677 www.ncbi.nlm.nih.gov/pubmed/21323541 www.ncbi.nlm.nih.gov/pubmed/25827949 www.ncbi.nlm.nih.gov/pubmed/26316170 www.ncbi.nlm.nih.gov/pubmed/17558408 www.ncbi.nlm.nih.gov/pubmed/18204098 www.ncbi.nlm.nih.gov/pubmed/20694011 www.ncbi.nlm.nih.gov/pubmed/21408207 www.ncbi.nlm.nih.gov/pubmed/20190752 www.ncbi.nlm.nih.gov/pubmed/24768677 www.ncbi.nlm.nih.gov/pubmed/21323541 www.ncbi.nlm.nih.gov/pubmed/25827949 www.ncbi.nlm.nih.gov/pubmed/26316170 www.ncbi.nlm.nih.gov/pubmed/17558408 www.ncbi.nlm.nih.gov/pubmed/18204098 www.ncbi.nlm.nih.gov/pubmed/20694011 www.ncbi.nlm.nih.gov/pubmed/21408207 www.ncbi.nlm.nih.gov/pubmed/20190752 www.ncbi.nlm.nih.gov/pubmed/24768677 www.ncbi.nlm.nih.gov/pubmed/21323541 www.ncbi.nlm.nih.gov/pubmed/25827949 www.ncbi.nlm.nih.gov/pubmed/26316170 |
| HLA-DQB1 |     2| rs9269386  | rs3021304     |  1.00| Vogt-Koyanagi-Harada disease                                                                                                                           | www.ncbi.nlm.nih.gov/pubmed/25108386                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DRB1 |     0| rs3104412  | rs1063355     |  1.00| ulcerative colitis                                                                                                                                     | www.ncbi.nlm.nih.gov/pubmed/24837172                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DRB1 |     1| rs9270700  | rs9270984     |  1.00| systemic lupus erythematosus                                                                                                                           | www.ncbi.nlm.nih.gov/pubmed/23273568                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| HLA-DRB1 |     1| rs9270700  | rs9270986     |  1.00| Myasthenia gravis                                                                                                                                      | www.ncbi.nlm.nih.gov/pubmed/25643325                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |

Trans-eQTLs
===========

-   Approximate pass as described on QTLtools website

| gene\_name | variant\_id |  variant\_chr|  variant\_pos|  log10\_nom\_pval|
|:-----------|:------------|-------------:|-------------:|-----------------:|
| HLA-DQA1   | rs3997880   |             2|     178431654|             31.06|
| HLA-DQA1   | rs145339215 |             2|     178431521|             30.62|
| HLA-DQA1   | rs3997878   |             2|     178431899|             27.23|
| HLA-DQA1   | rs141354030 |             2|     178451007|             20.85|

HLA lineages
------------

![](./qtls/plots/lineage_and_effects.png)

### F-test: is there a difference between lineages?

#### traditional ANOVA

| locus |   df|        F| p.value  |
|:------|----:|--------:|:---------|
| A     |   14|   10.523| 7.49e-22 |
| B     |   25|    4.859| 2.80e-13 |
| C     |   13|   27.398| 2.62e-54 |
| DPB1  |   21|    7.950| 8.72e-22 |
| DQA1  |    5|  107.821| 1.16e-84 |
| DQB1  |    4|   37.589| 1.52e-28 |
| DRB1  |   12|   27.992| 4.05e-52 |

#### Welch ANOVA

| locus |  num.df|  denom.df|        F| p.value  |
|:------|-------:|---------:|--------:|:---------|
| A     |      14|    55.246|   10.936| 2.53e-11 |
| B     |      23|    63.121|    6.118| 5.29e-09 |
| C     |      12|    98.242|   19.612| 4.86e-21 |
| DPB1  |      18|    35.687|    5.332| 1.06e-05 |
| DQA1  |       4|   213.228|  105.181| 2.54e-49 |
| DQB1  |       4|   257.192|   42.934| 1.43e-27 |
| DRB1  |      12|   100.543|   43.054| 3.61e-34 |
