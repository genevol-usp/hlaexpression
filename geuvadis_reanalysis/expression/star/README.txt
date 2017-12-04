# How to run pipeline

$ qsub map_and_quantify_round1.pbs

$ Rscript process_quants.R 1

$ ./write_custom_index.sh

$ qsub map_and_quantify_round2.pbs

$ Rscript process_quants.R 2

$ ./compile_quantifications.sh ./quantifications_2

$ Rscript write_quantifications_bed.R ./quantifications_2
