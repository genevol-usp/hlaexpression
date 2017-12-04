# How to run pipeline

$ qsub map_and_quantify_round1.pbs

$ Rscript process_quants.R 1

$ ./write_custom_index.sh

$ qsub map_and_quantify_round2.pbs

$ Rscript process_quants.R 2

$ Rscript make_quantifications_table.R quantifications_2
