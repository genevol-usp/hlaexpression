JOBS=8

parallel --gnu -j $JOBS ./run_hlatx.sh {} :::: samples_phase3_ena.txt

rm -r results/gencounts results/hla results/index results/map \
  results/xgencounts results/xhla results/xmap
