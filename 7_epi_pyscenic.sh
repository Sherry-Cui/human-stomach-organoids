# D16 epi
export SCENIC_DATA="/data2/scenic"
export EPI_DATA="$SCENIC_DATA/epi"
export DB_DATA="$SCENIC_DATA/cistarget_database/human_hg38"

nohup pyscenic grn \
--num_workers 20 \
--output adj.sample.tsv \
--method grnboost2 \
$EPI_DATA/sample.loom \
$DB_DATA/hs_hgnc_tfs.txt & 

nohup pyscenic ctx \
adj.sample.tsv $DB_DATA/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
--annotations_fname $DB_DATA/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname $EPI_DATA/sample.loom \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 20 \
--mask_dropouts &

export OMP_NUM_THREADS=1

nohup pyscenic aucell \
$EPI_DATA/sample.loom \
$EPI_DATA/reg.csv \
--output out_SCENIC.loom \
--num_workers 1 &