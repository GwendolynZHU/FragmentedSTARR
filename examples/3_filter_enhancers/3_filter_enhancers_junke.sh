FULL=distal.full
TATA_DOWN=distal.tata.down
TATA_UP=distal.tata.up
PAUSE_DOWN=distal.pause.down
PAUSE_UP=distal.pause.up

FULL_DATA=data/Junke_architecture/K562_full.bed
TATA_DOWN_DATA=data/Junke_architecture/K562_delta_TATA_downstream.bed
TATA_UP_DATA=data/Junke_architecture/K562_delta_TATA_upstream.bed
PAUSE_DOWN_DATA=data/Junke_architecture/K562_delta_pause_downstream.bed
PAUSE_UP_DATA=data/Junke_architecture/K562_delta_pause_upstream.bed

GROCAP_PL=/fs/cbsuhy01/storage/jz855/Reference/K562/GRO_cap/K562_GROcap_hg38_aligned_pl.bw
GROCAP_MN=/fs/cbsuhy01/storage/jz855/Reference/K562/GRO_cap/K562_GROcap_hg38_aligned_mn.bw

python src/3_filter_analysis.py \
    --outdir output/Junke_architecture \
    --design ${FULL} ${TATA_DOWN} ${TATA_UP} ${PAUSE_DOWN} ${PAUSE_UP} \
    --starr deep-ATAC-STARR \
    --resolution 5 \
    --either \
    --padj 0.05 \
    --annotate_extra \
    --original_bed ${FULL_DATA} ${TATA_DOWN_DATA} ${TATA_UP_DATA} ${PAUSE_DOWN_DATA} ${PAUSE_UP_DATA} \
    --annotate_bigwig \
    --bigwig ${GROCAP_PL} ${GROCAP_MN} \
    --core_promoter
