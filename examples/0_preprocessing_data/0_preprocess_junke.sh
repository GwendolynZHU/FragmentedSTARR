GROCAP_DIVERGENT=data/GRO_cap_PINTS_qpeak_call_117/K562_GROcap_hg38_1.1.7_qpeak_calls_1_divergent_peaks_element_60bp.bed
GENCODE_TSS_500=/fs/cbsuhy01/storage/jz855/Reference/hg38/proximal_v45/promoter_1kbp_protein_coding_TSS_centered_gencode_v45.bed

FULL=data/Junke_architecture/K562_full.bed
TATA_DOWN=data/Junke_architecture/K562_delta_TATA_downstream.bed
TATA_UP=data/Junke_architecture/K562_delta_TATA_upstream.bed

PAUSE_DOWN=data/Junke_architecture/K562_delta_pause_downstream.bed
PAUSE_UP=data/Junke_architecture/K562_delta_pause_upstream.bed

# Filter the divergent peaks
# - Remove peaks that are proximal (90% within 500 bp) to TSS
bedtools intersect -v -a $GROCAP_DIVERGENT -b $GENCODE_TSS_500 -f 0.9 > $GROCAP_DIVERGENT.tmp

# Generate the divergent peaks for partial elements
awk 'BEGIN{OFS="\t"}{$2 = $7 + 32; print}' $GROCAP_DIVERGENT.tmp > $TATA_UP.tmp
awk 'BEGIN{OFS="\t"}{$3 = $8 - 32; print}' $GROCAP_DIVERGENT.tmp > $TATA_DOWN.tmp

awk 'BEGIN{OFS="\t"}{$2 = $7 - 17; print}' $GROCAP_DIVERGENT.tmp > $PAUSE_UP.tmp
awk 'BEGIN{OFS="\t"}{$3 = $8 + 17; print}' $GROCAP_DIVERGENT.tmp > $PAUSE_DOWN.tmp

# Filter:
# - Remove peaks with negative coordinates
# - Remove peaks shorter than 40 bp
# - Remove peaks in chrM or random chromosomes
awk '!($1 == "chrM" || ($3 - $2 < 40) || $2 < 0 || index($1, "random") > 0)' $TATA_UP.tmp > $TATA_UP
awk '!($1 == "chrM" || ($3 - $2 < 40) || $2 < 0 || index($1, "random") > 0)' $TATA_DOWN.tmp > $TATA_DOWN

awk '!($1 == "chrM" || ($3 - $2 < 40) || $2 < 0 || index($1, "random") > 0)' $PAUSE_UP.tmp > $PAUSE_UP
awk '!($1 == "chrM" || ($3 - $2 < 40) || $2 < 0 || index($1, "random") > 0)' $PAUSE_DOWN.tmp > $PAUSE_DOWN

awk '!($1 == "chrM" || ($3 - $2 < 40) || $2 < 0 || index($1, "random") > 0)' $GROCAP_DIVERGENT.tmp > $FULL

# Remove temporary files
rm $TATA_UP.tmp $TATA_DOWN.tmp $PAUSE_UP.tmp $PAUSE_DOWN.tmp $GROCAP_DIVERGENT.tmp