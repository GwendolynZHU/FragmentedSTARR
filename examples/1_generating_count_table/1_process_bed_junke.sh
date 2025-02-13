FULL=data/Junke_architecture/K562_full.bed

TATA_DOWN=data/Junke_architecture/K562_delta_TATA_downstream.bed
TATA_UP=data/Junke_architecture/K562_delta_TATA_upstream.bed

PAUSE_DOWN=data/Junke_architecture/K562_delta_pause_downstream.bed
PAUSE_UP=data/Junke_architecture/K562_delta_pause_upstream.bed

python src/1_process_bed.py \
    --input_path ${FULL} \
    --file_source deep-ATAC-STARR \
    --design distal.full \
    --outdir output/Junke_architecture \
    --resolution 5

python src/1_process_bed.py \
    --input_path ${TATA_DOWN} \
    --file_source deep-ATAC-STARR \
    --design distal.tata.down \
    --outdir output/Junke_architecture \
    --resolution 5

python src/1_process_bed.py \
    --input_path ${TATA_UP} \
    --file_source deep-ATAC-STARR \
    --design distal.tata.up \
    --outdir output/Junke_architecture \
    --resolution 5

python src/1_process_bed.py \
    --input_path ${PAUSE_DOWN} \
    --file_source deep-ATAC-STARR \
    --design distal.pause.down \
    --outdir output/Junke_architecture \
    --resolution 5

python src/1_process_bed.py \
    --input_path ${PAUSE_UP} \
    --file_source deep-ATAC-STARR \
    --design distal.pause.up \
    --outdir output/Junke_architecture \
    --resolution 5