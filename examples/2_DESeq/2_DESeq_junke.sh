FULL=distal.full
TATA_DOWN=distal.tata.down
TATA_UP=distal.tata.up
PAUSE_DOWN=distal.pause.down
PAUSE_UP=distal.pause.up

./src/2_DESeq_STARR.r \
    --output output/Junke_architecture \
    --design ${FULL} ${TATA_DOWN} ${TATA_UP} ${PAUSE_DOWN} ${PAUSE_UP} \
    --starr deep-ATAC-STARR \
    --resolution 5 \
    --cutoff 5