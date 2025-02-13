FULL=distal.full
TATA_DOWN=distal.tata.down
TATA_UP=distal.tata.up
PAUSE_DOWN=distal.pause.down
PAUSE_UP=distal.pause.up

./src/4_visualization.r \
    --output output/Junke_architecture \
    --baseline ${FULL} \
    --comparison ${TATA_DOWN} ${TATA_UP} ${PAUSE_DOWN} ${PAUSE_UP} \
    --double_sided \
    --starr deep-ATAC-STARR \
    --resolution 5