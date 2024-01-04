#! /bin/bash

export PROJECT_DIR=/home/cbaker/Hydro/Summa-openWQ-revised-install
export DATA_DIR=/home/cbaker/Hydro/SUMMAapptainer/case_studies/mizuroute_Great_Slave_Lake

docker run -d -it --name SUMMA-openWQ \
    --mount type=bind,source=${PROJECT_DIR},target=/code/Summa-OpenWQ \
    --mount type=bind,source=${DATA_DIR},target=/code/data \
    summa-openwq:latest