#! /bin/bash

export PROJECT_DIR=/home/cbaker/Hydro/Summa-openWQ

docker run -d -it --name SUMMA-openWQ --mount type=bind,source=${PROJECT_DIR},target=/code/Summa-OpenWQ \
    summa-openwq:latest