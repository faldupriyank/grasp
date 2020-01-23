#!/bin/bash

make clean
make cleansrc
make DEBUG= -j

ROUNDS=3
ROOT=31
export OMP_PROC_BIND=true

TIMESTAMP=$(date +%y-%m-%d-%H-%M-%S)
LOG=${DBG_ROOT}/logs
EXP_NAME="Run-All-APPs-On-Web-Google"
EXP_NAME=${EXP_NAME}-${TIMESTAMP}
O_DIR=${LOG}/${EXP_NAME}
mkdir -p ${O_DIR}

#for DATASET in web-Google soc-LiveJournal1 dbpedia-link twitter_rv twitter_mpi friendster pld-arc sd-arc road-raod-usa; do
for DATASET in web-Google; do
    echo "DATASET: ${DATASET}"
    for app in PageRankOrig PageRankOpt PageRankDeltaOrig PageRankDeltaOpt BellmanFordOrig-iters BellmanFordOpt-iters; do
        echo "app: ${app}"
        DEG="out"
        DEGREE_USED_FOR_REORDERING=0 #Out-Degree
        if [ "${app}" == "PageRankDelta" ] || [ "${app}" == "BellmanFordOrig" ] || [ "${app}" == "BellmanFordOrig-iters" ] || [ "${app}" == "BellmanFordOpt" ] || [ "${app}" == "BellmanFordOpt-iters" ]; then
            DEGREE_USED_FOR_REORDERING=1 #In-Degree
            DEG="in"
        fi
        echo "DEGREE_USED_FOR_REORDERING: ${DEGREE_USED_FOR_REORDERING}"
        for REORDERING_ALGO in 5; do
            echo "REORDERING_ALGO: ${REORDERING_ALGO}"
                make DEGREE_USED_FOR_REORDERING=${DEGREE_USED_FOR_REORDERING} REORDERING_ALGO=${REORDERING_ALGO} ROUNDS=${ROUNDS} ROOT=${ROOT} DATASET=${DATASET} MAP_FILE=${DBG_ROOT}/datasets/${DATASET}.${DEG}.map run-${app}  > ${O_DIR}/${app}.${DATASET}.${DEGREE_USED_FOR_REORDERING}.${REORDERING_ALGO}.log 2>&1
       done
    done
done
