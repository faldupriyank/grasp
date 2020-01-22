#!/bin/bash

REORDER=dbg
CSIZE=1
for d in web-Google; do
    echo $d
    for app in BellmanFordOpt PageRankOpt Radii BC PageRankDeltaOpt; do
        echo $app
        ext=cvgr
        if [[ ${app} == "BellmanFordOpt" ]]; then
            ext=cintgr
        fi
        echo $d
        for p in lru pin grasp belady; do
            echo $p
            make POLICY=${p} TRACE=${app}.${d}.${ext}.${REORDER}.lru.llc.trace CSIZE=${CSIZE} > ${DBG_ROOT}/logs/${app}.${d}.${ext}.${REORDER}.lru.llc.trace.${p}.${CSIZE}
            echo $done
            echo
        done
    done
done
