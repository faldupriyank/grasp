#CONFIGS set by Priyank
LONG=1 #vertex is 8 bytes long
EDGELONG= #edge is 4 bytes long
OPENMP=1 #use OPENMP for parallelism
VGR=1 #use VGR formatted graph datasets instead of defult adjacency list

ifdef LONG
INTT = -DLONG
endif

ifdef EDGELONG
INTE = -DEDGELONG
endif

ifdef PD
PD = -DPD
endif

ifdef BYTE
CODE = -DBYTE
else ifdef NIBBLE
CODE = -DNIBBLE
else
CODE = -DBYTERLE
endif

ifdef LOWMEM
MEM = -DLOWMEM
endif

ifdef ALIGNED
ALIGN = -DALIGNED
endif

ifdef VGR
VGR = -DVGR
endif

PCC = g++
PCFLAGS = -std=c++14 -fopenmp -march=native -O3 -DOPENMP $(INTT) $(INTE) $(CODE) $(PD) $(MEM) $(ALIGN) $(VGR) #-D_OUTPUT_ -D_SNIPER_

COMMON_FILES= ligra.h graph.h compressedVertex.h vertex.h utils.h IO.h parallel.h gettime.h timer.h index_map.h maybe.h sequence.h edgeMap_utils.h binary_search.h quickSort.h blockRadixSort.h transpose.h parseCommandLine.h byte.h byteRLE.h nibble.h byte-pd.h byteRLE-pd.h nibble-pd.h vertexSubset.h encoder.C decoder.C pvector.h dbg.h sim_api.h

ALL= BC BellmanFordOrig PageRankOrig PageRankOpt PageRankDeltaOpt PageRankDeltaOrig Radii BC-iters BellmanFordOrig-iters BellmanFordOpt BellmanFordOpt-iters

all: $(ALL)

% : %.C $(COMMON_FILES)
	$(PCC) $(PCFLAGS) -o $@ $<

$(COMMON_FILES):
	ln -s ../ligra/$@ .

.PHONY : clean

clean :
	rm -f *.o $(ALL)

cleansrc :
	rm -f *.o $(ALL)
	rm $(COMMON_FILES)


COMMON=
#REORDERING_ALGO
#ORIGINAL=0
#Random=1
#Sort=2
#HubSort=3
#HubCluster=4
#DBG=5
#HubSortDBG=6
#HubClusterDBG=7
#MAP=10
ifndef REORDERING_ALGO
REORDERING_ALGO=5
endif
ifndef DEGREE_USED_FOR_REORDERING
DEGREE_USED_FOR_REORDERING=0 #OUTDEGREE, set to 1 for INDEGREE
endif
ifndef MAXITERS
MAXITERS=100
endif
ifndef ROUNDS
ROUNDS=3
endif
ifndef RAND_GRAN
RAND_GRAN=1
endif
ifndef NUM_ROOTS
NUM_ROOTS=8
endif
ifndef ROOT
ROOT=31
endif
ifndef LIGRA_ROOT
	LIGRA_ROOT=.
endif
ifndef GPATH
	GPATH=${DBG_ROOT}/datasets
endif
ifndef DATASET
	DATASET=web-Google
endif

ifdef THREADS
COMMON+=-threads ${THREADS}
endif
ifdef MAP_FILE
COMMON+=-map_file ${MAP_FILE}
endif

ifeq ($(THREADS), 1)
NUMA=numactl -N 0 -i 0
TASKSET=taskset -c 0
else
NUMA=numactl -N all -i all
TASKSET=
endif

COMMON+=-degree_used_for_reordering ${DEGREE_USED_FOR_REORDERING}
COMMON+=-rounds ${ROUNDS}
COMMON+=-rand_gran ${RAND_GRAN}
COMMON+=-num_roots ${NUM_ROOTS}

CVGR=.cvgr
CINTGR=.cintgr
CSVGR=.csvgr

run-PageRankOrig: # PULL (OUT-DEGREE)
	${NUMA} ${TASKSET} ${DBG_ROOT}/apps/PageRankOrig -reordering_algo ${REORDERING_ALGO} -maxiters ${MAXITERS} -is_pagerank 1 ${COMMON} ${GPATH}/${DATASET}${CVGR}

run-PageRankOpt: # PULL (OUT-DEGREE)
	${NUMA} ${TASKSET} ${DBG_ROOT}/apps/PageRankOpt -reordering_algo ${REORDERING_ALGO} -maxiters ${MAXITERS} -is_pagerank 1 ${COMMON} ${GPATH}/${DATASET}${CVGR}

run-PageRankDeltaOrig: # PUSH (IN-DEGREE)
	${NUMA} ${TASKSET} ${DBG_ROOT}/apps/PageRankDeltaOrig -reordering_algo ${REORDERING_ALGO} -maxiters ${MAXITERS} -is_dense_write 1 ${COMMON} ${GPATH}/${DATASET}${CVGR}

run-PageRankDeltaOpt: # PUSH (IN-DEGREE)
	${NUMA} ${TASKSET} ${DBG_ROOT}/apps/PageRankDeltaOpt -reordering_algo ${REORDERING_ALGO} -maxiters ${MAXITERS} ${COMMON} ${GPATH}/${DATASET}${CVGR}

run-Radii: # PULL-PUSH (OUT-DEGREE)
	${NUMA} ${TASKSET} ${DBG_ROOT}/apps/Radii -reordering_algo ${REORDERING_ALGO} ${COMMON} ${GPATH}/${DATASET}${CVGR}

run-BellmanFordOrig: #PUSH (IN_DEGREE) WEIGHTED
	${NUMA} ${TASKSET} ${DBG_ROOT}/apps/BellmanFordOrig -reordering_algo ${REORDERING_ALGO} -is_dense_write 1 -r ${ROOT} ${COMMON} ${GPATH}/${DATASET}${CINTGR}

run-BellmanFordOrig-iters: #PUSH (IN_DEGREE) WEIGHTED
	${NUMA} ${TASKSET} ${DBG_ROOT}/apps/BellmanFordOrig-iters -reordering_algo ${REORDERING_ALGO} -is_dense_write 1 ${COMMON} ${GPATH}/${DATASET}${CINTGR}

run-BellmanFordOpt: #PUSH (IN_DEGREE) WEIGHTED
	${NUMA} ${TASKSET} ${DBG_ROOT}/apps/BellmanFordOpt -reordering_algo ${REORDERING_ALGO} -is_dense_write 1 -r ${ROOT} ${COMMON} ${GPATH}/${DATASET}${CINTGR}

run-BellmanFordOpt-iters: #PUSH (IN_DEGREE) WEIGHTED
	${NUMA} ${TASKSET} ${DBG_ROOT}/apps/BellmanFordOpt-iters -reordering_algo ${REORDERING_ALGO} -is_dense_write 1 ${COMMON} ${GPATH}/${DATASET}${CINTGR}

run-BC: # PULL-PUSH (OUT_DEGREE)
	${NUMA} ${TASKSET} ${DBG_ROOT}/apps/BC -reordering_algo ${REORDERING_ALGO} -r ${ROOT} ${COMMON} ${GPATH}/${DATASET}${CVGR}

run-BC-iters: # PULL-PUSH (OUT_DEGREE)
	${NUMA} ${TASKSET} ${DBG_ROOT}/apps/BC-iters -reordering_algo ${REORDERING_ALGO} ${COMMON} ${GPATH}/${DATASET}${CVGR}

