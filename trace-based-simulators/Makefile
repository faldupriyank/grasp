CXX=/usr/bin/g++
CFLAGS=
LDFLAGS=
PWD=$(shell pwd)

#size of cache in MB
CSIZE=1

# Supported values: grasp, lru, belady, pin
POLICY=grasp

# Path to a directory where the trace files are located
GDIR=${DBG_ROOT}/datasets

# Name of a trace file
TRACE=PageRankOpt.web-Google.cvgr.dbg.lru.llc.trace

#CFLAGS=-std=c++

run: ${POLICY}.bin
	${PWD}/${POLICY}.bin ${GDIR}/${TRACE} ${CSIZE}

${POLICY}.bin: ${POLICY}.cpp common.h cache.h
	${CXX} ${CFLAGS} ${LDFLAGS} -o ${POLICY}.bin ${POLICY}.cpp 

clean:
	rm -rf *.bin

help:
	@echo TRACE=${TRACE};\
		echo POLICY=grasp;\
		echo CSIZE=1;\
		echo GDIR=${MEMPHIS_DATASETS};\




