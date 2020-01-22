/*
 * @author Priyank Faldu <Priyank.Faldu@ed.ac.uk> <http://faldupriyank.com>
 *
 * Copyright 2020 The University of Edinburgh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _CACHE_H_
#define _CACHE_H_

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <map>
#include <cmath>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <random>
#include "common.h"

using namespace std;

struct block_t {
    uint64_t addr;
    int64_t time;
    int rrip;
    bool pin;
};

int M_RRIP_INIT;
long long unsigned cache_size;
int assoc;
int sets;
int num_blocks;
uint64_t total_addresses;
uint64_t* trace;
block_t** blocks; 
uint64_t misses;
uint64_t hits;
int set_shift;
int block_shift;

const int max_regions = 10; // Increase it as needed
int num_regions = 0;
region_t regions[max_regions];

const int max_border_regions = 2; // Increase it as needed
int num_border_regions = 0;
region_t border_regions[max_border_regions];

// File must contain at least six regions
// propertyA-0
// propertyA-n
// propertyA-f
// propertyB-0
// propertyB-n
// propertyB-f
void add_all_regions() {

    memset(regions, 0, sizeof(region_t) * max_regions);
    memset(border_regions, 0, sizeof(region_t) * max_border_regions);
    num_regions = 0;

    if ( magic_map["propertyA-f"] > 0 ) {
        regions[num_regions].min = magic_map[string("propertyA-0")];
        regions[num_regions].max = magic_map[string("propertyA-n")];
        strcpy(regions[num_regions].str, "propertyA");
        num_regions++;
        add_border_regions(&border_regions[num_border_regions++], "propertyA", magic_map["propertyA-f"]);
    }

    if ( magic_map["propertyB-f"] > 0 ) {
        regions[num_regions].min = magic_map[string("propertyB-0")];
        regions[num_regions].max = magic_map[string("propertyB-n")];
        strcpy(regions[num_regions].str, "propertyB");
        num_regions++;
        add_border_regions(&border_regions[num_border_regions++], "propertyB", magic_map["propertyB-f"]);
    }

    assert(num_regions <= max_regions);
    printf("%35s %20s %20s\n", "Property Region", "Size(bytes)", "Size(frac)");

    assert(num_border_regions <= max_border_regions);
    add_border_boundry(border_regions, num_border_regions, cache_size);
}


void read_file(FILE* fp, const string& filename) {
// File format
// #regions (8 bytes) --> number of regions
// Followed by #regions records, each record occupying 25+8=33 bytes.
//      Each record contains name of the region(25 bytes) and value(8 bytes). 
// Followed by #address (8 bytes) --> number of L2 misses
// Followed by a trace of addresses, 8 bytes each.

// Region name convention:
//  Region-1 --> lower bound address for the region
//  Region-n --> upper bound address for the region
//  Region-f --> percentage of LLC capacity to be allocated to this region
//  Any other region can be ignored.
    assert(fp!=NULL);
    uint64_t res=0;
    uint64_t len=0;
    res = fread(&len, 8, 1, fp);
    char key[25];
    uint64_t val;
    printf("%35s %20s\n", "Region", "Value");
    for ( int i = 0 ; i < len ; i++ ) {
        res = fread(&key, 1, 25, fp);
        res = fread(&val, 8, 1, fp);
        magic_map[string(key)] = val;
        printf("%35s %20lu\n", key, val);
    }
    printf("\n");
    
    add_all_regions();

    fread(&total_addresses, 8, 1, fp);

    printf("%35s %20lu\n", "Total Addresses:", total_addresses);

    trace = (uint64_t*) malloc(total_addresses * 8);
    uint64_t bytes = fread(trace, 8, total_addresses, fp);
    printf("%35s %20lu\n", "Elements read:", bytes);
    if ( bytes != total_addresses ) {
        printf("Could not read all addresses.\n");
    }   
    assert(trace);
}

void init_cache(int argc, char* argv[]) {
    if ( argc != 3 ) {
        printf("Wrong number of arguments! Usage: ./exe filename cache_size(in MB)\n");
        exit(-1);
    }

    cache_size = stoi(argv[2]) * 1024 * 1024; // cache size in MB
    num_blocks = cache_size / 64;
    assoc = 16;
    sets = num_blocks / assoc;
    blocks = (block_t**)malloc(sizeof(block_t*)*sets);
    assert(blocks);

    for ( int i = 0 ; i < sets ; i++ ) {
        blocks[i] = (block_t*)malloc(sizeof(block_t)*assoc);
        assert(blocks[i]);
    }

    for ( uint64_t i = 0 ; i < sets ; i++ ) {
        for ( int j = 0 ; j < assoc ; j++ ) {
            blocks[i][j].addr = 0;
            blocks[i][j].time = -1;
            blocks[i][j].pin = false;
            blocks[i][j].rrip = M_RRIP_INIT;
        }
    }
    set_shift = log2(sets);
    block_shift = log2(64);

    printf("%35s %20s\n", "Cache Property", "Value");
    printf("%35s %20llu\n", "cache size (MB)", cache_size / (1024 * 1024));
    printf("%35s %20llu\n", "cache size (Bytes)", cache_size);
    printf("%35s %20d\n", "sets", sets);
    printf("%35s %20d\n", "num_blocks", num_blocks);
    printf("%35s %20d\n", "Assoc", assoc);
    printf("%35s %20d\n", "set_shift", set_shift);
    printf("%35s %20d\n", "block_shift", block_shift);
    printf("\n");

    FILE* fp = fopen(argv[1], "rb");
    if ( !fp ) {
        printf("Could not open file: %s\n", argv[1]);
    }
    read_file(fp, argv[1]);
}

inline int get_set(uint64_t trace, int block_shift, int sets) {
    return (trace >> block_shift) % sets;
}

bool find_block(uint64_t tr, int block_shift, int sets, block_t** _block, block_t** _set, int* set_num = NULL) {
    int set_index = get_set(tr, block_shift, sets);
    block_t* set = &blocks[set_index][0];
    *_set = set;
    if ( set_num ) {
        *set_num = set_index;
    }

    bool hit = false;
    for ( int t = 0 ; t < assoc ; t++ ) {
        if ( set[t].addr == tr ) {
            hit = true;
            block_t* block = &set[t];
            *_block = block;
            break;
        }
    }
    return hit;
}

inline void update_compute_stats(uint64_t tr, bool hit) {
    if (!hit) {
        misses++;
    }
}

void display() {
    printf("\n------------------------------------------------------------------\n");
    printf("%35s %10lu\n", "total-accesses:", total_addresses);
    printf("%35s %10lu\n", "total-misses:", misses);
    printf("------------------------------------------------------------------\n");
    printf("%35s %10.6f\n", "miss-rate:", misses/(1.0*total_addresses));
    printf("------------------------------------------------------------------\n");
}

#endif
