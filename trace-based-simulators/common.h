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

#ifndef _COMMON_H_
#define _COMMON_H_

#include <map>

using namespace std;

map<string, uint64_t> magic_map;
struct region_t {
    uint64_t min;
    uint64_t max;
    uint64_t border_high_reuse;
    uint64_t border_moderate_reuse;
    char str[25]; // For display purpose only
};

inline bool is_in_high_reuse_region(const uint64_t addr, region_t regions[], const int num_regions) {
    for ( int i = 0 ; i < num_regions ; i++ ) {
        if ( (addr >= regions[i].min) && (addr < regions[i].border_high_reuse) ) {
            return true;
        }
    }
    return false;
}

inline bool is_in_moderate_reuse_region(const uint64_t addr, region_t regions[], const int num_regions) {
    for ( int i = 0 ; i < num_regions ; i++ ) {
        if ( (addr >= regions[i].border_high_reuse) && (addr < regions[i].border_moderate_reuse) ) {
            return true;
        }
    }
    return false;
}

int add_border_regions(region_t* regions, const string& str, int f = -1) {
    if ( f > 0 ) {
        regions->min = magic_map[str+"-0"];
        regions->max = magic_map[str+"-n"];
        strcpy(regions->str, str.c_str());
        regions->border_high_reuse = f;
    }
    return 0;
}

int add_border_boundry(region_t regions[], int num_regions, int cache_size) {
    long long unsigned int total_capacity = cache_size;
    long long unsigned int total_used = 0;
    const float factor = 1;
    cache_size = cache_size * factor;
    for ( int i = 0 ; i < num_regions; i++ ) {
        int f = regions[i].border_high_reuse * cache_size / 100.0;
        printf("%35s %20d %20f\n", regions[i].str, f, f * 1.0 / cache_size);
        regions[i].border_high_reuse = regions[i].min + (f);
        if ( regions[i].border_high_reuse > regions[i].max ) {
            regions[i].border_high_reuse = regions[i].max;
        }
        regions[i].border_high_reuse += 8; // size of a ptr
        regions[i].border_moderate_reuse = regions[i].min + (2*f);
        if ( regions[i].border_moderate_reuse > regions[i].max ) {
            regions[i].border_moderate_reuse = regions[i].max;
        }
        regions[i].border_moderate_reuse += 8;
        total_used += (regions[i].border_high_reuse - regions[i].min - 8);
    }
    printf("\n");
    printf("%35s %20llu\n", "Cache Capcity (Bytes)", total_capacity);
    printf("%35s %20llu\n", "Used Capacity (Bytes)", total_used);
    printf("%35s %20f\n", "Utilization (Frac)", total_used * 100.0 / total_capacity);
    printf("\n");
}

#endif

