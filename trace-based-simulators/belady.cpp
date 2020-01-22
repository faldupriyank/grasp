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

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <map>
#include <cmath>

#include "cache.h"

using namespace std;

int main(int argc, char* argv[]) {

    init_cache(argc, argv);
    printf("\n\n%35s %20s\n\n", "Belady is initialized!", "");

    uint64_t* timestamp = (uint64_t*) malloc(total_addresses * 8);
    assert(timestamp);
    std::memset(timestamp, 1, sizeof(timestamp[0])*total_addresses);

    uint64_t* timestamp_set = (uint64_t*) malloc(total_addresses * 8);
    assert(timestamp_set);
    std::memset(timestamp_set, 1, sizeof(timestamp_set[0])*total_addresses);
    typedef std::map<uint64_t, uint64_t> map_t; 
    typedef map_t::iterator map_it;
    map_t time_map;

    for (int64_t i = total_addresses-1; i >=0 ; i-- ) {
        {
            map_it it;
            it = time_map.find(trace[i]);
            if ( it != time_map.end()  ) {
                timestamp[i] = it->second; 
                it->second = i;
            } else {
                time_map[trace[i]] = i;
            }
        }
    }

    printf("\n\n%35s %20s\n\n", "Belady future index is generated. Running second pass now.", "");

    for ( uint64_t i = 0 ; i < total_addresses ; i++ ) {
        uint64_t& tr = trace[i];

        block_t* block;
        block_t* set;
        bool hit = find_block(tr, block_shift, sets, &block, &set);
        update_compute_stats(tr, hit);

        if ( !hit ) {
            // Miss
            uint64_t max = set[0].time;
            int max_index = 0;
            for ( int t = 1 ; t < assoc ; t++ ) {
                if ( set[t].time > max ) {
                    max = set[t].time;
                    max_index = t;
                }
            }
            //max_index is to be replaced
            set[max_index].time = timestamp[i];
            set[max_index].addr = tr;
        } else {
            // Hit
            block->time = timestamp[i];
        }
    }

    display();
    return 0;
}



