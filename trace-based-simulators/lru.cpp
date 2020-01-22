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

    printf("\n\n%35s %20s\n", "LRU is initialized!","");

    for ( uint64_t i = 0 ; i < total_addresses ; i++ ) {
        uint64_t& tr = trace[i];

        block_t* block;
        block_t* set;
        bool hit = find_block(tr, block_shift, sets, &block, &set);
        update_compute_stats(tr, hit);

        if ( !hit ) {
            // Miss
            int min = set[0].time;
            int min_index = 0;
            for ( int t = 1 ; t < assoc ; t++ ) {
                if ( set[t].time < min ) {
                    min = set[t].time;
                    min_index = t;
                }
            }
            //min_index is to be replaced
            set[min_index].time = i;
            set[min_index].addr = tr;
        } else {
            // Hit
            block->time = i;
        }
    }

    display();
    return 0;
}
