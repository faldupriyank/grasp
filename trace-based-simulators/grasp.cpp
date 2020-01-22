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
#include <random>

#include "cache.h"

using namespace std;

extern int M_RRIP_INIT;
int main(int argc, char* argv[]) {

    const int num_bits_rrip = 3;
    const int M_RRIP = pow(2,num_bits_rrip)-1;
    const int I_RRIP = M_RRIP-1;
    const int P_RRIP = 1;
    const int H_RRIP = 0;
    M_RRIP_INIT = M_RRIP;

    init_cache(argc, argv);

    printf("\n\n%35s %20s\n", "Grasp is initialized!","");
    printf("%35s %20d\n", "Num bits per cache block", num_bits_rrip);
    printf("%35s %20d\n", "Worst Insertion Position", M_RRIP);
    printf("%35s %20d\n", "Intermediate Insertion Position", I_RRIP);
    printf("%35s %20d\n", "Priority Insertion Position", P_RRIP);
    printf("%35s %20d\n", "Priority Hit-Promotion Position", H_RRIP);

    for ( uint64_t i = 0 ; i < total_addresses ; i++ ) {
        uint64_t& tr = trace[i];

        block_t* block;
        block_t* set;
        bool hit = find_block(tr, block_shift, sets, &block, &set);
        update_compute_stats(tr, hit);

        bool in_high_range = is_in_high_reuse_region(tr, border_regions, num_border_regions);
        if ( !hit ) {
            int min_index = 0;
            int max_rrip = set[0].rrip;
            for ( int t = 1 ; t < assoc ; t++ ) {
                if ( set[t].rrip > max_rrip ) {
                    max_rrip = set[t].rrip; 
                    min_index = t;
                }
            }

            if ( max_rrip < M_RRIP ) {
                int diff = M_RRIP - max_rrip;
                for ( int t = 0 ; t < assoc ; t++ ) {
                    set[t].rrip += diff;
                    assert(set[t].rrip <= M_RRIP);
                }
            }

            if ( in_high_range ) {
                set[min_index].rrip = P_RRIP;
            } else {
                bool is_just_outside = is_in_moderate_reuse_region(tr, border_regions, num_border_regions);
                if ( is_just_outside ) {
                    set[min_index].rrip = I_RRIP;
                } else {
                    set[min_index].rrip = M_RRIP;
                }
            }
            set[min_index].addr = tr;
        } else {
            // hit
            if ( in_high_range ) {
                block->rrip = H_RRIP;
            } else {
                if ( block->rrip > 0 ) {
                    block->rrip--;
                }
            }
        }
    }

    display();
    return 0;
}
