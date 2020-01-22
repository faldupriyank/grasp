/*
 * @author Priyank Faldu <Priyank.Faldu@ed.ac.uk> <http://faldupriyank.com>
 *
 * Copyright 2019 The University of Edinburgh
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

#include <algorithm>
#include <cassert>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <parallel/algorithm>
#include <random>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h> 
#include <unistd.h>
#include <vector>

using namespace std;

//#define _DEBUG_OUTPUT_

//el --> edge list of the form (src, dst) in text file
//bel --> edge list of the form (src, dst) in binary file
//vgr --> binary csr format with no weight for edges
//cvgr --> binary csr format with no weight for edges (no self or redundant edges)
//csvgr --> binary csr format with no weight for edges (no self or redundant edges)
//         graph is symmetric -- so for every edge (u, v) there also exist (v, u)
//cintgr --> binary csr format with int weight for edges (no slef or redundant edges)

//edges are assumed to take 8 bytes and vertices (and edge weights) are assumed to take 4 bytes in binary file

//all *gr files contain a header of 24 bytes as follows:
//number of vertices
//number of edges
//major number
//minor number
//gr files are implemented based on the implementation from Galois.

string supported_lists[][2] =
    {
        { "el", "bel" },

        { "el", "vgr" },
        { "el", "cvgr" },
        { "el", "csvgr" },
        { "el", "cintgr" },

        { "bel", "vgr" },
        { "bel", "cvgr" },
        { "bel", "csvgr" },
        { "bel", "cintgr" },

        { "el", "print" },
        { "bel", "print" },
        { "vgr", "print" },
        { "cvgr", "print" },
        { "csvgr", "print" },
        { "cintgr", "print" },
    };


struct edge_t {
    uint32_t src;
    uint32_t dst;
    edge_t() {};
    edge_t(uint32_t s, uint32_t d): src(s), dst(d) {}
    edge_t(const edge_t& e) {
        src = e.src; dst = e.dst;
    }

    bool operator <(const edge_t& e) const { 
        return ((this->src < e.src) || (this->src == e.src && this->dst < e.dst));
    }

    bool operator ==(const edge_t& e) const { 
        return (this->src == e.src && this->dst == e.dst);
    }
};

struct edge_less {
    bool operator () (const edge_t& lhs, const edge_t& rhs) const {
        return (lhs.src < rhs.src);
    }
};

int read_el_file(const string& str, vector<edge_t>& edges, uint64_t& num_vert, uint64_t& num_edges) {

    edges.clear();
    ifstream ifs(str.c_str(), std::ifstream::in);
    if ( !ifs.good() ) {
        cout << "Unable to open : " << str << endl;
        return -1;
    }

    num_vert = 0;
    while (1) {
        uint32_t src, dst;
        ifs >> src >> dst;
        if (ifs.good()) {
            edges.push_back(edge_t(src, dst));
            if ( src > num_vert ) {
                num_vert = src;
            }
            if ( dst > num_vert ) {
                num_vert = dst;
            }
        } else {
            break;
        }
    }
    ifs.close();

    num_edges = edges.size();
    if (num_edges > 0) {
        num_vert++;
    }

    return 0;
}

char* get_file_mapping(const string& str, size_t& file_size) {
	struct stat sb;
	int fd = open(str.c_str(), O_RDONLY);
	if(fd < 0) {
		printf("Cannot open %s\n", str.c_str());
        return NULL;
	}
	fstat(fd, &sb);
	printf("#FileSize: %lu\n", (uint64_t)sb.st_size);
    file_size = sb.st_size;

	char* temp= (char*) mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
	assert(temp != NULL);

	close(fd);
	return temp;
}

int read_bel_file(const string& str, vector<edge_t>& edges, uint64_t& num_vert, uint64_t& num_edges) {

    edges.clear();

    size_t file_size = 0;
	edge_t* _edges = (edge_t*) get_file_mapping(str, file_size);
    assert((file_size % sizeof(edge_t)) == 0);
    
    num_edges = file_size/sizeof(edge_t);
	edges.resize(num_edges);
    num_vert = 0;

#pragma omp parallel for schedule(static) reduction(max: num_vert)
	for ( uint64_t i = 0 ; i < num_edges ; i++ ) {
		edges[i] = _edges[i];

        if ( edges[i].src > num_vert ) {
            num_vert = edges[i].src;
        }
        if ( edges[i].dst > num_vert ) {
            num_vert = edges[i].dst;
        }
	}

    if ( num_edges > 0 ) {
        num_vert++;
    }

    munmap(_edges, file_size);
    return 0;
}

int read_gr_file(const string& input_format, const string& input_file, uint64_t** _offsets, uint32_t** _dest_list, uint32_t** _weights, uint64_t** num_vert, uint64_t** num_edges) {
    size_t file_size = 0;
	char* temp = get_file_mapping(input_file, file_size);
	uint64_t* major;
	uint64_t* minor;

	major = (uint64_t*)(temp);
	temp += sizeof(uint64_t);

	minor = (uint64_t*)(temp);
	temp += sizeof(uint64_t);

	*num_vert = (uint64_t*)(temp);
	temp += sizeof(uint64_t);

	*num_edges = (uint64_t*)(temp);
	temp += sizeof(uint64_t);

    if ( input_format == "vgr" || input_format == "cvgr" || input_format == "csvgr" ) {
        assert(*major == 1);
        assert(*minor == 0);
    } else if ( input_format == "cintgr" ) {
        assert(*major == 1);
        assert(*minor == 4);
    }

	*_offsets = (uint64_t*)(temp);
	temp += ((**num_vert)*sizeof(uint64_t));

	*_dest_list = (uint32_t*)(temp);
	temp += ((**num_edges)*sizeof(uint32_t));
	if ( (**num_edges) % 2 == 1 ) {
		temp += sizeof(uint32_t);
	}

    if ( *major == 1 && *minor == 4 ) {
		*_weights = (uint32_t*)(temp);
    }
    
    return 0;
}

int generate_random_weights(uint32_t** _weigths, const uint64_t& num_edges, const uint32_t& min, const uint32_t& max) {
    uint32_t* weights = new uint32_t[num_edges];
    assert(weights);
    *_weigths = weights;

    default_random_engine generator;
    uniform_int_distribution<uint32_t> distribution(min,max);
    for ( uint64_t i = 0 ; i < num_edges ; i++ ) {
        weights[i] = distribution(generator);
    }

#ifdef _DEBUG_OUTPUT_
    cout << "Edge weights" << endl;
    for ( uint64_t i = 0 ; i < num_edges ; i++ ) {
        cout << i << "\t" << weights[i] << endl;
    }
#endif
    return 0;
}

int check_supported_in_out_format_pairs(const string& in, const string& out) {
    size_t len = sizeof(supported_lists) / sizeof(supported_lists[0]); 

    bool found = false;
    for ( size_t i = 0 ; i < len ; i++ ) {
        if ( supported_lists[i][0] == in && supported_lists[i][1] == out ) {
            found = true;
            break;
        }
    }

    if ( !found ) {
        cout << "Input format: " << in
            << " Output format: " << out
            << " - The combination of input/output format is not supported yet!"
            << endl;

        cout << "Total supported format conversion pairs: " << len << endl;
        for ( size_t i = 0 ; i < len ; i++ ) {
            cout << supported_lists[i][0] << " to " << supported_lists[i][1] << endl;
        }
        return -1;
    }

    return 0;
}

int make_symmetric_edges(const vector<edge_t>& edges, vector<edge_t>& symmetric_edges, const uint64_t& num_vert, const uint64_t& num_edges) {
    cout << "Adding symmetric edges" << endl;
    assert(edges.size() == num_edges);
    symmetric_edges.clear();
    symmetric_edges.resize(num_edges * 2);

#pragma omp parallel for
    for ( uint64_t i = 0 ; i < num_edges ; i++ ) {
        symmetric_edges[i] = edges[i];
        symmetric_edges[i+num_edges] = edge_t(edges[i].dst, edges[i].src);
	}

    return 0;
}

int make_edges_clean(const vector<edge_t>& edges, vector<edge_t>& cleaned_up_edges, const uint64_t& num_vert, const uint64_t& num_edges) {
    cout << "Cleaning up edges (removing self and redundant edges)" << endl;
    assert(edges.size() == num_edges);
    
    cleaned_up_edges.clear();

    edge_t* sorted_edges = new edge_t[num_edges];
    assert(sorted_edges);

#pragma omp parallel for
    for ( uint64_t i = 0 ; i < num_edges ; i++ ) {
        sorted_edges[i] = edges[i];
	}

    //sorts based on src and dst
	__gnu_parallel::sort<edge_t*>(&sorted_edges[0], &sorted_edges[num_edges]); 

#ifdef _DEBUG_OUTPUT_
    cout << endl;
    cout << "edges in original order" << endl;
	for ( uint64_t i = 0 ; i < num_edges ; i++ ) {
		cout << edges[i].src << "\t" << edges[i].dst << endl;
	}

    cout << "edges sorted by both source and desitnation" << endl;
	for ( uint64_t i = 0 ; i < num_edges ; i++ ) {
		cout << sorted_edges[i].src << "\t" << sorted_edges[i].dst << endl;
	}
    cout << endl;
#endif

    uint64_t num_self_edges = 0;
    uint64_t num_redundant_edges = 0;
    if ( sorted_edges[0].src != sorted_edges[0].dst ) {
        cleaned_up_edges.push_back(sorted_edges[0]);
    } else {
        num_self_edges++;
    }

    for ( uint64_t i = 1 ; i < num_edges ; i++ ) {
        if ( sorted_edges[i].src == sorted_edges[i].dst ) {
            num_self_edges++;
            continue;
        }

        if ( sorted_edges[i] == sorted_edges[i-1] ) {
            num_redundant_edges++;
            continue;
        }
        cleaned_up_edges.push_back(sorted_edges[i]);
    }
    assert(cleaned_up_edges.size() <= num_edges);

    cout << "Removed " << num_self_edges << " self_edges" << endl;
    cout << "Removed " << num_redundant_edges << " redundant_edges" << endl;


#ifdef _DEBUG_OUTPUT_
    cout << endl;

    cout << "edges after cleaned up" << endl;
	for ( uint64_t i = 0 ; i < cleaned_up_edges.size() ; i++ ) {
		cout << cleaned_up_edges[i].src << "\t" << cleaned_up_edges[i].dst << endl;
	}
#endif

    delete[] sorted_edges;
    return 0;
}

void print_edgelist_from_gr(const uint64_t* offsets, const uint32_t* dest_list, const uint32_t* weights, const uint64_t& num_vert, const uint64_t& num_edges) {

    cout << "offsets, size:" << num_vert << endl;
    for ( uint64_t v = 0 ; v < num_vert ; v++ ) {
        cout << v << " " << offsets[v] << endl;
    }

    cout << "dest_list, size:" << num_edges << endl;
    for ( uint64_t e = 0 ; e < num_edges ; e++ ) {
        cout << e << " " << dest_list[e] << endl;
    }

    uint64_t current_vert = 0;
    cout << endl << "edges" << endl;
    if ( weights ) {
        for ( uint64_t t = 0 ; t < offsets[0] && t < num_edges ; t++ ) {
            cout << 0 << " " << dest_list[t] << " " << weights[t] << endl;
        }
        for ( uint64_t i = 1 ; i < num_vert ; i++ ) {
            for ( uint64_t t = offsets[i-1] ; t < offsets[i] && t < num_edges ; t++ ) {
                cout << i << " " << dest_list[t] << " " << weights[t] << endl;
            }
        }
    } else {
        for ( uint64_t t = 0 ; t < offsets[0] && t < num_edges ; t++ ) {
            cout << 0 << " " << dest_list[t] << endl;
        }
        for ( uint64_t i = 1 ; i < num_vert ; i++ ) {
            for ( uint64_t t = offsets[i-1] ; t < offsets[i] && t < num_edges ; t++ ) {
                cout << i << " " << dest_list[t] << endl;
            }
        }
    }
}

int read_input_file(const string& input_format, const string& input_file, vector<edge_t>& edges, uint64_t& num_vert, uint64_t& num_edges) {

    num_vert = 0;
    num_edges = 0;
    edges.clear();

    int status = 0;
    cout << "Started reading from : " << input_file << endl;
    if ( input_format == "el" ) {
        status = read_el_file(input_file, edges, num_vert, num_edges);
        assert(edges.size() == num_edges);
    } else if ( input_format == "bel" ) {
        status = read_bel_file(input_file, edges, num_vert, num_edges);
        assert(edges.size() == num_edges);
	} else if ( input_format == "vgr" || input_format == "csvgr" || input_format == "cintgr" || input_format == "cvgr" ) {
		uint64_t* offsets = NULL;
        uint32_t* dest_list = NULL;
        uint32_t* weights = NULL;
        uint64_t* _num_vert = NULL;
        uint64_t* _num_edges = NULL;
        status = read_gr_file(input_format, input_file, &offsets, &dest_list, &weights, &_num_vert, &_num_edges);
        if ( status != 0 ) {
            return status;
        }
        assert(offsets != NULL);
        assert(dest_list != NULL);
        assert(_num_vert != NULL);
        assert(_num_edges != NULL);
        if ( input_format == "cintgr" ) {
            assert(weights != NULL);
        }
        num_vert = *_num_vert;
        num_edges = *_num_edges;
#ifdef _DEBUG_OUTPUT_
        print_edgelist_from_gr(offsets, dest_list, weights, num_vert, num_edges);
#endif
    } else {
        return -1;
    }

    cout << "Sucessfully read file with " << num_vert << " vertices and " << num_edges << " edges." << endl;
    return status;
}

int write_bel_file(const string& str, const vector<edge_t>& edges, const uint64_t& num_vert, const uint64_t& num_edges) {
    
    ofstream fout(str, ios::out | ios::binary);
    if (!fout.good()) {
        cout << "Unable to open : " << str << endl;
        return -1;
    }

    fout.write((char*)&edges[0], edges.size() * sizeof(edge_t));
    fout.close();

    return 0;
}

int convert_to_vgr(const vector<edge_t>& edges, const uint64_t& num_vert, const uint64_t& num_edges, uint64_t** _offsets, uint32_t** _dest_list) {

    edge_t* sorted_edges = new edge_t[num_edges];
    assert(sorted_edges);

    uint64_t* offsets = new uint64_t[num_vert+1];
    assert(offsets);
    *_offsets = offsets;

    uint32_t* dest_list = new uint32_t[num_edges];
    assert(dest_list);
    *_dest_list = dest_list;

#pragma omp parallel for
    for ( uint64_t i = 0 ; i < num_edges ; i++ ) {
        sorted_edges[i] = edges[i];
	}

	__gnu_parallel::stable_sort<edge_t*>(&sorted_edges[0], &sorted_edges[num_edges], edge_less());

    uint64_t v_id = 0;
    uint64_t e_id = 0;
    for ( ; e_id < num_edges; v_id++ ) {
        offsets[v_id] = e_id;
        if ( v_id == sorted_edges[e_id].src ) {
            for ( ; e_id < num_edges; e_id++ ) {
                if ( v_id == sorted_edges[e_id].src ) {
                    dest_list[e_id] = sorted_edges[e_id].dst;
                } else {
                    break;
                }
            }
        }
    }
    while ( v_id <= num_vert ) {
        offsets[v_id] = e_id;
        v_id++;
    }
    assert(offsets[0] == 0);
    assert(e_id == num_edges);

#ifdef _DEBUG_OUTPUT_
    cout << endl;

    cout << "edges in original order" << endl;
	for ( uint64_t i = 0 ; i < num_edges ; i++ ) {
		cout << edges[i].src << "\t" << edges[i].dst << endl;
	}

    cout << "edges sorted by only source" << endl;
	for ( uint64_t i = 0 ; i < num_edges ; i++ ) {
		cout << sorted_edges[i].src << "\t" << sorted_edges[i].dst << endl;
	}

    cout << endl << "prefix sum" << endl;
    for ( uint64_t i = 0 ; i <= num_vert ; i++ ) {
        cout << i << "\t" << offsets[i] << endl;
    }

    cout << endl << "dest list" << endl;
    for ( uint64_t i = 0 ; i < num_edges ; i++ ) {
        cout << i << "\t" << dest_list[i] << endl;
    }
    cout << endl;
#endif

    delete[] sorted_edges;
    return 0;
}

int write_gr_file(const string& str, uint64_t num_vert, uint64_t num_edges, uint64_t* offsets, uint32_t* dest_list, uint32_t* weights) {

    assert(offsets != NULL);
    assert(dest_list != NULL);

    if ( weights == NULL ) {
        cout << "Started writing vgr to : " << str << endl;
    } else {
        cout << "Started writing cintgr to : " << str << endl;
    }

    ofstream ofs(str, std::ios::binary|std::ofstream::trunc);
    uint64_t ver_major = 1;
    uint64_t ver_minor = (weights == NULL) ? 0 : 4;
    ofs.write((char*)&ver_major, 8);
    ofs.write((char*)&ver_minor, 8);
    ofs.write((char*)&num_vert, 8);
    ofs.write((char*)&num_edges, 8);
    ofs.write((char*)&offsets[1], sizeof(uint64_t)*(num_vert));
    ofs.write((char*)dest_list, sizeof(uint32_t)*num_edges);
    if ( num_edges % 2 == 1 ) {
        uint32_t x=0;
        ofs.write((char*)&x, sizeof(uint32_t));
    }

    if (weights != NULL ) {
        ofs.write((char*)weights, sizeof(uint32_t)*num_edges);
    }

    ofs.close();
    return 0;
}

int write_output_file(const string& output_format, const string& output_file, const vector<edge_t>& edges, const uint64_t& num_vert, const uint64_t& num_edges) {

    if (output_format != "print") {
        cout << "Started converting to " << output_format << " format." << endl;
    }
    int status = 0;
    if ( output_format == "bel" ) {
        status = write_bel_file(output_file, edges, num_vert, num_edges);
    } else if ( output_format == "vgr" ) {
        uint64_t* offsets = NULL;
        uint32_t* dest_list = NULL;
        status = convert_to_vgr(edges, num_vert, num_edges, &offsets, &dest_list);
        if (status != 0) {
            return status;
        }
        assert(offsets != NULL);
        assert(dest_list != NULL);

        status = write_gr_file(output_file, num_vert, num_edges, offsets, dest_list, NULL);
    } else if ( output_format == "cvgr" ) {
        vector<edge_t> cleaned_up_edges;
        status = make_edges_clean(edges, cleaned_up_edges, num_vert, num_edges);
        if ( status != 0 ) {
            return status;
        }

        uint64_t* offsets = NULL;
        uint32_t* dest_list = NULL;
        status = convert_to_vgr(cleaned_up_edges, num_vert, cleaned_up_edges.size(), &offsets, &dest_list);
        if (status != 0) {
            return status;
        }
        assert(offsets != NULL);
        assert(dest_list != NULL);

        status = write_gr_file(output_file, num_vert, cleaned_up_edges.size(), offsets, dest_list, NULL);
    } else if ( output_format == "csvgr" ) {
        vector<edge_t> symmetric_edges;
        status = make_symmetric_edges(edges, symmetric_edges, num_vert, num_edges);
        if ( status != 0 ) {
            return status;
        }

        vector<edge_t> cleaned_up_edges;
        status = make_edges_clean(symmetric_edges, cleaned_up_edges, num_vert, symmetric_edges.size());
        if ( status != 0 ) {
            return status;
        }
        symmetric_edges.clear();

        uint64_t* offsets = NULL;
        uint32_t* dest_list = NULL;
        status = convert_to_vgr(cleaned_up_edges, num_vert, cleaned_up_edges.size(), &offsets, &dest_list);
        if ( status != 0 ) {
            return status;
        }
        assert(offsets != NULL);
        assert(dest_list != NULL);

        status = write_gr_file(output_file, num_vert, cleaned_up_edges.size(), offsets, dest_list, NULL);
    } else if ( output_format == "cintgr" ) {
        vector<edge_t> cleaned_up_edges;
        status = make_edges_clean(edges, cleaned_up_edges, num_vert, num_edges);
        if ( status != 0 ) {
            return status;
        }

        uint64_t* offsets = NULL;
        uint32_t* dest_list = NULL;
        status = convert_to_vgr(cleaned_up_edges, num_vert, cleaned_up_edges.size(), &offsets, &dest_list);
        if ( status != 0 ) {
            return status;
        }
        assert(offsets != NULL);
        assert(dest_list != NULL);
        
        uint32_t* weights = NULL;
        generate_random_weights(&weights, cleaned_up_edges.size(), 1, 100);
        assert(weights != NULL);
        status = write_gr_file(output_file, num_vert, cleaned_up_edges.size(), offsets, dest_list, weights);
    } else if ( output_format == "print" ) {
        cout << "Number of vert: " << num_vert << " edges: " << num_edges << endl;
    } else {
        return -1;
    }

    if ( status == 0 && output_format != "print" ) {
        cout << "Successfully wrote to " << output_file << "." << endl;
    }
    return status;
}

int main(int argc, char* argv[]) {
    if ( argc != 4 ) {
        cout << "Usage: ./exe input_format output_format file(including full path, but without ext)" << endl;
        return -1;
    }

    string input_format = string(argv[1]);
    string output_format = string(argv[2]);
    string file = string(argv[3]);

    string input_file = file + "." + input_format;
    string output_file = file + "." + output_format;
   
    int status = check_supported_in_out_format_pairs(input_format, output_format);
    if ( status != 0 ) {
        return status;
    }
    cout << endl << "Starting conversion from " << input_format << " " << output_format << endl;
    uint64_t num_vert = 0;
    uint64_t num_edges = 0;
    vector<edge_t> edges;

    status = read_input_file(input_format, input_file, edges, num_vert, num_edges);
    if ( status != 0 ) {
        return status;
    }

    status = write_output_file(output_format, output_file, edges, num_vert, num_edges);
    if ( status != 0 ) {
        return status;
    }

    cout << endl;
    return 0;
}
