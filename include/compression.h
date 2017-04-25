#pragma once

#define _XOPEN_SOURCE 500
#define _POSIX_C_SOURCE 200809L

#include <sys/stat.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <fcntl.h>
#include <unistd.h>

#include <Judy.h>

#include <bft/kseq.h>
KSEQ_INIT(int, read)

#include <bft/insertNode.h>
#include <bft/printMemory.h>
#include <bft/fasta.h>
#include <bft/replaceAnnotation.h>
#include <bft/write_to_disk.h>
#include <bft/marking.h>

#include <bft/list.h>
#include <bft/xxhash.h>

#include <bft/bft.h>
#include <bft/kseq.h>

#include <bft/intersection.h>

#define MAX_SIZE_LIST_PART 1000
#define DECOMP_NB_MAX_LOAD_PART 500

#define SIZE_BLOCK_BYTES 64
#define FOPEN_MAX_HARD 1024

typedef struct{
    int nb_hash;
    int seed_hash1;
    int seed_hash2;
    int nb_blocks;
    int size_block_bits;
    int size_block_bytes;
    uint64_t nb_bits_bf;
    uint64_t nb_bytes_bf;
    double fp_rate_bf;
    uint8_t* bf;
} BloomFilter;

typedef struct{
    int id_read;
    int position;
    char mismatch_char;
} Mismatch;

typedef struct{
    int occ_count;
    int nb_mismatches;
    Mismatch* list_mismatches;
    int64_t* id_occ;
    char* read;
} Read_count;

typedef struct{
    int position;
    int length;
    int rev;
    int occ_count;
    int nb_mismatches;
    Mismatch* list_mismatches;
    int64_t* id_occ;
} Pos_length_occ;

inline Read_count* create_Read_count();
inline Pos_length_occ* create_pos_length_occ();
inline void set_Pos_length_occ(Pos_length_occ* pos_length_occ, int pos, int length, int occ_count,
                               int rev, int nb_mismatches, int64_t* id_occ, Mismatch* list_mismatches);
inline Pos_length_occ* create_set_Pos_length_occ(int pos, int length, int occ_count, int rev, int nb_mismatches,
                                                 int64_t* id_occ, Mismatch* list_mismatches);
inline Pos_length_occ* create_set_Pos_length_occ_from_Read_count(int pos, int length, int rev,
                                                                 const Read_count* read_count);

inline Read_count* create_Read_count(){

    Read_count* read_count = calloc(1, sizeof(Read_count));
    ASSERT_NULL_PTR(read_count, "create_Read_count()\n");

    read_count->id_occ = NULL;
    read_count->list_mismatches = NULL;
    read_count->read = NULL;

    return read_count;
}

inline Pos_length_occ* create_pos_length_occ(){

    Pos_length_occ* pos_length_occ = calloc(1, sizeof(Pos_length_occ));
    ASSERT_NULL_PTR(pos_length_occ, "create_Read_count()\n");

    pos_length_occ->id_occ = NULL;
    pos_length_occ->list_mismatches = NULL;

    return pos_length_occ;
}

inline void set_Pos_length_occ(Pos_length_occ* pos_length_occ, int pos, int length, int occ_count,
                               int rev, int nb_mismatches, int64_t* id_occ, Mismatch* list_mismatches){

    ASSERT_NULL_PTR(pos_length_occ, "set_Pos_length_occ()\n");

    pos_length_occ->position = pos;
    pos_length_occ->length = length;
    pos_length_occ->rev = rev;
    pos_length_occ->occ_count = occ_count;
    pos_length_occ->id_occ = id_occ;
    pos_length_occ->list_mismatches = list_mismatches;
    pos_length_occ->nb_mismatches = nb_mismatches;

    return;
}

inline Pos_length_occ* create_set_Pos_length_occ(int pos, int length, int occ_count, int rev, int nb_mismatches,
                                                 int64_t* id_occ, Mismatch* list_mismatches){

    Pos_length_occ* pos_length_occ = malloc(sizeof(Pos_length_occ));
    ASSERT_NULL_PTR(pos_length_occ, "binning_reads()\n")

    pos_length_occ->position = pos;
    pos_length_occ->length = length;
    pos_length_occ->rev = rev;
    pos_length_occ->occ_count = occ_count;
    pos_length_occ->id_occ = id_occ;
    pos_length_occ->list_mismatches = list_mismatches;
    pos_length_occ->nb_mismatches = nb_mismatches;

    return pos_length_occ;
}

inline Pos_length_occ* create_set_Pos_length_occ_from_Read_count(int pos, int length, int rev, const Read_count* read_count){

    Pos_length_occ* pos_length_occ = malloc(sizeof(Pos_length_occ));
    ASSERT_NULL_PTR(pos_length_occ, "binning_reads()\n")

    pos_length_occ->position = pos;
    pos_length_occ->length = length;
    pos_length_occ->rev = rev;
    pos_length_occ->occ_count = read_count->occ_count;
    pos_length_occ->id_occ = read_count->id_occ;
    pos_length_occ->list_mismatches = read_count->list_mismatches;
    pos_length_occ->nb_mismatches = read_count->nb_mismatches;

    return pos_length_occ;
}

inline int custom_atoi(const char *old_str, char** new_str, char end) {

    int x = 0;

    while (*old_str != end) x = (x*10) + (*old_str++ - '0');

    *new_str = old_str;

    return x;
};

inline int custom_strstr_dna_triplet(const char *str, int len, bool upper_char, bool check_last_triplet_only){

    const uint32_t t_a = upper_char ? 0x414141 : 0x616161;
    const uint32_t t_c = upper_char ? 0x434343 : 0x636363;
    const uint32_t t_g = upper_char ? 0x474747 : 0x676767;
    const uint32_t t_t = upper_char ? 0x545454 : 0x747474;

    if (check_last_triplet_only){

        const uint32_t tmp = (str[len-3] << 16) | (str[len-2] << 8) | str[len-1];

        if ((tmp == t_a) || (tmp == t_c) || (tmp == t_g) || (tmp == t_t)) return len-3;
    }
    else{

        uint32_t tmp = (str[len-1] << 8) | str[len-2];

        for (int i = len-3; i >= 0; i--){

            tmp = ((tmp & 0xffff) << 8) | str[i];
            if ((tmp == t_a) || (tmp == t_c) || (tmp == t_g) || (tmp == t_t)) return i;
        }
    }

    return -1;
}


Pos_length_occ* parse_pos_length_occ(char** meta_to_parse, bool pair_ended);

uint32_t compress_FASTx_files(char* filename_read_1, char* filename_read_2, bool pair_ended,
                             int size_seed, int size_minimizer, int size_kmer, uint32_t prev_nb_part,
                             char* output_prefix, bool compress_shift, bool recycle_paths,
                             BFT_Root* root_no_iupac, BFT_Root* root_iupac);

void decompress_FASTx_file(char* output_prefix, BFT_Root* graph_no_iupac, BFT_Root* graph_iupac,
                           bool pair_ended, int size_kmer, int size_seed);

size_t insert_kmer_into_bf_from_graph(BFT_kmer* kmer, BFT* graph, va_list args);
size_t insert_kmer_into_bf(BloomFilter* bloom_filter, char* kmer, int length_kmer);

BloomFilter* create_BloomFilter(uint64_t nb_elem, double fp_rate_max);
void free_BloomFilter(BloomFilter* bloom_filter);

int64_t binning_reads(char* filename_reads_mate_1, char* filename_reads_mate_2, bool pair_ended, char* prefix_bin_name, int size_minimizer);
void create_super_reads(char* prefix_bin_name, bool pair_ended, int size_minimizer);
char** subdivide_bin(char* prefix_bin_name, uint64_t id_bin, int size_minimizer, bool delete_id_bin_file);
void reassemble_bins_to_super_bins(char* prefix_bin_name, uint64_t id_super_bin, char** filenames_bin, int size_minimizer);
uint64_t create_spanning_super_reads(char* prefix_bin_name, bool pair_ended, int64_t nb_reads_per_file, int size_minimizer);

int string_cmp(const void *a, const void *b);
int pos_length_occ_cmp(const void *a, const void *b);
int pos_length_occ_cmp_rev(const void *a, const void *b);
int mismatch_cmp(const void *a, const void *b);
int mismatch_cmp_rev(const void *a, const void *b);

void min_simReads_max_simkmers2(char* filename_read_1, char* filename_read_2, bool pair_ended,
                               char* output_prefix, int size_kmer, int size_seed, int size_minimizer,
                               bool compress_shift, BFT_Root* graph_no_iupac, BFT_Root* graph_iupac);

uint32_t compressKmers_from_KmerFiles_bis2(char* output_prefix, int size_seed, int size_kmer, uint32_t prev_nb_parts);
uint32_t compressKmers_from_KmerFiles_bis3(char* output_prefix, int size_seed, int size_kmer, uint32_t prev_nb_parts);

//void compress_annotations_BFT_disk(BFT_Root* bft, char* filename_bft);

void insert_kmers_partitions2(char* output_prefix, int size_seed, int size_kmer, bool recycle_paths,
                                  BFT_Root* root_no_IUPAC, BFT_Root* root_IUPAC);

void insert_partition(resultPresence* res, BFT_annotation* bft_annot, BFT* root, UC* uc, uint32_t part, uint32_t size_part_bytes);

void insert_kmers_list(List* list_kmers, BFT_annotation* bft_annot, BFT_Root* root_iupac, int lvl_root_iupac,
                       BFT_Root* root_no_iupac, int lvl_root_no_iupac, uint32_t part, uint32_t size_part_bytes);

void serialize_subpaths_recycling(char* filename_subpaths_recycling, Pvoid_t* recycl_subpaths, int size_seed);

uint32_t unserialize_subpaths_recycling(FILE* file_part_recycl2, Pvoid_t* recycled_parts, int size_seed,
                                    int64_t load_from_pos_start, int64_t nb_parts_to_load, bool seeds_already_loaded);

void sort_subpaths_recycling(Pvoid_t* recycled_parts, int size_seed, bool sort_and_replace);

uint32_t* recycle_sub_paths(char* seed, uint32_t recycled_part, Pvoid_t* recycl_subpaths);

void create_subgraph_decomp(char* output_prefix, FILE* file_partitions, Pvoid_t* prefix_part_recycl,
                            BFT_Root** graph_no_iupac, BFT_Root** graph_iupac, int size_seed);
size_t load_subgraph(BFT_kmer* kmer, BFT* graph, va_list args);
void cdbg_traverse_overlap_insert(BFT* bft_no_iupac, BFT* bft_iupac, BFT* insert_bft_no_iupac, BFT* insert_bft_iupac,
                                char* overlap, uint32_t partition, uint32_t size_partition);
size_t insert_kmer_with_partition(BFT_kmer* kmer, BFT* graph, va_list args);

bool cdbg_get_successor_overlap(BFT* bft, bool graph_is_iupac, char** overlap, uint32_t partition);
size_t kmer_has_partition(BFT_kmer* kmer, BFT* graph, va_list args);

void decompress2(char* output_prefix, BFT_Root* graph_no_iupac, BFT_Root* graph_iupac, bool pair_ended, int size_kmer, int size_seed);

bool init_extract_subseq_shifted(char* filename_shifts, FILE** file_shift);
uint8_t get_header_next_subseq_shifted(FILE* file_shift);
char* extract_next_subseq_shifted(FILE* file_shift, bool compress_shift, bool is_iupac, int length);
void deinit_extract_subseq_shifted(FILE** file_shift);

void extract_reads(char* output_prefix, bool paired_end);
void reorder_paired_reads(char* output_prefix);

void reorder_reads(char* filename_src, char* filename_dest, int size_minimizer, int size_window);

void write_meta_info_file(const char* output_prefix, bool paired, const char* version_tool);
bool read_meta_info_file(char* input_prefix);

void append_file(char* filename_base, char* filename_to_append);
void write_graph_info_file(char* graph, long int offset_graph_iupac, long int offset_info_file, uint32_t nb_parts,
                           int length_overlap, int pos_path_recycling, int threshold_path_recycling, const char* version_tool);
