#include <libgen.h>

#include <jemalloc/jemalloc.h>
const char* malloc_conf = "narenas:1,tcache:false,lg_dirty_mult:8,lg_chunk:22";

#include <Judy.h>

#include <bft/printMemory.h>
#include <bft/write_to_disk.h>
#include <bft/snippets.h>
#include <bft/getRSS.h>

#include "fileio.h"
#include "lz_utils.h"
#include "compression.h"

#define MINIMIZER_MIN 8
#define MINIMIZER_MAX 15

#define OVERLAP_MIN 9
#define OVERLAP_MAX 15

void print_cli(){

    printf("Usage:\n");

    printf("\ndarrc <action parameter> <action-specific parameters> <general parameters>\n");

    printf("\n<action parameter>:\n");

    printf("\n\t-c\t[--compress]\t\tcompression\n");
	printf("\t-u\t[--update]\t\tupdate\n");
	printf("\t-d\t[--decompress]\t\tdecompression\n");
    printf("\n\t-v\t[--version]\t\tprint version info\n");
	printf("\t-h\t[--help]\t\tprint help info\n");

    printf("\n<compression parameters>:\n");

	printf("\n\t-k\t[--kmer]\targ\tlength of k-mers, must be either 18, 27, 36 (default), 45, 54 or 63\n");
	printf("\t-o\t[--overlap]\targ\tlength of k-mers overlap, must be between 8 and 15 (default: 11)\n");

	printf("\n<compression and update parameters>:\n");

	printf("\n\t-min\t[--minimizer]\targ\tlength of minimizers, must be between 8 and 15 (default: 9)\n");
	printf("\t-mis\t[--mismatch]\targ\tnumber of mismatches allowed during merging (default: 5)\n");
	printf("\t-1\t[--mate1]\targ\tinput FASTA/Q file: single-end reads or first mate of paired-end reads\n");
	printf("\t-2\t[--mate2]\targ\tinput FASTA/Q file: second mate of paired-end reads\n");
	printf("\t-l1\t[--listmate1]\targ\tlist of input FASTA/Q files: single-end reads or first mate of paired-end reads\n");
	printf("\t-l2\t[--listmate2]\targ\tlist of input FASTA/Q files: second mate of paired-end reads\n");

    printf("\n<decompression parameters>: None\n");

    printf("\n<general parameters>:\n");

	printf("\n\t-g\t[--graph]\targ\tgraph filename prefix (output for compression, input otherwise)\n");
	printf("\t-m\t[--meta]\targ\tmeta filename prefix (output for compression/update, input otherwise)\n");
	printf("\t-lm\t[--listmeta]\targ\tlist of meta filename prefixes (output for compression/update, input otherwise)\n");
	printf("\t-dir\t[--directory]\targ\tcompression and update:\tdirectory for temporary files\n");
	printf("\t\t\t\t\tdecompression:\t\tdirectory for temporary and output files\n");
	printf("\n");

    return;
}

int main(int argc, char *argv[])
{

    struct timeval tval_before, tval_after, tval_last, tval_result;
    gettimeofday(&tval_before, NULL);
    tval_last = tval_before;

    BFT_Root* root_no_iupac = NULL;
    BFT_Root* root_iupac = NULL;

    FILE* file_input = NULL;
    FILE* file_input_tmp = NULL;

    int i = 2;
    int nb_files_2_read_mate1 = -1;
    int nb_files_2_read_mate2 = -1;
    int nb_meta_prefix = -1;

    int length_kmer = 36;
    int length_overlap = 11;
    int length_minimizer = 9;

    int nb_mismatch = 5;

    int treshold_no_path_recycling = 10;
    int pos_trigger_path_recycling = 0;
    //int interval_comp_ratio = 20;

    int nb_meta_stream = 17;

    int len_tmp_graph, len_tmp_meta, len_tmp_str, z;

    int treshold_no_path_recycling_tmp = treshold_no_path_recycling;

    long int off_graph_iupac, off_info;

    uint32_t nb_parts = 0;

    //double min_inc_comp_ratio = 0.5 / 100;

    size_t size_buffer = SIZE_BUFFER;

    bool input_is_single_file = false;
    bool input_is_multiple_file = false;
    bool compress = false;
    bool update = false;
    bool pair_ended = false;
    bool recycle_paths = true;
    bool compress_shifts = true;

    char* tmp = NULL;
    char* buffer = NULL;
    char* graph_filename = NULL;
    char* tmp_dir = NULL;
    char* basename_meta = NULL;
    char* basename_meta_tmp = NULL;
    char* basename_graph = NULL;
    char* basename_graph_tmp = NULL;
    char* output_graph = NULL;
    char* output_meta = NULL;
    char* left_mate = NULL;
    char* right_mate = NULL;

    const char* version = "0.1";

    const char* header1 = "__/\\\\\\\\\\\\\\\\\\\\\\\\________/\\\\\\\\\\\\\\\\\\_______/\\\\\\\\\\\\\\\\\\________/\\\\\\\\\\\\\\\\\\____________/\\\\\\\\\\\\\\\\\\_        \n";
    const char* header2 = " _\\/\\\\\\////////\\\\\\____/\\\\\\\\\\\\\\\\\\\\\\\\\\___/\\\\\\///////\\\\\\____/\\\\\\///////\\\\\\_______/\\\\\\////////__       \n";
    const char* header3 = "  _\\/\\\\\\______\\//\\\\\\__/\\\\\\/////////\\\\\\_\\/\\\\\\_____\\/\\\\\\___\\/\\\\\\_____\\/\\\\\\_____/\\\\\\/___________      \n";
    const char* header4 = "   _\\/\\\\\\_______\\/\\\\\\_\\/\\\\\\_______\\/\\\\\\_\\/\\\\\\\\\\\\\\\\\\\\\\/____\\/\\\\\\\\\\\\\\\\\\\\\\/_____/\\\\\\_____________     \n";
    const char* header5 = "    _\\/\\\\\\_______\\/\\\\\\_\\/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\_\\/\\\\\\//////\\\\\\____\\/\\\\\\//////\\\\\\____\\/\\\\\\_____________    \n";
    const char* header6 = "     _\\/\\\\\\_______\\/\\\\\\_\\/\\\\\\/////////\\\\\\_\\/\\\\\\____\\//\\\\\\___\\/\\\\\\____\\//\\\\\\___\\//\\\\\\____________   \n";
    const char* header7 = "      _\\/\\\\\\_______/\\\\\\__\\/\\\\\\_______\\/\\\\\\_\\/\\\\\\_____\\//\\\\\\__\\/\\\\\\_____\\//\\\\\\___\\///\\\\\\__________  \n";
    const char* header8 = "       _\\/\\\\\\\\\\\\\\\\\\\\\\\\/___\\/\\\\\\_______\\/\\\\\\_\\/\\\\\\______\\//\\\\\\_\\/\\\\\\______\\//\\\\\\____\\////\\\\\\\\\\\\\\\\\\_ \n";
    const char* header9 = "        _\\////////////_____\\///________\\///__\\///________\\///__\\///________\\///________\\/////////__\n";

    const char* ext_meta = ".meta.darrc";
    const char* ext_graph = ".graph.darrc";

    const char* ext_graph_iupac = "_graph_iupac";
    const char* ext_graph_no_iupac = "_graph_no_iupac";

    char** meta_filenames = NULL;
    char** meta_prefix = NULL;

    char** paths_and_names_mate1 = malloc(sizeof(char*));
    ASSERT_NULL_PTR(paths_and_names_mate1, "main()\n")

    char** paths_and_names_mate2 = malloc(sizeof(char*));
    ASSERT_NULL_PTR(paths_and_names_mate2, "main()\n")

    paths_and_names_mate1[0] = NULL;
    paths_and_names_mate2[0] = NULL;

    printf("\n%s%s%s%s%s%s%s%s%s\n\n", header1, header2, header3, header4, header5, header6, header7, header8, header9);

    if (argc == 1) print_cli();
    else {

        if ((strcmp("-v", argv[1]) == 0) || (strcmp("--version", argv[1]) == 0)){
            printf("%s\n", version);
            exit(EXIT_SUCCESS);
        }
        else if ((strcmp("-h", argv[1]) == 0) || (strcmp("--help", argv[1]) == 0)){
            print_cli();
            exit(EXIT_SUCCESS);
        }

        buffer = calloc(size_buffer, sizeof(char));
        ASSERT_NULL_PTR(buffer, "main()");

        if ((strcmp("-c", argv[1]) == 0) || (strcmp("--compress", argv[1]) == 0) || (strcmp("-u", argv[1]) == 0) || (strcmp("--update", argv[1]) == 0)){

            compress = true;
            if ((strcmp("-u", argv[1]) == 0) || (strcmp("--update", argv[1]) == 0)) update =  true;

            for (; i < argc; i += 2){

                if ((strcmp("-k", argv[i]) == 0) || (strcmp("--kmer", argv[i]) == 0)){

                    if (update) ERROR("-k or --kmer is not a valid update parameter.\n");

                    length_kmer = atoi(argv[i+1]);

                    if (length_kmer <= 0) ERROR("Length k of k-mers must be superior 0\n")
                    else if (length_kmer > KMER_LENGTH_MAX) ERROR("Length k of k-mers must be inferior or equal to 63\n")
                    else if (length_kmer % NB_CHAR_SUF_PREF) ERROR("Length k of k-mers must be a multiple of 9\n");
                }
                else if ((strcmp("-o", argv[i]) == 0) || (strcmp("--overlap", argv[i]) == 0)){

                    if (update) ERROR("-o or --overlap is not a valid update parameter.\n");

                    length_overlap = atoi(argv[i+1]);

                    if (length_overlap < OVERLAP_MIN) ERROR("Length of overlap must be superior or equal to 9\n")
                    else if (length_overlap > OVERLAP_MAX) ERROR("Length of overlap must be inferior or equal to 15\n")
                }
                else if ((strcmp("-min", argv[i]) == 0) || (strcmp("--minimizer", argv[i]) == 0)){

                    length_minimizer = atoi(argv[i+1]);

                    if (length_minimizer < MINIMIZER_MIN) ERROR("Length of minimizer must be superior or equal to 8\n")
                    else if (length_minimizer > MINIMIZER_MAX) ERROR("Length of minimizer must be inferior or equal to 15\n")
                }
                else if ((strcmp("-mis", argv[i]) == 0) || (strcmp("--mismatch", argv[i]) == 0)){

                    if ((nb_mismatch = atoi(argv[i+1])) < 0) ERROR("Number of mismatches must be superior or equal to 0\n")
                }
                else if ((strcmp("-1", argv[i]) == 0) || (strcmp("--mate1", argv[i]) == 0)){

                    if (input_is_multiple_file) ERROR("-1 parameter cannot be used with -l1 and/or -l2 parameter.\n");

                    nb_files_2_read_mate1 = 1;
                    input_is_single_file = true;

                    file_input = fopen(argv[i+1], "r");
                    ASSERT_NULL_PTR(file_input,"Cannot open input file for mate 1.\n")

                    fclose(file_input);

                    paths_and_names_mate1[0] = malloc((strlen(argv[i+1]) + 1) * sizeof(char));
                    ASSERT_NULL_PTR(paths_and_names_mate1[0], "main()\n")

                    strcpy(paths_and_names_mate1[0], argv[i+1]);
                }
                else if ((strcmp("-2", argv[i]) == 0) || (strcmp("--mate2", argv[i]) == 0)){

                    if (input_is_multiple_file) ERROR("-2 parameter cannot be used with -l1 and/or -l2 parameter.\n");

                    nb_files_2_read_mate2 = 1;
                    input_is_single_file = true;

                    file_input = fopen(argv[i+1], "r");
                    ASSERT_NULL_PTR(file_input, "Cannot open input file for mate 2.\n")

                    fclose(file_input);

                    paths_and_names_mate2[0] = malloc((strlen(argv[i+1]) + 1) * sizeof(char));
                    ASSERT_NULL_PTR(paths_and_names_mate2[0],"main()\n")

                    strcpy(paths_and_names_mate2[0], argv[i+1]);
                }
                else if ((strcmp("-l1", argv[i]) == 0) || (strcmp("--listmate1", argv[i]) == 0)){

                    if (input_is_single_file) ERROR("-l1 parameter cannot be used with -1 and/or -2 parameter.\n");

                    nb_files_2_read_mate1 = 0;
                    input_is_multiple_file = true;

                    file_input = fopen(argv[i+1], "r");
                    ASSERT_NULL_PTR(file_input, "Cannot open input list of files for mate 1.\n")

                    while (getline(&buffer, &size_buffer, file_input) != -1){

                        if (buffer[strlen(buffer) - 1] == '\n') buffer[strlen(buffer) - 1] = '\0';

                        if (strlen(buffer)){

                            if ((file_input_tmp = fopen(buffer, "r")) == NULL){
                                fprintf(stderr,"Cannot open file %s specified on line %d of the list of files for mate 1.\n", buffer, nb_files_2_read_mate1);
                                exit(EXIT_FAILURE);
                            }

                            fclose(file_input_tmp);
                        }
                        else if ((nb_files_2_read_mate2 != -1) && ((nb_files_2_read_mate1 + 1) > nb_files_2_read_mate2)){
                            ERROR("Different number of lines in the two mate file lists.\n");
                        }
                        else if ((nb_files_2_read_mate2 != -1) && (paths_and_names_mate2[nb_files_2_read_mate1]) == NULL){
                            fprintf(stderr,"No file specified on line %d of the lists of mate files.\n", nb_files_2_read_mate1);
                            exit(EXIT_FAILURE);
                        }

                        if ((nb_files_2_read_mate2 != -1) && ((nb_files_2_read_mate1 + 1) > nb_files_2_read_mate2)){
                            ERROR("Different number of lines in the two mate file lists.\n");
                        }

                        paths_and_names_mate1 = realloc(paths_and_names_mate1, (nb_files_2_read_mate1 + 1) * sizeof(char*));
                        ASSERT_NULL_PTR(paths_and_names_mate1, "main()\n")

                        if (strlen(buffer)){

                            paths_and_names_mate1[nb_files_2_read_mate1] = malloc((strlen(buffer) + 1) * sizeof(char));
                            ASSERT_NULL_PTR(paths_and_names_mate1[nb_files_2_read_mate1], "main()\n")

                            strcpy(paths_and_names_mate1[nb_files_2_read_mate1], buffer);
                        }
                        else paths_and_names_mate1[nb_files_2_read_mate1] = NULL;

                        nb_files_2_read_mate1++;
                    }

                    fclose(file_input);
                }
                else if ((strcmp("-l2", argv[i]) == 0) || (strcmp("--listmate2", argv[i]) == 0)){

                    if (input_is_single_file) ERROR("-l2 parameter cannot be used with -1 and/or -2 parameter.\n");

                    nb_files_2_read_mate2 = 0;
                    input_is_multiple_file = true;

                    file_input = fopen(argv[i+1], "r");
                    ASSERT_NULL_PTR(file_input, "Cannot open input list of files for mate 2.\n")

                    while (getline(&buffer, &size_buffer, file_input) != -1){

                        if (buffer[strlen(buffer) - 1] == '\n') buffer[strlen(buffer) - 1] = '\0';

                        if (strlen(buffer)){

                            if ((file_input_tmp = fopen(buffer, "r")) == NULL){
                                fprintf(stderr,"Cannot open file %s specified on line %d of the list of files for mate 2.\n", buffer, nb_files_2_read_mate2);
                                exit(EXIT_FAILURE);
                            }

                            fclose(file_input_tmp);
                        }
                        else if ((nb_files_2_read_mate1 != -1) && ((nb_files_2_read_mate2 + 1) > nb_files_2_read_mate1)){
                            ERROR("Different number of lines in the two mate file lists.\n");
                        }
                        else if ((nb_files_2_read_mate1 != -1) && (paths_and_names_mate1[nb_files_2_read_mate2]) == NULL){
                            fprintf(stderr,"No file specified on line %d of the lists of mate files.\n", nb_files_2_read_mate2);
                            exit(EXIT_FAILURE);
                        }

                        if ((nb_files_2_read_mate1 != -1) && ((nb_files_2_read_mate2 + 1) > nb_files_2_read_mate1))
                            ERROR("Different number of lines in the two mate file lists.\n");

                        paths_and_names_mate2 = realloc(paths_and_names_mate2, (nb_files_2_read_mate2 + 1) * sizeof(char*));
                        ASSERT_NULL_PTR(paths_and_names_mate1, "main()\n")

                        if (strlen(buffer)){

                            paths_and_names_mate2[nb_files_2_read_mate2] = malloc((strlen(buffer) + 1) * sizeof(char));
                            ASSERT_NULL_PTR(paths_and_names_mate2[nb_files_2_read_mate2], "main()\n")

                            strcpy(paths_and_names_mate2[nb_files_2_read_mate2], buffer);
                        }
                        else paths_and_names_mate2[nb_files_2_read_mate2] = NULL;

                        nb_files_2_read_mate2++;
                    }

                    fclose(file_input);
                }
                else break;
            }
        }
        else if ((strcmp("-d", argv[1]) == 0) || (strcmp("--decompress", argv[1]) == 0)){

        }
        else {
            fprintf(stderr, "Unrecognized first parameter used: %s.\n", argv[1]);
            exit(EXIT_FAILURE);
        }

        for (; i < argc; i += 2){

            if ((strcmp("-m", argv[i]) == 0) || (strcmp("--meta", argv[i]) == 0)){

                if (input_is_multiple_file) ERROR("-m parameter cannot be used with -l1 and/or -l2 parameter, use -lm instead.\n");

                if (strlen(argv[i+1]) == 0) ERROR("A meta filename prefix must be specified after parameter -m.\n")

                nb_meta_prefix = 1;

                meta_prefix = malloc(sizeof(char*));
                ASSERT_NULL_PTR(meta_prefix, "main()\n");

                meta_prefix[0] = malloc((strlen(argv[i+1]) + strlen(ext_meta) + 1) * sizeof(char));
                ASSERT_NULL_PTR(meta_prefix[0], "main()\n");

                strcpy(meta_prefix[0], argv[i+1]);
                strcpy(&(meta_prefix[0][strlen(meta_prefix[0])]), ext_meta);

                if (compress){

                    file_input = fopen(meta_prefix[0], "wb");
                    ASSERT_NULL_PTR(file_input, "Cannot create file starting with meta prefix.\n");

                    if (fwrite(buffer, sizeof(char), 1, file_input) != 1) ERROR("Cannot write in file starting with meta prefix\n");

                    fclose(file_input);

                    if (remove(meta_prefix[0]) == -1) ERROR("Cannot delete file starting with meta prefix.\n")
                }
                else {

                    file_input = fopen(meta_prefix[0], "rb");
                    ASSERT_NULL_PTR(file_input, "Cannot read file starting with meta prefix.\n");

                    if (fread(buffer, sizeof(uint8_t), 1, file_input) != 1) ERROR("Cannot read file starting with meta prefix\n");

                    fclose(file_input);
                }
            }
            else if ((strcmp("-lm", argv[i]) == 0) || (strcmp("--list_meta", argv[i]) == 0)){

                if (input_is_single_file) ERROR("-lm parameter cannot be used with -1 and/or -2 parameter, use -m instead.\n");

                nb_meta_prefix = 0;

                file_input = fopen(argv[i+1], "r");
                ASSERT_NULL_PTR(file_input, "Cannot open input list of meta filename prefixes.\n")

                while (getline(&buffer, &size_buffer, file_input) != -1){

                    if (buffer[strlen(buffer) - 1] == '\n') buffer[strlen(buffer) - 1] = '\0';
                    strcpy(&buffer[strlen(buffer)], ext_meta);

                    if (!compress){

                        if ((file_input_tmp = fopen(buffer, "rb")) == NULL){
                            fprintf(stderr,"Cannot open file %s with prefix specified on line %d of the list of meta filename prefixes.\n", buffer, nb_meta_prefix);
                            exit(EXIT_FAILURE);
                        }

                        fclose(file_input_tmp);
                    }
                    else if (strlen(argv[i+1]) == 0){
                        fprintf(stderr,"Empty line on line %d of the list of meta filename prefixes.\n", nb_meta_prefix);
                        exit(EXIT_FAILURE);
                    }

                    meta_prefix = realloc(meta_prefix, (nb_meta_prefix + 1) * sizeof(char*));
                    ASSERT_NULL_PTR(meta_prefix, "main()\n")

                    meta_prefix[nb_meta_prefix] = malloc((strlen(buffer) + strlen(ext_meta) + 1) * sizeof(char));
                    ASSERT_NULL_PTR(meta_prefix[nb_meta_prefix], "main()\n")

                    strcpy(meta_prefix[nb_meta_prefix], buffer);

                    nb_meta_prefix++;
                }

                fclose(file_input);
            }
            else if ((strcmp("-g", argv[i]) == 0) || (strcmp("--graph", argv[i]) == 0)){

                if (strlen(argv[i+1]) == 0) ERROR("A graph filename prefix must be specified after parameter -g.\n")

                graph_filename = malloc((strlen(argv[i+1]) + strlen(ext_graph) + 1) * sizeof(char));
                ASSERT_NULL_PTR(graph_filename, "main()\n");

                strcpy(graph_filename, argv[i+1]);
                strcpy(&graph_filename[strlen(graph_filename)], ext_graph);

                if (compress){

                    file_input = fopen("test_file", "wb");
                    ASSERT_NULL_PTR(file_input, "Cannot create file starting with graph prefix.\n");

                    if (fwrite(buffer, sizeof(char), 1, file_input) != 1) ERROR("Cannot write in file starting with graph prefix\n");

                    fclose(file_input);

                    if (remove("test_file") == -1) ERROR("Cannot delete file starting with graph prefix.\n")
                }
                else {

                    file_input = fopen(graph_filename, "rb");
                    ASSERT_NULL_PTR(file_input, "Cannot read file starting with graph prefix.\n");

                    if (fread(buffer, sizeof(uint8_t), 1, file_input) != 1) ERROR("Cannot read file starting with graph prefix\n");

                    fclose(file_input);
                }
            }
            else if ((strcmp("-dir", argv[i]) == 0) || (strcmp("--directory", argv[i]) == 0)){

                len_tmp_str = strlen(argv[i+1]);

                if (len_tmp_str == 0) ERROR("A directory must be specified after parameter -dir.\n")

                tmp_dir = malloc((len_tmp_str + 2) * sizeof(char));
                ASSERT_NULL_PTR(tmp_dir, "main()\n");

                strcpy(tmp_dir, argv[i+1]);

                if (tmp_dir[len_tmp_str - 1] != '/'){
                    tmp_dir[len_tmp_str] = '/';
                    tmp_dir[len_tmp_str + 1] = '\0';
                }

                strcpy(buffer, tmp_dir);
                strcpy(&buffer[strlen(buffer)], "test_file");

                file_input = fopen(buffer, "wb");
                ASSERT_NULL_PTR(file_input, "Cannot create file in directory.\n");

                if (fwrite(buffer, sizeof(char), 1, file_input) != 1) ERROR("Cannot write in directory\n");

                fclose(file_input);

                if (remove(buffer) == -1) ERROR("Cannot delete file in directory.\n")
            }
        }

        free(buffer);

        if (tmp_dir == NULL){

            tmp_dir = malloc(sizeof(char));
            ASSERT_NULL_PTR(tmp_dir, "main()\n");

            tmp_dir[0] = '\0';
        }

        if (graph_filename == NULL) ERROR("A graph prefix name as input or output is required (-g).\n");
        if (meta_prefix == NULL) ERROR("A meta prefix name as input or output is required (-m).\n");
        if (compress && (nb_files_2_read_mate1 <= 0)) ERROR("An input file or list of input files is required (-1 or -l1) to compress or update.\n");

        if (compress){

            file_input = fopen(graph_filename, "rb");

            if (file_input != NULL){

                fclose(file_input);

                basename_graph_tmp = malloc((strlen(graph_filename) + 1) * sizeof(char));
                ASSERT_NULL_PTR(basename_graph_tmp, "main()")

                strcpy(basename_graph_tmp, graph_filename);
                basename_graph = basename(basename_graph_tmp);
                tmp = strchr(basename_graph, '.');
                if (tmp != NULL) *tmp = '\0';
                len_tmp_graph = strlen(tmp_dir) + strlen(basename_graph);

                output_graph = malloc((len_tmp_graph + 50) * sizeof(char));
                ASSERT_NULL_PTR(output_graph, "main()")

                strcpy(output_graph, tmp_dir);
                strcpy(&output_graph[strlen(output_graph)], basename_graph);
                strcpy(&output_graph[strlen(output_graph)], ext_graph_no_iupac);

                FIO_decompressFilename(output_graph, graph_filename, NULL);

                root_no_iupac = read_BFT_Root(output_graph);

                file_input = fopen(output_graph, "rb");
                ASSERT_NULL_PTR(file_input, "main()\n");

                fseek(file_input, 0 - ((long int) (2 * (sizeof(long int) + sizeof(int)) + sizeof(int) + sizeof(uint32_t))), SEEK_END);

                fread(&nb_parts, sizeof(uint32_t), 1, file_input);
                fread(&length_overlap, sizeof(int), 1, file_input);
                fread(&off_graph_iupac, sizeof(long int), 1, file_input);
                fread(&off_info, sizeof(long int), 1, file_input);
                fread(&pos_trigger_path_recycling, sizeof(int), 1, file_input);
                fread(&treshold_no_path_recycling_tmp, sizeof(int), 1, file_input);

                fclose(file_input);

                root_iupac = read_BFT_Root_offset(output_graph, off_graph_iupac);
                if (remove(output_graph)) printf("Warning: Could not remove temporary file.\n");

                free(basename_graph_tmp);
                free(output_graph);
            }

            for (z = 0; z < MAX(nb_files_2_read_mate1, nb_files_2_read_mate2); z++){

                pair_ended = false;
                left_mate = NULL;
                right_mate = NULL;

                if ((paths_and_names_mate1 == NULL) || (paths_and_names_mate1[z] == NULL)) left_mate = paths_and_names_mate2[z];
                else {

                    left_mate = paths_and_names_mate1[z];

                    if ((paths_and_names_mate2 != NULL) && (paths_and_names_mate2[z] != NULL)){

                        pair_ended = true;
                        right_mate = paths_and_names_mate2[z];
                    }
                }

                printf("\nCompressing %s", left_mate);
                if (pair_ended) printf(" and %s", right_mate);
                printf("\n");

                if (z){

                    basename_graph_tmp = malloc((strlen(graph_filename) + 1) * sizeof(char));
                    ASSERT_NULL_PTR(basename_graph_tmp, "main()")

                    strcpy(basename_graph_tmp, graph_filename);
                    basename_graph = basename(basename_graph_tmp);
                    tmp = strchr(basename_graph, '.');
                    if (tmp != NULL) *tmp = '\0';
                    len_tmp_graph = strlen(tmp_dir) + strlen(basename_graph);

                    output_graph = malloc((len_tmp_graph + 50) * sizeof(char));
                    ASSERT_NULL_PTR(output_graph, "main()")

                    strcpy(output_graph, tmp_dir);
                    strcpy(&output_graph[strlen(output_graph)], basename_graph);
                    strcpy(&output_graph[strlen(output_graph)], ext_graph_no_iupac);

                    root_no_iupac = read_BFT_Root(output_graph);
                    if (remove(output_graph)) printf("Warning: Could not remove temporary file.\n");

                    strcpy(output_graph, tmp_dir);
                    strcpy(&output_graph[strlen(output_graph)], basename_graph);
                    strcpy(&output_graph[strlen(output_graph)], ext_graph_iupac);

                    root_iupac = read_BFT_Root(output_graph);
                    if (remove(output_graph)) printf("Warning: Could not remove temporary file.\n");

                    free(basename_graph_tmp);
                    free(output_graph);
                }

                basename_meta_tmp = malloc((strlen(meta_prefix[z]) + 1) * sizeof(char));
                ASSERT_NULL_PTR(basename_meta_tmp, "main()")

                meta_filenames = malloc(nb_meta_stream * sizeof(char*));
                ASSERT_NULL_PTR(meta_filenames, "main()")

                strcpy(basename_meta_tmp, meta_prefix[z]);
                basename_meta = basename(basename_meta_tmp);
                tmp = strchr(basename_meta, '.');
                if (tmp != NULL) *tmp = '\0';
                len_tmp_meta = strlen(tmp_dir) + strlen(basename_meta);

                output_meta = malloc((len_tmp_meta + 50) * sizeof(char));
                ASSERT_NULL_PTR(output_meta, "main()")

                strcpy(output_meta, tmp_dir);
                strcpy(&output_meta[strlen(tmp_dir)], basename_meta);

                if (z == pos_trigger_path_recycling){

                    nb_parts = compress_FASTx_files(left_mate, right_mate, pair_ended,
                                                    length_overlap, length_minimizer, length_kmer, nb_parts, output_meta,
                                                    compress_shifts, !recycle_paths, root_no_iupac, root_iupac);

                    treshold_no_path_recycling_tmp *= 2;
                    pos_trigger_path_recycling += treshold_no_path_recycling_tmp;
                }
                else{

                    nb_parts = compress_FASTx_files(left_mate, right_mate, pair_ended,
                                                    length_overlap, length_minimizer, length_kmer, nb_parts, output_meta,
                                                    compress_shifts, recycle_paths, root_no_iupac, root_iupac);
                }

                write_meta_info_file(output_meta, pair_ended, version);

                printf("[LZMA compressing meta-data]\n\n");

                for (int l = 0; l < nb_meta_stream; l++){

                    meta_filenames[l] = malloc((len_tmp_meta + 50) * sizeof(char));
                    ASSERT_NULL_PTR(meta_filenames[l], "main()")

                    strcpy(meta_filenames[l], output_meta);
                }

                strcpy(&(meta_filenames[0][len_tmp_meta]), "_info");
                strcpy(&(meta_filenames[1][len_tmp_meta]), "_runs");
                strcpy(&(meta_filenames[2][len_tmp_meta]), "_partitions");
                strcpy(&(meta_filenames[3][len_tmp_meta]), "_read_info");
                strcpy(&(meta_filenames[4][len_tmp_meta]), "_recycl_parts1");
                strcpy(&(meta_filenames[5][len_tmp_meta]), "_recycl_parts2");
                strcpy(&(meta_filenames[6][len_tmp_meta]), "_shifts");
                strcpy(&(meta_filenames[7][len_tmp_meta]), "_span_sup_reads_pos");
                strcpy(&(meta_filenames[8][len_tmp_meta]), "_span_sup_reads_length");
                strcpy(&(meta_filenames[9][len_tmp_meta]), "_span_sup_reads_rev");
                strcpy(&(meta_filenames[10][len_tmp_meta]), "_span_sup_reads_occ");
                strcpy(&(meta_filenames[11][len_tmp_meta]), "_span_sup_reads_pos_mis");
                strcpy(&(meta_filenames[12][len_tmp_meta]), "_span_sup_reads_char_mis");
                strcpy(&(meta_filenames[13][len_tmp_meta]), "_span_sup_reads_id_mis");
                strcpy(&(meta_filenames[14][len_tmp_meta]), "_span_sup_reads_nb_mis");
                strcpy(&(meta_filenames[15][len_tmp_meta]), "_span_sup_reads_mate");
                strcpy(&(meta_filenames[16][len_tmp_meta]), "_span_sup_reads_mate_info");

                encode_lzma(meta_filenames, nb_meta_stream - (!pair_ended * 2), meta_prefix[z], COMP_LVL_MAX);

                for (int l = 0; l < nb_meta_stream - (!pair_ended * 2); l++)
                    if (remove(meta_filenames[l])) printf("Warning: Could not remove temporary file.\n");

                basename_graph_tmp = malloc((strlen(graph_filename) + 1) * sizeof(char));
                ASSERT_NULL_PTR(basename_graph_tmp, "main()")

                strcpy(basename_graph_tmp, graph_filename);
                basename_graph = basename(basename_graph_tmp);
                tmp = strchr(basename_graph, '.');
                if (tmp != NULL) *tmp = '\0';
                len_tmp_graph = strlen(tmp_dir) + strlen(basename_graph);

                output_graph = malloc((len_tmp_graph + 50) * sizeof(char));
                ASSERT_NULL_PTR(output_graph, "main()")

                strcpy(output_graph, tmp_dir);
                strcpy(&output_graph[strlen(output_graph)], basename_graph);
                strcpy(&output_graph[strlen(output_graph)], ext_graph_no_iupac);

                strcpy(&(meta_filenames[0][len_tmp_meta]), ext_graph_no_iupac);

                rename(meta_filenames[0], output_graph);

                strcpy(output_graph, tmp_dir);
                strcpy(&output_graph[strlen(output_graph)], basename_graph);
                strcpy(&output_graph[strlen(output_graph)], ext_graph_iupac);

                strcpy(&(meta_filenames[0][len_tmp_meta]), ext_graph_iupac);

                rename(meta_filenames[0], output_graph);

                free(basename_graph_tmp);
                free(output_graph);

                for (int l = 0; l < nb_meta_stream; l++) free(meta_filenames[l]);

                gettimeofday(&tval_after, NULL);
                time_spent(&tval_last, &tval_after, &tval_result);

                printf("Elapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
                printf("Peak of memory: %llu mb\n\n", ((unsigned long long int)getPeakRSS())/1024);

                tval_last = tval_after;
            }

            printf("[ZSTD compressing G-DBG]\n\n");

            basename_graph_tmp = malloc((strlen(graph_filename) + 1) * sizeof(char));
            ASSERT_NULL_PTR(basename_graph_tmp, "main()")

            strcpy(basename_graph_tmp, graph_filename);
            basename_graph = basename(basename_graph_tmp);
            tmp = strchr(basename_graph, '.');
            if (tmp != NULL) *tmp = '\0';
            len_tmp_graph = strlen(tmp_dir) + strlen(basename_graph);

            output_graph = malloc((len_tmp_graph + 50) * sizeof(char));
            ASSERT_NULL_PTR(output_graph, "main()")

            strcpy(output_graph, tmp_dir);
            strcpy(&output_graph[strlen(output_graph)], basename_graph);
            strcpy(&output_graph[strlen(output_graph)], ext_graph_no_iupac);

            char* output_graph_append = malloc((len_tmp_graph + 50) * sizeof(char));
            ASSERT_NULL_PTR(output_graph_append, "main()")

            strcpy(output_graph_append, tmp_dir);
            strcpy(&output_graph_append[strlen(output_graph_append)], basename_graph);
            strcpy(&output_graph_append[strlen(output_graph_append)], ext_graph_iupac);

            file_input = fopen(output_graph, "ab");
            ASSERT_NULL_PTR(file_input, "main()")

            off_graph_iupac = ftell(file_input);

            fclose(file_input);

            append_file(output_graph, output_graph_append);
            if (remove(output_graph_append)) printf("Warning: Could not remove temporary file.\n");

            file_input = fopen(output_graph, "ab");
            ASSERT_NULL_PTR(file_input, "main()")

            off_info = ftell(file_input);

            fclose(file_input);

            write_graph_info_file(output_graph, off_graph_iupac, off_info, nb_parts, length_overlap,
                                  pos_trigger_path_recycling, treshold_no_path_recycling_tmp, version);

            remove(graph_filename);

            FIO_compressFilename(graph_filename, output_graph, NULL, 1);
            if (remove(output_graph)) printf("Warning: Could not remove temporary file.\n");

            free(output_graph_append);
            free(basename_graph_tmp);
            free(output_graph);
        }
        else{

            basename_graph_tmp = malloc((strlen(graph_filename) + 1) * sizeof(char));
            ASSERT_NULL_PTR(basename_graph_tmp, "main()")

            strcpy(basename_graph_tmp, graph_filename);
            basename_graph = basename(basename_graph_tmp);
            tmp = strchr(basename_graph, '.');
            if (tmp != NULL) *tmp = '\0';
            len_tmp_graph = strlen(tmp_dir) + strlen(basename_graph);

            output_graph = malloc((len_tmp_graph + 50) * sizeof(char));
            ASSERT_NULL_PTR(output_graph, "main()")

            strcpy(output_graph, tmp_dir);
            strcpy(&output_graph[strlen(output_graph)], basename_graph);
            strcpy(&output_graph[strlen(output_graph)], ext_graph_no_iupac);

            FIO_decompressFilename(output_graph, graph_filename, NULL);

            root_no_iupac = read_BFT_Root(output_graph);

            file_input = fopen(output_graph, "rb");
            ASSERT_NULL_PTR(file_input, "main()\n");

            fseek(file_input, 0 - ((long int) (2 * (sizeof(long int) + sizeof(int)) + sizeof(int) + sizeof(uint32_t))), SEEK_END);

            fread(&nb_parts, sizeof(uint32_t), 1, file_input);
            fread(&length_overlap, sizeof(int), 1, file_input);
            fread(&off_graph_iupac, sizeof(long int), 1, file_input);

            fclose(file_input);

            root_iupac = read_BFT_Root_offset(output_graph, off_graph_iupac);
            if (remove(output_graph)) printf("Warning: Could not remove temporary file.\n");

            free(basename_graph_tmp);
            free(output_graph);

            for (z = 0; z < nb_meta_prefix; z++){

                basename_meta_tmp = malloc((strlen(meta_prefix[z]) + 1) * sizeof(char));
                ASSERT_NULL_PTR(basename_meta_tmp, "main()")

                meta_filenames = malloc(nb_meta_stream * sizeof(char*));
                ASSERT_NULL_PTR(meta_filenames, "main()")

                strcpy(basename_meta_tmp, meta_prefix[z]);
                basename_meta = basename(basename_meta_tmp);
                tmp = strchr(basename_meta, '.');
                if (tmp != NULL) *tmp = '\0';
                len_tmp_meta = strlen(tmp_dir) + strlen(basename_meta);

                output_meta = malloc((len_tmp_meta + 50) * sizeof(char));
                ASSERT_NULL_PTR(output_meta, "main()")

                strcpy(output_meta, tmp_dir);
                strcpy(&output_meta[strlen(tmp_dir)], basename_meta);

                for (int l = 0; l < nb_meta_stream; l++){

                    meta_filenames[l] = malloc((len_tmp_meta + 50) * sizeof(char));
                    ASSERT_NULL_PTR(meta_filenames[l], "main()")

                    strcpy(meta_filenames[l], output_meta);
                }

                strcpy(&(meta_filenames[0][len_tmp_meta]), "_info");
                strcpy(&(meta_filenames[1][len_tmp_meta]), "_runs");
                strcpy(&(meta_filenames[2][len_tmp_meta]), "_partitions");
                strcpy(&(meta_filenames[3][len_tmp_meta]), "_read_info");
                strcpy(&(meta_filenames[4][len_tmp_meta]), "_recycl_parts1");
                strcpy(&(meta_filenames[5][len_tmp_meta]), "_recycl_parts2");
                strcpy(&(meta_filenames[6][len_tmp_meta]), "_shifts");
                strcpy(&(meta_filenames[7][len_tmp_meta]), "_span_sup_reads_pos");
                strcpy(&(meta_filenames[8][len_tmp_meta]), "_span_sup_reads_length");
                strcpy(&(meta_filenames[9][len_tmp_meta]), "_span_sup_reads_rev");
                strcpy(&(meta_filenames[10][len_tmp_meta]), "_span_sup_reads_occ");
                strcpy(&(meta_filenames[11][len_tmp_meta]), "_span_sup_reads_pos_mis");
                strcpy(&(meta_filenames[12][len_tmp_meta]), "_span_sup_reads_char_mis");
                strcpy(&(meta_filenames[13][len_tmp_meta]), "_span_sup_reads_id_mis");
                strcpy(&(meta_filenames[14][len_tmp_meta]), "_span_sup_reads_nb_mis");
                strcpy(&(meta_filenames[15][len_tmp_meta]), "_span_sup_reads_mate");
                strcpy(&(meta_filenames[16][len_tmp_meta]), "_span_sup_reads_mate_info");

                decode_lzma(meta_prefix[z], 1, meta_filenames);

                pair_ended = read_meta_info_file(output_meta);
                if (remove(meta_filenames[0])) printf("Warning: Could not remove temporary file.\n");

                decode_lzma(meta_prefix[z], nb_meta_stream - (!pair_ended * 2), meta_filenames);

                decompress_FASTx_file(output_meta, root_no_iupac, root_iupac, pair_ended, length_kmer, length_overlap);

                for (int l = 0; l < nb_meta_stream - (!pair_ended * 2); l++){
                    if (remove(meta_filenames[l])) printf("Warning: Could not remove temporary file.\n");
                }

                for (int l = 0; l < nb_meta_stream; l++) free(meta_filenames[l]);

                free(meta_filenames);
                free(basename_meta_tmp);
                free(output_meta);

                gettimeofday(&tval_after, NULL);
                time_spent(&tval_last, &tval_after, &tval_result);

                printf("\nElapsed time: %ld.%06ld s\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
                printf("Peak of memory: %llu mb\n", ((unsigned long long int)getPeakRSS())/1024);

                tval_last = tval_after;
            }

            freeBFT_Root(root_iupac);
            freeBFT_Root(root_no_iupac);
        }

        for (int l = 0; l < nb_meta_prefix; l++){

            free(meta_prefix[l]);

            if (compress){
                free(paths_and_names_mate1[l]);
                free(paths_and_names_mate2[l]);
            }
        }

        free(meta_prefix);
        free(paths_and_names_mate1);
        free(paths_and_names_mate2);

        free(tmp_dir);
        free(graph_filename);
    }

    return EXIT_SUCCESS;
}
