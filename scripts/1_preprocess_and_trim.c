#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <zlib.h>
#include <glob.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>

// Configuration constants
#define UMI1_LEN 7
#define UMI2_LEN 7
#define MAX_LINE_LEN 1024
#define MAX_SEQ_LEN 512
#define MAX_PATH_LEN 512

// TRA/TRB structure patterns
#define PRE_UMI1_TRA "GACTCTGATGACGACGCACA"
#define LINKER_FWD_TRA "GTACACGCTGGATCCGACTTGTAGA"
#define FLANK_TRA_SEQ "TACTCTGCTGATACCGATGC"

#define PRE_UMI1_TRB "GCATCGGTATCAGCAGAGTA"
#define LINKER_REV_TRB "TCTACAAGTCGGATCCAGCGTGTAC"
#define FLANK_TRB_SEQ "TGTGCGTCGTCATCAGAGTC"

// Progress tracking
typedef struct {
    long processed_pairs;
    long tra_pairs;
    long trb_pairs;
    long read_limit;
    time_t start_time;
} progress_t;

// FASTQ record structure
typedef struct {
    char header[MAX_LINE_LEN];
    char sequence[MAX_SEQ_LEN];
    char plus[MAX_LINE_LEN];
    char quality[MAX_SEQ_LEN];
} fastq_record_t;

// Function prototypes
void show_usage(const char *program_name);
int find_fastq_pair(const char *directory, char *r1_file, char *r2_file, char *base_name);
void reverse_complement(const char *seq, char *rc_seq);
int read_fastq_record(gzFile fp, fastq_record_t *record);
char *find_pattern(const char *sequence, const char *pattern);
int extract_umi_and_trim(const char *sequence, const char *pre_umi, const char *linker, 
                        const char *flank, char *umi1, char *umi2, char *trimmed_seq, int *found_rc);
void update_progress(progress_t *prog, int force_update);
int create_directory(const char *path);

int main(int argc, char *argv[]) {
    char input_dir[MAX_PATH_LEN] = "";
    char output_dir[MAX_PATH_LEN] = "PairTCR_results/1_preprocess_and_trim_output";
    char output_prefix[256] = "";
    long read_limit = 100000;
    
    // Parse command line arguments
    int opt;
    static struct option long_options[] = {
        {"limit", required_argument, 0, 'n'},
        {"output_prefix", required_argument, 0, 'o'},
        {"outdir", required_argument, 0, 'd'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };
    
    while ((opt = getopt_long(argc, argv, "n:o:d:h", long_options, NULL)) != -1) {
        switch (opt) {
            case 'n':
                read_limit = atol(optarg);
                break;
            case 'o':
                strncpy(output_prefix, optarg, sizeof(output_prefix) - 1);
                break;
            case 'd':
                strncpy(output_dir, optarg, sizeof(output_dir) - 1);
                break;
            case 'h':
                show_usage(argv[0]);
                return 0;
            default:
                show_usage(argv[0]);
                return 1;
        }
    }
    
    if (optind >= argc) {
        fprintf(stderr, "Error: Input directory required\n");
        show_usage(argv[0]);
        return 1;
    }
    
    strncpy(input_dir, argv[optind], sizeof(input_dir) - 1);
    
    // Find FASTQ pair
    char r1_file[MAX_PATH_LEN], r2_file[MAX_PATH_LEN], base_name[256];
    if (find_fastq_pair(input_dir, r1_file, r2_file, base_name) != 0) {
        return 1;
    }
    
    // Set output prefix if not provided
    if (strlen(output_prefix) == 0) {
        strncpy(output_prefix, base_name, sizeof(output_prefix) - 1);
        printf("Using '%s' as output prefix.\n", base_name);
    }
    
    // Create output directory
    if (create_directory(output_dir) != 0) {
        fprintf(stderr, "Error creating output directory: %s\n", output_dir);
        return 1;
    }
    
    // Construct output file paths
    char tra_r1_out[MAX_PATH_LEN], tra_r2_out[MAX_PATH_LEN];
    char trb_r1_out[MAX_PATH_LEN], trb_r2_out[MAX_PATH_LEN];
    
    snprintf(tra_r1_out, sizeof(tra_r1_out), "%s/%s_TRA_1.fq.gz", output_dir, output_prefix);
    snprintf(tra_r2_out, sizeof(tra_r2_out), "%s/%s_TRA_2.fq.gz", output_dir, output_prefix);
    snprintf(trb_r1_out, sizeof(trb_r1_out), "%s/%s_TRB_1.fq.gz", output_dir, output_prefix);
    snprintf(trb_r2_out, sizeof(trb_r2_out), "%s/%s_TRB_2.fq.gz", output_dir, output_prefix);
    
    // Open input files
    gzFile r1_in = gzopen(r1_file, "r");
    gzFile r2_in = gzopen(r2_file, "r");
    if (!r1_in || !r2_in) {
        fprintf(stderr, "Error opening input files\n");
        return 1;
    }
    
    // Open output files
    gzFile tra_r1_out_fp = gzopen(tra_r1_out, "w");
    gzFile tra_r2_out_fp = gzopen(tra_r2_out, "w");
    gzFile trb_r1_out_fp = gzopen(trb_r1_out, "w");
    gzFile trb_r2_out_fp = gzopen(trb_r2_out, "w");
    
    if (!tra_r1_out_fp || !tra_r2_out_fp || !trb_r1_out_fp || !trb_r2_out_fp) {
        fprintf(stderr, "Error opening output files\n");
        return 1;
    }
    
    // Initialize progress tracking
    progress_t progress = {0, 0, 0, read_limit, time(NULL)};
    
    printf("Starting processing...\n");
    printf("Input R1: %s\n", r1_file);
    printf("Input R2: %s\n", r2_file);
    printf("Read limit: %ld\n", read_limit);
    printf("Output directory: %s\n", output_dir);
    printf("Processing...\n");
    
    // Main processing loop
    fastq_record_t r1_record, r2_record;
    char r1_rc[MAX_SEQ_LEN], r2_rc[MAX_SEQ_LEN];
    char umi1[UMI1_LEN + 1], umi2[UMI2_LEN + 1];
    char trimmed_seq[MAX_SEQ_LEN], trimmed_qual[MAX_SEQ_LEN];
    
    while (progress.processed_pairs < read_limit) {
        // Read records
        if (read_fastq_record(r1_in, &r1_record) != 0 || 
            read_fastq_record(r2_in, &r2_record) != 0) {
            break;
        }
        
        progress.processed_pairs++;
        update_progress(&progress, 0);
        
        int read_type = 0; // 0=none, 1=TRA, 2=TRB
        int found_rc = 0;
        
        // Check R1 for TRA pattern
        if (extract_umi_and_trim(r1_record.sequence, PRE_UMI1_TRA, LINKER_FWD_TRA, 
                                FLANK_TRA_SEQ, umi1, umi2, trimmed_seq, &found_rc)) {
            read_type = 1; // TRA
            
            // Create modified headers
            char r1_header_mod[MAX_LINE_LEN], r2_header_mod[MAX_LINE_LEN];
            snprintf(r1_header_mod, sizeof(r1_header_mod), "%s UMI:TRA:%s_%s%s", 
                    r1_record.header, umi1, umi2, found_rc ? ":RC" : "");
            snprintf(r2_header_mod, sizeof(r2_header_mod), "%s UMI:TRA:%s_%s%s", 
                    r2_record.header, umi1, umi2, found_rc ? ":RC" : "");
            
            // Calculate trimmed quality
            int trim_pos = strlen(r1_record.sequence) - strlen(trimmed_seq);
            if (found_rc) {
                trim_pos = strlen(trimmed_seq);
                strncpy(trimmed_qual, r1_record.quality, trim_pos);
                trimmed_qual[trim_pos] = '\0';
            } else {
                strcpy(trimmed_qual, r1_record.quality + trim_pos);
            }
            
            // Write TRA output
            if (strlen(trimmed_seq) > 0 && strlen(r2_record.sequence) > 0) {
                gzprintf(tra_r1_out_fp, "%s\n%s\n%s\n%s\n", 
                        r1_header_mod, trimmed_seq, r1_record.plus, trimmed_qual);
                gzprintf(tra_r2_out_fp, "%s\n%s\n%s\n%s\n", 
                        r2_header_mod, r2_record.sequence, r2_record.plus, r2_record.quality);
                progress.tra_pairs++;
            }
        }
        
        // Check R2 for TRB pattern (only if not TRA)
        if (read_type == 0) {
            found_rc = 0;
            if (extract_umi_and_trim(r2_record.sequence, PRE_UMI1_TRB, LINKER_REV_TRB, 
                                    FLANK_TRB_SEQ, umi1, umi2, trimmed_seq, &found_rc)) {
                read_type = 2; // TRB
                
                // Create modified headers
                char r1_header_mod[MAX_LINE_LEN], r2_header_mod[MAX_LINE_LEN];
                snprintf(r1_header_mod, sizeof(r1_header_mod), "%s UMI:TRB:%s_%s%s", 
                        r1_record.header, umi1, umi2, found_rc ? ":RC" : "");
                snprintf(r2_header_mod, sizeof(r2_header_mod), "%s UMI:TRB:%s_%s%s", 
                        r2_record.header, umi1, umi2, found_rc ? ":RC" : "");
                
                // Calculate trimmed quality
                int trim_pos = strlen(r2_record.sequence) - strlen(trimmed_seq);
                if (found_rc) {
                    trim_pos = strlen(trimmed_seq);
                    strncpy(trimmed_qual, r2_record.quality, trim_pos);
                    trimmed_qual[trim_pos] = '\0';
                } else {
                    strcpy(trimmed_qual, r2_record.quality + trim_pos);
                }
                
                // Write TRB output
                if (strlen(r1_record.sequence) > 0 && strlen(trimmed_seq) > 0) {
                    gzprintf(trb_r1_out_fp, "%s\n%s\n%s\n%s\n", 
                            r1_header_mod, r1_record.sequence, r1_record.plus, r1_record.quality);
                    gzprintf(trb_r2_out_fp, "%s\n%s\n%s\n%s\n", 
                            r2_header_mod, trimmed_seq, r2_record.plus, trimmed_qual);
                    progress.trb_pairs++;
                }
            }
        }
    }
    
    // Final progress update
    update_progress(&progress, 1);
    
    // Close files
    gzclose(r1_in);
    gzclose(r2_in);
    gzclose(tra_r1_out_fp);
    gzclose(tra_r2_out_fp);
    gzclose(trb_r1_out_fp);
    gzclose(trb_r2_out_fp);
    
    // Print summary
    printf("\n--- Processing Summary ---\n");
    printf("Processed %ld read pairs (limit was %ld).\n", progress.processed_pairs, read_limit);
    printf("TRA pairs identified (UMI added, R1 trimmed to downstream): %ld\n", progress.tra_pairs);
    printf("TRB pairs identified (UMI added, R2 trimmed to downstream): %ld\n", progress.trb_pairs);
    printf("Output files written to directory: %s\n", output_dir);
    
    return 0;
}

void show_usage(const char *program_name) {
    printf("Usage: %s [OPTIONS] input_dir\n", program_name);
    printf("Process paired-end FASTQ files to identify TRA/TRB structures\n\n");
    printf("Options:\n");
    printf("  -n, --limit LIMIT        Maximum number of read pairs to process (default: 100000)\n");
    printf("  -o, --output_prefix PREFIX  Prefix for output files\n");
    printf("  -d, --outdir DIR         Output directory (default: PairTCR_results/1_preprocess_and_trim_output)\n");
    printf("  -h, --help               Show this help message\n");
}

int find_fastq_pair(const char *directory, char *r1_file, char *r2_file, char *base_name) {
    glob_t glob_result;
    char pattern[MAX_PATH_LEN];
    
    snprintf(pattern, sizeof(pattern), "%s/*_1.fq.gz", directory);
    
    if (glob(pattern, 0, NULL, &glob_result) != 0 || glob_result.gl_pathc == 0) {
        fprintf(stderr, "Error: No R1 file (*_1.fq.gz) found in directory: %s\n", directory);
        globfree(&glob_result);
        return 1;
    }
    
    strcpy(r1_file, glob_result.gl_pathv[0]);
    globfree(&glob_result);
    
    // Extract base name
    char *filename = strrchr(r1_file, '/');
    filename = filename ? filename + 1 : r1_file;
    strcpy(base_name, filename);
    
    char *ext_pos = strstr(base_name, "_1.fq.gz");
    if (ext_pos) {
        *ext_pos = '\0';
    }
    
    // Construct R2 filename
    snprintf(r2_file, MAX_PATH_LEN, "%s/%s_2.fq.gz", directory, base_name);
    
    // Check if R2 file exists
    if (access(r2_file, F_OK) != 0) {
        fprintf(stderr, "Error: Corresponding R2 file not found: %s\n", r2_file);
        return 1;
    }
    
    printf("Found FASTQ pair: %s_1.fq.gz, %s_2.fq.gz\n", base_name, base_name);
    return 0;
}

void reverse_complement(const char *seq, char *rc_seq) {
    int len = strlen(seq);
    for (int i = 0; i < len; i++) {
        switch (seq[len - 1 - i]) {
            case 'A': rc_seq[i] = 'T'; break;
            case 'T': rc_seq[i] = 'A'; break;
            case 'C': rc_seq[i] = 'G'; break;
            case 'G': rc_seq[i] = 'C'; break;
            case 'N': rc_seq[i] = 'N'; break;
            default: rc_seq[i] = 'N'; break;
        }
    }
    rc_seq[len] = '\0';
}

int read_fastq_record(gzFile fp, fastq_record_t *record) {
    if (gzgets(fp, record->header, sizeof(record->header)) == NULL) return 1;
    if (gzgets(fp, record->sequence, sizeof(record->sequence)) == NULL) return 1;
    if (gzgets(fp, record->plus, sizeof(record->plus)) == NULL) return 1;
    if (gzgets(fp, record->quality, sizeof(record->quality)) == NULL) return 1;
    
    // Remove newlines
    record->header[strcspn(record->header, "\n")] = '\0';
    record->sequence[strcspn(record->sequence, "\n")] = '\0';
    record->plus[strcspn(record->plus, "\n")] = '\0';
    record->quality[strcspn(record->quality, "\n")] = '\0';
    
    return 0;
}

char *find_pattern(const char *sequence, const char *pattern) {
    return strstr(sequence, pattern);
}

int extract_umi_and_trim(const char *sequence, const char *pre_umi, const char *linker, 
                        const char *flank, char *umi1, char *umi2, char *trimmed_seq, int *found_rc) {
    char rc_sequence[MAX_SEQ_LEN];
    const char *search_seq = sequence;
    *found_rc = 0;
    
    // Build pattern to search for
    char pattern[256];
    snprintf(pattern, sizeof(pattern), "%s", pre_umi);
    
    char *match_pos = find_pattern(search_seq, pattern);
    
    if (!match_pos) {
        // Try reverse complement
        reverse_complement(sequence, rc_sequence);
        search_seq = rc_sequence;
        match_pos = find_pattern(search_seq, pattern);
        if (match_pos) {
            *found_rc = 1;
        }
    }
    
    if (!match_pos) {
        return 0; // Pattern not found
    }
    
    // Extract UMIs and find complete pattern
    char *pos = match_pos + strlen(pre_umi);
    
    // Extract UMI1
    strncpy(umi1, pos, UMI1_LEN);
    umi1[UMI1_LEN] = '\0';
    pos += UMI1_LEN;
    
    // Check linker
    if (strncmp(pos, linker, strlen(linker)) != 0) {
        return 0;
    }
    pos += strlen(linker);
    
    // Extract UMI2
    strncpy(umi2, pos, UMI2_LEN);
    umi2[UMI2_LEN] = '\0';
    pos += UMI2_LEN;
    
    // Check flank
    if (strncmp(pos, flank, strlen(flank)) != 0) {
        return 0;
    }
    pos += strlen(flank);
    
    // Check for A or T
    if (*pos != 'A' && *pos != 'T') {
        return 0;
    }
    pos++;
    
    // Extract trimmed sequence
    if (*found_rc) {
        // For RC, take sequence before the pattern start in original sequence
        int pattern_start = strlen(sequence) - (pos - search_seq);
        strncpy(trimmed_seq, sequence, pattern_start);
        trimmed_seq[pattern_start] = '\0';
    } else {
        // For forward, take sequence after the pattern end
        int pos_in_orig = pos - search_seq;
        strcpy(trimmed_seq, sequence + pos_in_orig);
    }
    
    return 1; // Success
}

void update_progress(progress_t *prog, int force_update) {
    static time_t last_update = 0;
    time_t now = time(NULL);
    
    if (force_update || (now - last_update) >= 1) {
        double elapsed = difftime(now, prog->start_time);
        double rate = prog->processed_pairs / (elapsed > 0 ? elapsed : 1);
        double eta = (prog->read_limit - prog->processed_pairs) / (rate > 0 ? rate : 1);
        
        printf("\rProcessing: %ld/%ld pairs (%.1f%%) | TRA: %ld | TRB: %ld | Rate: %.0f pairs/s | ETA: %.0fs   ", 
               prog->processed_pairs, prog->read_limit, 
               (100.0 * prog->processed_pairs) / prog->read_limit,
               prog->tra_pairs, prog->trb_pairs, rate, eta);
        fflush(stdout);
        last_update = now;
    }
}

int create_directory(const char *path) {
    char temp_path[MAX_PATH_LEN];
    char *p = NULL;
    size_t len;
    
    strncpy(temp_path, path, sizeof(temp_path) - 1);
    len = strlen(temp_path);
    
    if (temp_path[len - 1] == '/') {
        temp_path[len - 1] = '\0';
    }
    
    for (p = temp_path + 1; *p; p++) {
        if (*p == '/') {
            *p = '\0';
            if (mkdir(temp_path, 0755) != 0 && errno != EEXIST) {
                return -1;
            }
            *p = '/';
        }
    }
    
    if (mkdir(temp_path, 0755) != 0 && errno != EEXIST) {
        return -1;
    }
    
    return 0;
} 