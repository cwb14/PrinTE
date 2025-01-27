#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <deque>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <mutex>
#include <thread>
#include <random>
#include <algorithm>
#include <functional>
#include <atomic>

using namespace std;

// Function to parse the genome FASTA file
map<string, deque<char>> parse_fasta(const string& file_path) {
    map<string, deque<char>> sequences;
    ifstream infile(file_path);
    if(!infile) {
        cerr << "Error opening file: " << file_path << endl;
        exit(1);
    }
    string line;
    string current_header;
    while(getline(infile, line)) {
        if(line.empty()) continue;
        if(line[0] == '>') {
            current_header = line.substr(1);
            size_t pos = current_header.find(' ');
            if(pos != string::npos) {
                current_header = current_header.substr(0, pos);
            }
            sequences[current_header] = deque<char>();
        } else {
            // Append sequence
            if(current_header.empty()) {
                cerr << "Error: Sequence without header in fasta file." << endl;
                exit(1);
            }
            for(char c : line) {
                sequences[current_header].push_back(toupper(c));
            }
        }
    }
    infile.close();
    cout << "Completed reading " << sequences.size() << " sequences from " << file_path << endl;
    return sequences;
}

// Function to parse the LTR FASTA file and categorize sequences into 'solo' and 'intact'
void parse_ltr_fasta(const string& file_path, vector<string>& solo_sequences, vector<string>& intact_sequences) {
    ifstream infile(file_path);
    if(!infile) {
        cerr << "Error opening file: " << file_path << endl;
        exit(1);
    }
    cout << "Reading LTR fasta file: " << file_path << endl;
    string line;
    string current_seq;
    string current_type;
    while(getline(infile, line)) {
        if(line.empty()) continue;
        if(line[0] == '>') {
            // Process previous sequence
            if(!current_seq.empty() && !current_type.empty()) {
                if(current_type == "solo") {
                    solo_sequences.push_back(current_seq);
                } else if(current_type == "intact") {
                    intact_sequences.push_back(current_seq);
                }
                current_seq.clear();
            }
            // Determine type from header
            string header = line.substr(1);
            if(header.find("_solo") != string::npos) {
                current_type = "solo";
            } else if(header.find("_intact") != string::npos) {
                current_type = "intact";
            } else {
                cerr << "Warning: Unknown LTR type in header '" << header << "'. Skipping this sequence." << endl;
                current_type.clear();
            }
        } else {
            if(!current_type.empty()) {
                for(char c : line) {
                    current_seq += toupper(c);
                }
            }
        }
    }
    // Add the last sequence
    if(!current_seq.empty() && !current_type.empty()) {
        if(current_type == "solo") {
            solo_sequences.push_back(current_seq);
        } else if(current_type == "intact") {
            intact_sequences.push_back(current_seq);
        }
    }
    infile.close();
    cout << "Completed reading " << solo_sequences.size() << " solo and " << intact_sequences.size() << " intact LTR sequences from " << file_path << endl;
}

// Function to get cumulative chromosome lengths for weighted random selection
void get_cumulative_lengths(const map<string, deque<char>>& genome_sequences, vector<string>& chromosomes, vector<size_t>& cumulative_lengths) {
    size_t cumulative = 0;
    for(const auto& kv : genome_sequences) {
        chromosomes.push_back(kv.first);
        cumulative += kv.second.size();
        cumulative_lengths.push_back(cumulative);
    }
}

// Function to write the modified genome to an output FASTA file
void write_genome(const string& output_file, const map<string, deque<char>>& genome_sequences) {
    ofstream outfile(output_file);
    if(!outfile) {
        cerr << "Error opening output file: " << output_file << endl;
        exit(1);
    }
    for(const auto& kv : genome_sequences) {
        outfile << ">" << kv.first << endl;
        const deque<char>& seq = kv.second;
        size_t seq_len = seq.size();
        size_t i = 0;
        while(i < seq_len) {
            size_t line_len = min(static_cast<size_t>(60), seq_len - i);
            for(size_t j = 0; j < line_len; ++j) {
                outfile << seq[i + j];
            }
            outfile << endl;
            i += line_len;
        }
    }
    outfile.close();
}

int main(int argc, char* argv[]) {
    // Declare variables for command-line arguments
    string genome_file;
    string ltr_file;
    double rate = 0.0;
    int generations = 0;
    string solo_intact_ratio_str = "1:1";
    string output_file = "genome_modified.fa";
    int num_threads = 1;

    // Parse command-line arguments
    for(int i = 1; i < argc; ++i) {
        if(strcmp(argv[i], "--genome") == 0 && i+1 < argc) {
            genome_file = argv[++i];
        } else if(strcmp(argv[i], "--ltr") == 0 && i+1 < argc) {
            ltr_file = argv[++i];
        } else if(strcmp(argv[i], "--rate") == 0 && i+1 < argc) {
            rate = atof(argv[++i]);
        } else if(strcmp(argv[i], "--generations") == 0 && i+1 < argc) {
            generations = atoi(argv[++i]);
        } else if(strcmp(argv[i], "--solo_intact_ratio") == 0 && i+1 < argc) {
            solo_intact_ratio_str = argv[++i];
        } else if(strcmp(argv[i], "--output") == 0 && i+1 < argc) {
            output_file = argv[++i];
        } else if(strcmp(argv[i], "--threads") == 0 && i+1 < argc) {
            num_threads = atoi(argv[++i]);
        } else {
            cerr << "Unknown argument: " << argv[i] << endl;
            return 1;
        }
    }

    // Check required arguments
    if(genome_file.empty() || ltr_file.empty() || rate == 0.0 || generations == 0) {
        cerr << "Usage: " << argv[0] << " --genome <file> --ltr <file> --rate <float> --generations <int> [--solo_intact_ratio <string>] [--output <file>] [--threads <int>]" << endl;
        return 1;
    }

    // Parse solo_intact_ratio
    double solo_ratio = 1.0;
    double intact_ratio = 1.0;
    size_t pos = solo_intact_ratio_str.find(':');
    if(pos != string::npos) {
        string solo_str = solo_intact_ratio_str.substr(0, pos);
        string intact_str = solo_intact_ratio_str.substr(pos+1);
        solo_ratio = atof(solo_str.c_str());
        intact_ratio = atof(intact_str.c_str());
    } else {
        cerr << "Error parsing --solo_intact_ratio '" << solo_intact_ratio_str << "'" << endl;
        return 1;
    }

    if(solo_ratio < 0 || intact_ratio < 0) {
        cerr << "Ratios must be non-negative." << endl;
        return 1;
    }

    if(solo_ratio == 0 && intact_ratio == 0) {
        cerr << "At least one ratio must be greater than zero." << endl;
        return 1;
    }

    cout << "Parsed solo_intact_ratio: solo=" << solo_ratio << ", intact=" << intact_ratio << endl;

    // Read genome sequences
    cout << "Reading genome file: " << genome_file << endl;
    map<string, deque<char>> genome_sequences = parse_fasta(genome_file);

    // Read LTR sequences
    vector<string> solo_ltr_list;
    vector<string> intact_ltr_list;
    parse_ltr_fasta(ltr_file, solo_ltr_list, intact_ltr_list);

    // Check that we have sequences
    if(solo_ltr_list.empty() && solo_ratio > 0) {
        cerr << "Error: No 'solo' LTR sequences found, but solo_ratio > 0." << endl;
        exit(1);
    }
    if(intact_ltr_list.empty() && intact_ratio > 0) {
        cerr << "Error: No 'intact' LTR sequences found, but intact_ratio > 0." << endl;
        exit(1);
    }

    // Calculate total genome length
    size_t genome_length = 0;
    for(const auto& kv : genome_sequences) {
        genome_length += kv.second.size();
    }
    cout << "Total genome length: " << genome_length << " bases." << endl;

    // Calculate total insertions
    int total_insertions = static_cast<int>(rate * genome_length * generations);
    cout << "Total insertions to simulate: " << total_insertions << endl;

    // Determine the prefix for the txt file by removing .fa or .fasta from output_file
    string output_prefix = output_file;
    size_t fa_pos = output_prefix.rfind(".fa");
    size_t fasta_pos = output_prefix.rfind(".fasta");
    if(fa_pos != string::npos && fa_pos == output_prefix.length() - 3) {
        output_prefix = output_prefix.substr(0, fa_pos);
    }
    else if(fasta_pos != string::npos && fasta_pos == output_prefix.length() - 6) {
        output_prefix = output_prefix.substr(0, fasta_pos);
    }
    // Else, keep the original name without modification

    // Create the txt filename
    string txt_filename = output_prefix + ".txt";

    // Write the required message to the txt file
    ofstream txt_file(txt_filename);
    if(!txt_file) {
        cerr << "Error creating TXT file: " << txt_filename << endl;
        exit(1);
    }
    txt_file << "Simulating " << generations << " generations. Inserting " << total_insertions << " unique LTRs." << endl;
    txt_file.close();
    cout << "Written simulation details to " << txt_filename << endl;

    if(total_insertions == 0) {
        cerr << "Warning: Total insertions calculated as 0. Please check rate and generations." << endl;
        exit(1);
    }

    // Prepare for insertion site selection
    vector<string> chromosomes;
    vector<size_t> cumulative_lengths;
    get_cumulative_lengths(genome_sequences, chromosomes, cumulative_lengths);

    // Create per-chromosome mutexes
    map<string, mutex> chromosome_mutexes;
    for(const auto& chrom : chromosomes) {
        chromosome_mutexes[chrom];  // Default constructor creates a mutex
    }

    // Calculate insertion probabilities
    double total_weight = solo_ratio + intact_ratio;
    double prob_solo = solo_ratio / total_weight;
    double prob_intact = intact_ratio / total_weight;
    cout << "Insertion probabilities: solo=" << prob_solo << ", intact=" << prob_intact << endl;

    // For random number generation
    random_device rd;
    mt19937 gen(rd());

    // Shared flag for first insertion verification
    atomic<bool> first_insertion_done(false);

    // Create threads
    vector<thread> threads;

    int insertions_per_thread = total_insertions / num_threads;
    int remainder = total_insertions % num_threads;

    cout << "Starting insertions with " << num_threads << " threads..." << endl;

    for(int t = 0; t < num_threads; ++t) {
        int num_insertions = insertions_per_thread;
        if(t < remainder) {
            num_insertions += 1;
        }
        threads.push_back(thread([&, num_insertions, t]() {
            // Thread function
            // Initialize thread-local random number generators
            mt19937 thread_gen(rd());
            uniform_int_distribution<size_t> genome_dist(1, cumulative_lengths.back());

            // For solo/intact selection
            uniform_real_distribution<double> prob_dist(0.0, 1.0);

            // For progress reporting
            int report_interval = max(1, num_insertions / 10);

            for(int i = 0; i < num_insertions; ++i) {
                // Randomly select a position in the genome
                size_t rand_num = genome_dist(thread_gen);

                // Binary search to find the chromosome
                auto it = lower_bound(cumulative_lengths.begin(), cumulative_lengths.end(), rand_num);
                int idx = distance(cumulative_lengths.begin(), it);
                string chrom = chromosomes[idx];

                // Position within chromosome
                size_t chrom_start = (idx == 0) ? 0 : cumulative_lengths[idx - 1];
                size_t pos_in_chrom = rand_num - chrom_start - 1;

                // Decide whether to insert solo or intact
                string insertion_type;
                double rand_prob = prob_dist(thread_gen);
                if(solo_ratio > 0 && intact_ratio > 0) {
                    if(rand_prob < prob_solo) {
                        insertion_type = "solo";
                    } else {
                        insertion_type = "intact";
                    }
                } else if(solo_ratio > 0) {
                    insertion_type = "solo";
                } else if(intact_ratio > 0) {
                    insertion_type = "intact";
                } else {
                    cerr << "Error: Unable to determine insertion type." << endl;
                    continue;
                }

                // Select LTR sequence
                string ltr_seq;
                if(insertion_type == "solo") {
                    uniform_int_distribution<size_t> solo_dist(0, solo_ltr_list.size() - 1);
                    ltr_seq = solo_ltr_list[solo_dist(thread_gen)];
                } else if(insertion_type == "intact") {
                    uniform_int_distribution<size_t> intact_dist(0, intact_ltr_list.size() - 1);
                    ltr_seq = intact_ltr_list[intact_dist(thread_gen)];
                }

                // Acquire mutex for chromosome
                {
                    lock_guard<mutex> lock(chromosome_mutexes[chrom]);

                    deque<char>& seq = genome_sequences[chrom];

                    if(pos_in_chrom < 5 || pos_in_chrom >= seq.size()) {
                        continue;  // Invalid position
                    }

                    // Extract motif
                    string motif;
                    for(size_t k = pos_in_chrom - 5; k < pos_in_chrom; ++k) {
                        motif += seq[k];
                    }

                    // Construct insertion
                    string insertion = ltr_seq + motif;

                    // If first insertion, perform TSD verification
                    if(!first_insertion_done.load()) {
                        size_t upstream_start = (pos_in_chrom >= 15) ? pos_in_chrom - 15 : 0;
                        size_t downstream_end = min(pos_in_chrom + 10, seq.size());

                        string before_insertion(seq.begin() + upstream_start, seq.begin() + downstream_end);

                        // Insert temporarily
                        deque<char> temp_seq = seq;
                        temp_seq.insert(temp_seq.begin() + pos_in_chrom, insertion.begin(), insertion.end());

                        downstream_end += insertion.size();
                        string after_insertion(temp_seq.begin() + upstream_start, temp_seq.begin() + downstream_end);

                        cout << "\n=== TSD Functionality Verification ===" << endl;
                        cout << "Chromosome: " << chrom << endl;
                        cout << "Insertion Position: " << pos_in_chrom << endl;
                        cout << "\n-- Before Insertion --" << endl;
                        cout << before_insertion << endl;
                        cout << "\n-- After Insertion --" << endl;
                        cout << after_insertion << endl;
                        cout << "=== End of Verification ===\n" << endl;

                        first_insertion_done.store(true);
                    }

                    // Insert into sequence
                    seq.insert(seq.begin() + pos_in_chrom, insertion.begin(), insertion.end());
                }

                // Report progress
                if((i + 1) % report_interval == 0) {
                    cout << "Thread " << t << ": Inserted " << (i + 1) << " / " << num_insertions << " insertions." << endl;
                }
            }
        }));
    }

    // Wait for threads to complete
    for(auto& th : threads) {
        th.join();
    }

    cout << "All insertions completed." << endl;

    // Write modified genome to output file
    cout << "Writing modified genome to " << output_file << "..." << endl;
    write_genome(output_file, genome_sequences);
    cout << "Finished writing the modified genome." << endl;

    return 0;
}
