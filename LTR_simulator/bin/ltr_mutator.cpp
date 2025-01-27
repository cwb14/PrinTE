#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <unordered_map>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <mutex>
#include <omp.h>

// Function to display usage information
void print_usage() {
    std::cout << "Usage:\n"
              << "  mutator -fasta <input_fasta> -rate <mutation_rate> -generations <num_generations> -mode <0|1> -threads <num_threads>\n"
              << "\nArguments:\n"
              << "  -fasta        Path to the input multi-FASTA file.\n"
              << "  -rate         Mutation rate per generation (e.g., 1.3e-8).\n"
              << "  -generations  Number of generations to simulate (e.g., 10000).\n"
              << "  -mode         Mutation mode:\n"
              << "                  0 - Fixed number of generations per sequence.\n"
              << "                  1 - Random number of generations per sequence (0 to max).\n"
              << "  -threads      Number of threads to use for parallel processing.\n";
}

// Function to parse command-line arguments
bool parse_arguments(int argc, char* argv[], std::string &fasta_path, double &mutation_rate, int &generations, int &mode, int &threads) {
    if (argc != 11) {
        print_usage();
        return false;
    }

    for (int i = 1; i < argc; i += 2) {
        std::string arg = argv[i];
        if (arg == "-fasta") {
            fasta_path = argv[i + 1];
        } else if (arg == "-rate") {
            try {
                mutation_rate = std::stod(argv[i + 1]);
            } catch (const std::invalid_argument&) {
                std::cerr << "Error: Invalid mutation rate provided.\n";
                return false;
            }
        } else if (arg == "-generations") {
            try {
                generations = std::stoi(argv[i + 1]);
                if (generations < 0) {
                    std::cerr << "Error: Number of generations must be non-negative.\n";
                    return false;
                }
            } catch (const std::invalid_argument&) {
                std::cerr << "Error: Invalid number of generations provided.\n";
                return false;
            }
        } else if (arg == "-mode") {
            try {
                mode = std::stoi(argv[i + 1]);
                if (mode != 0 && mode != 1) {
                    std::cerr << "Error: Mode must be either 0 or 1.\n";
                    return false;
                }
            } catch (const std::invalid_argument&) {
                std::cerr << "Error: Invalid mode provided.\n";
                return false;
            }
        } else if (arg == "-threads") {
            try {
                threads = std::stoi(argv[i + 1]);
                if (threads <= 0) {
                    std::cerr << "Error: Number of threads must be a positive integer.\n";
                    return false;
                }
            } catch (const std::invalid_argument&) {
                std::cerr << "Error: Invalid number of threads provided.\n";
                return false;
            }
        } else {
            std::cerr << "Error: Unknown argument " << arg << "\n";
            print_usage();
            return false;
        }
    }
    return true;
}

// Function to load sequences from a FASTA file
bool load_fasta(const std::string &fasta_path, std::vector<std::pair<std::string, std::string>> &sequences) {
    std::ifstream infile(fasta_path);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open file " << fasta_path << "\n";
        return false;
    }

    std::string line, header, sequence;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!header.empty()) {
                sequences.emplace_back(header, sequence);
                sequence.clear();
            }
            header = line.substr(1); // Remove '>'
        } else {
            sequence += line;
        }
    }
    // Add the last sequence
    if (!header.empty()) {
        sequences.emplace_back(header, sequence);
    }

    infile.close();
    return true;
}

// Function to write mutated sequences to a FASTA file
bool write_fasta(const std::string &output_path, const std::vector<std::pair<std::string, std::string>> &sequences) {
    std::ofstream outfile(output_path);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open output file " << output_path << "\n";
        return false;
    }

    for (const auto &seq_pair : sequences) {
        outfile << ">" << seq_pair.first << "\n";
        // Wrap sequences at 60 characters per line for readability
        const std::string &seq = seq_pair.second;
        for (size_t i = 0; i < seq.length(); i += 60) {
            outfile << seq.substr(i, 60) << "\n";
        }
    }

    outfile.close();
    return true;
}

// Function to create substitution dictionary
std::unordered_map<char, std::vector<char>> get_substitution_dict() {
    std::unordered_map<char, std::vector<char>> subs_dict;
    subs_dict['A'] = {'C', 'G', 'T'};
    subs_dict['C'] = {'A', 'G', 'T'};
    subs_dict['G'] = {'A', 'C', 'T'};
    subs_dict['T'] = {'A', 'C', 'G'};
    subs_dict['a'] = {'c', 'g', 't'};
    subs_dict['c'] = {'a', 'g', 't'};
    subs_dict['g'] = {'a', 'c', 't'};
    subs_dict['t'] = {'a', 'c', 'g'};
    return subs_dict;
}

// Mutex for synchronized printing
std::mutex print_mutex;

// Function to extract base name from filename (removes the last extension)
std::string get_base_filename(const std::string &filename) {
    size_t last_dot = filename.find_last_of(".");
    if (last_dot == std::string::npos) {
        return filename; // No extension found
    } else {
        return filename.substr(0, last_dot);
    }
}

int main(int argc, char* argv[]) {
    // Variables to store command-line arguments
    std::string fasta_path;
    double mutation_rate;
    int generations;
    int mode;
    int threads;

    // Parse command-line arguments
    if (!parse_arguments(argc, argv, fasta_path, mutation_rate, generations, mode, threads)) {
        return EXIT_FAILURE;
    }

    // Load sequences from FASTA file
    std::vector<std::pair<std::string, std::string>> sequences;
    std::cout << "Loading sequences from " << fasta_path << "...\n";
    auto start_time = std::chrono::high_resolution_clock::now();
    if (!load_fasta(fasta_path, sequences)) {
        return EXIT_FAILURE;
    }
    size_t total_sequences = sequences.size();
    std::cout << "Total sequences loaded: " << total_sequences << "\n";

    // Calculate Genome_size: sum of all sequence lengths
    unsigned long long genome_size = 0;
    for (const auto &seq_pair : sequences) {
        genome_size += seq_pair.second.length();
    }

    // Prepare substitution dictionary
    std::unordered_map<char, std::vector<char>> substitution_dict = get_substitution_dict();

    // Prepare vector to hold mutated sequences
    std::vector<std::pair<std::string, std::string>> mutated_sequences(total_sequences);

    // Initialize total number of genomic mutations
    unsigned long long total_mutations_sum = 0;

    // Set number of threads for OpenMP
    omp_set_num_threads(threads);
    std::cout << "Starting mutation simulation with " << threads << " threads...\n";

    // Parallel processing of sequences
    #pragma omp parallel reduction(+: total_mutations_sum)
    {
        // Each thread has its own RNG
        std::mt19937 rng;
        std::random_device rd;
        rng.seed(rd());

        // Distributions
        std::uniform_real_distribution<> dist_real(0.0, 1.0);
        std::uniform_int_distribution<> dist_sub(0, 2); // 3 substitution options

        // Iterate through each sequence
        #pragma omp for schedule(dynamic)
        for (size_t idx = 0; idx < total_sequences; ++idx) {
            const std::string &header = sequences[idx].first;
            std::string sequence = sequences[idx].second; // Make a copy for mutation
            size_t seq_length = sequence.length();

            // Determine number of generations based on mode
            int num_generations;
            if (mode == 0) {
                num_generations = generations;
            } else { // mode == 1
                std::uniform_int_distribution<> dist_num_gen(0, generations);
                num_generations = dist_num_gen(rng);
            }

            // Calculate total mutations: Poisson with λ = generations * mutation_rate * sequence_length
            double lambda = static_cast<double>(num_generations) * mutation_rate * static_cast<double>(seq_length);
            std::poisson_distribution<> poisson_dist(lambda);
            int total_mutations = poisson_dist(rng);

            // Accumulate total mutations
            total_mutations_sum += total_mutations;

            // Debug: Print sequence ID and number of mutations
            {
                std::lock_guard<std::mutex> guard(print_mutex);
                std::cout << "Sequence ID: " << header << " | Generations: " << num_generations 
                          << " | Total Mutations: " << total_mutations << "\n";
            }

            // Apply mutations
            if (total_mutations > 0) {
                std::uniform_int_distribution<size_t> dist_pos(0, seq_length - 1);
                for (int mut = 0; mut < total_mutations; ++mut) {
                    // Randomly select a position to mutate
                    size_t pos = dist_pos(rng);
                    char current_nuc = sequence[pos];

                    // Check if current nucleotide is mutable
                    auto it = substitution_dict.find(current_nuc);
                    if (it != substitution_dict.end()) {
                        const std::vector<char> &options = it->second;
                        int sub_index = dist_sub(rng);
                        char new_nuc = options[sub_index];
                        sequence[pos] = new_nuc;

                        // No per-mutation print statements to improve efficiency
                    }
                    // If nucleotide is not A, C, G, T (e.g., N), skip mutation
                    // Recurrent mutations are allowed, so positions can be selected multiple times
                }
            }

            // Store mutated sequence
            mutated_sequences[idx] = std::make_pair(header, sequence);

            // Progress indicator every 10,000 sequences
            if ((idx + 1) % 10000 == 0) {
                std::lock_guard<std::mutex> guard(print_mutex);
                std::cout << "Mutated " << (idx + 1) << " / " << total_sequences << " sequences...\n";
            }
        }
    }

    std::cout << "Mutation simulation completed.\n";

    // Define output FASTA file name: <input_fasta_base>_<generations>.fa
    std::string base_filename = get_base_filename(fasta_path);
    std::string output_fasta = base_filename + "_" + std::to_string(generations) + ".fa";

    // Write mutated sequences to output FASTA file
    std::cout << "Writing mutated sequences to " << output_fasta << "...\n";
    if (!write_fasta(output_fasta, mutated_sequences)) {
        return EXIT_FAILURE;
    }
    std::cout << "All mutated sequences have been written successfully.\n";

    // If mode is 0, write the supplementary txt file
    if (mode == 0) {
        std::string output_txt = base_filename + "_" + std::to_string(generations) + ".txt";
        std::ofstream txtfile(output_txt);
        if (!txtfile.is_open()) {
            std::cerr << "Error: Unable to open supplementary file " << output_txt << "\n";
            return EXIT_FAILURE;
        }

        // Calculate fractions
        double fraction_mutated_bases = static_cast<double>(total_mutations_sum) / static_cast<double>(genome_size);
        double expected_divergence = fraction_mutated_bases * 2.0;

        // Write to the supplementary file
        txtfile << "Fraction of mutated bases: " << fraction_mutated_bases << "\n";
        txtfile << "Expected divergence: " << expected_divergence << "\n";

        txtfile.close();
        std::cout << "Supplementary file " << output_txt << " has been written successfully.\n";
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - start_time;
    std::cout << "Total time taken: " << diff.count() << " seconds.\n";

    return EXIT_SUCCESS;
}
