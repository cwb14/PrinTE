#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <unordered_map>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cmath>
#include <omp.h>
#include <sstream>
#include <algorithm>
#include <mutex>
#include <stdexcept>

// --- Data structures for BED and segments ---
struct BedRegion {
    size_t start, end;
    std::string attr, extra;
};

struct Segment {
    size_t start, end;
    int generations;  // if negative, means random up to G
    Segment(size_t s, size_t e, int g) : start(s), end(e), generations(g) {}
};

// --- Usage message ---
void print_usage() {
    std::cerr << "Usage: mutator -fasta <in.fa> -rate <u> -generations <G> -mode <0|1|2|3>"
                 " -threads <T> -seed <S> -out_prefix <P> [-bed <regions.bed>]\n";
}

// --- Parse args ---
bool parse_args(int argc, char* argv[],
                std::string &fasta, double &mu,
                int &G, int &mode, int &threads,
                int &seed, std::string &prefix,
                std::string &bedfile)
{
    if (argc != 15 && argc != 17) {
        print_usage();
        return false;
    }
    for (int i = 1; i + 1 < argc; i += 2) {
        std::string a = argv[i], v = argv[i+1];
        if      (a == "-fasta")       fasta    = v;
        else if (a == "-rate")        mu       = std::stod(v);
        else if (a == "-generations") G        = std::stoi(v);
        else if (a == "-mode")        mode     = std::stoi(v);
        else if (a == "-threads")     threads  = std::stoi(v);
        else if (a == "-seed")        seed     = std::stoi(v);
        else if (a == "-out_prefix")  prefix   = v;
        else if (a == "-bed")         bedfile  = v;
        else {
            std::cerr << "Unknown argument: " << a << "\n";
            print_usage();
            return false;
        }
    }
    if ((mode == 2 || mode == 3) && bedfile.empty()) {
        std::cerr << "Error: -bed required for mode " << mode << "\n";
        return false;
    }
    return true;
}

// --- Load BED into map: chrom -> vector<BedRegion> ---
bool load_bed(const std::string &bedfile,
              std::unordered_map<std::string,
                std::vector<BedRegion>> &bed_map)
{
    std::ifstream in(bedfile);
    if (!in) {
        std::cerr << "Cannot open BED file: " << bedfile << "\n";
        return false;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        BedRegion r;
        std::string na, strand;
        if (!(ss >> r.attr >> r.start >> r.end >> r.attr >> na >> strand)) {
            std::cerr << "Skipping malformed BED line: " << line << "\n";
            continue;
        }
        ss >> r.extra;  // optional
        // only keep valid intervals
        if (r.end <= r.start) continue;
        bed_map[r.attr].push_back(r);
    }
    for (auto &kv : bed_map) {
        auto &v = kv.second;
        std::sort(v.begin(), v.end(), [](auto &a, auto &b){
            return a.start < b.start;
        });
    }
    return true;
}

// --- Load FASTA ---
bool load_fasta(const std::string &fasta,
                std::vector<std::pair<std::string,std::string>> &seqs)
{
    std::ifstream in(fasta);
    if (!in) {
        std::cerr << "Cannot open FASTA file: " << fasta << "\n";
        return false;
    }
    std::string line, hdr, seq;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!hdr.empty() && !seq.empty()) {
                seqs.emplace_back(hdr, seq);
            }
            hdr = line.substr(1);
            seq.clear();
        } else {
            seq += line;
        }
    }
    if (!hdr.empty() && !seq.empty()) {
        seqs.emplace_back(hdr, seq);
    }
    return true;
}

// --- Write FASTA ---
bool write_fasta(const std::string &out,
                 const std::vector<std::pair<std::string,std::string>> &seqs)
{
    std::ofstream o(out);
    if (!o) {
        std::cerr << "Cannot write FASTA file: " << out << "\n";
        return false;
    }
    for (const auto &p : seqs) {
        o << ">" << p.first << "\n";
        const auto &s = p.second;
        for (size_t i = 0; i < s.size(); i += 60) {
            o << s.substr(i, 60) << "\n";
        }
    }
    return true;
}

// --- Substitution map ---
auto get_subs() {
    static const std::unordered_map<char,std::vector<char>> m = {
        {'A',{'C','G','T'}}, {'C',{'A','G','T'}}, {'G',{'A','C','T'}}, {'T',{'A','C','G'}},
        {'a',{'c','g','t'}}, {'c',{'a','g','t'}}, {'g',{'a','c','t'}}, {'t',{'a','c','g'}}
    };
    return m;
}

int main(int argc, char* argv[]) {
    std::string fasta, prefix, bedfile;
    double mu = 0.0;
    int G = 0, mode = 0, threads = 1, seed = 0;

    if (!parse_args(argc, argv, fasta, mu, G, mode, threads, seed, prefix, bedfile)) {
        return EXIT_FAILURE;
    }

    if (threads < 1) threads = 1;
    int max_threads = omp_get_max_threads();
    if (threads > max_threads) {
        std::cerr << "Requested " << threads << " threads, but only "
                  << max_threads << " available. Using " << max_threads << ".\n";
        threads = max_threads;
    }
    omp_set_dynamic(0);
    omp_set_num_threads(threads);

    std::unordered_map<std::string, std::vector<BedRegion>> bed_map;
    if (mode == 2 || mode == 3) {
        std::cout << "Loading BED " << bedfile << "...\n";
        if (!load_bed(bedfile, bed_map)) return EXIT_FAILURE;
        std::cout << "  " << bed_map.size() << " chromosomes indexed.\n";
    }

    std::cout << "Loading FASTA " << fasta << "...\n";
    auto t0 = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<std::string,std::string>> seqs;
    if (!load_fasta(fasta, seqs)) return EXIT_FAILURE;
    size_t N = seqs.size();
    if (N == 0) {
        std::cerr << "No sequences loaded.\n";
        return EXIT_FAILURE;
    }
    std::cout << "  " << N << " sequences.\n";

    unsigned long long genome_size = 0;
    for (const auto &p : seqs) genome_size += p.second.size();

    auto subs = get_subs();
    std::vector<std::pair<std::string,std::string>> out_seqs(N);
    unsigned long long total_muts = 0;

    #pragma omp parallel for schedule(static) reduction(+:total_muts)
    for (size_t i = 0; i < N; ++i) {
        std::mt19937 rng(seed + omp_get_thread_num());
        std::uniform_int_distribution<> randG(0, G);
        std::uniform_int_distribution<> pickSub(0, 2);

        const auto &hdr = seqs[i].first;
        auto seq = seqs[i].second;
        size_t L = seq.size();

        // Build segments
        std::vector<Segment> segs;
        if (mode < 2) {
            int g = (mode == 0 ? G : randG(rng));
            if (L > 0)
                segs.emplace_back(0, L, g);
        } else {
            size_t prev = 0;
            auto it = bed_map.find(hdr);
            if (it != bed_map.end()) {
                for (const auto &r : it->second) {
                    if (r.start > prev && prev < L)
                        segs.emplace_back(prev, std::min(r.start, L), G);
                    int g = G;
                    if (mode == 2) {
                        g = (r.attr.rfind("gene", 0) == 0 ? G : randG(rng));
                    } else {
                        g = (r.extra == "new_TE" ? randG(rng) : G);
                    }
                    if (r.start < L && r.end > r.start)
                        segs.emplace_back(r.start, std::min(r.end, L), g);
                    prev = std::min(r.end, L);
                }
            }
            if (prev < L)
                segs.emplace_back(prev, L, G);
        }

        // Mutate each segment
        unsigned long long seq_muts = 0;
        for (const auto &S : segs) {
            if (S.end <= S.start) continue;
            double lambda = double(S.generations) * mu * double(S.end - S.start);
            std::poisson_distribution<> pd(lambda);
            int mcount = pd(rng);
            seq_muts += mcount;
            total_muts += mcount;
            if (mcount > 0) {
                std::uniform_int_distribution<size_t> pickPos(S.start, S.end - 1);
                for (int m = 0; m < mcount; ++m) {
                    size_t pos = pickPos(rng);
                    if (pos >= L) continue;
                    auto it2 = subs.find(seq[pos]);
                    if (it2 != subs.end())
                        seq[pos] = it2->second[pickSub(rng)];
                }
            }
        }

        out_seqs[i] = { hdr, std::move(seq) };
        if ((i+1) % 10000 == 0 && omp_get_thread_num() == 0)
            std::cout << "  processed " << (i+1) << "/" << N << " sequences\n";
    }

    // Write FASTA
    std::string outfa = prefix + ".fa";
    std::cout << "Writing " << outfa << "...\n";
    if (!write_fasta(outfa, out_seqs)) return EXIT_FAILURE;

    // Write metrics
    std::string outtxt = prefix + ".txt";
    std::ofstream ot(outtxt);
    if (!ot) {
        std::cerr << "Cannot write metrics file: " << outtxt << "\n";
        return EXIT_FAILURE;
    }
    double lambda_obs = double(total_muts) / double(genome_size);
    double U = genome_size * (1.0 - std::exp(-lambda_obs));
    double rec = double(total_muts) - U;
    double mu1 = U / double(genome_size);
    double mu2 = mu1 * 2.0;
    ot << "Genome size: " << genome_size << "\n"
       << "Total muts: " << total_muts << "\n"
       << "Recurrent est: " << rec << "\n"
       << "Muts/site: " << mu1 << "\n"
       << "Muts/site*2: " << mu2 << "\n";

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt = t1 - t0;
    std::cout << "Done in " << dt.count() << "s\n";
    return EXIT_SUCCESS;
}
