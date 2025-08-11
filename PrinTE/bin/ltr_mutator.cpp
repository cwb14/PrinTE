/**********************************************************************
 * mutator.cpp  —  point-mutation simulator with BED-aware segments
 *                random generations now follow an exponential decay
 *                distribution (λ = 1.4 × 10⁻⁶) instead of uniform.
 *
 *   linux: g++ -std=c++17 -fopenmp -O3 -o mutator mutator.cpp
 *
 *
 *   arm64: # Step 1: Compile only (no linking)
 *   $(brew --prefix llvm)/bin/clang++ -std=c++17 -O3 -fopenmp \
 *   -I$(brew --prefix libomp)/include \
 *   -c TESS/PrinTE/bin/ltr_mutator.cpp -o mutator.o
 *
 *   # Step 2: Link manually with static libomp
 *   $(brew --prefix llvm)/bin/clang++ -O3 mutator.o \
 *   $(brew --prefix libomp)/lib/libomp.a \
 *   -lc++ -lm -o TESS/PrinTE/bin/ltr_mutator_mac
 *
 *********************************************************************/

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <random>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// ────────────────────────────────────────────────────────────────────
// Data structures
// ────────────────────────────────────────────────────────────────────
struct BedRegion {
    size_t start, end;
    std::string chrom;
    std::string extra;
};

struct Segment {
    size_t start, end;
    int generations;
    Segment(size_t s, size_t e, int g) : start(s), end(e), generations(g) {}
};

// ────────────────────────────────────────────────────────────────────
void print_usage() {
    std::cerr <<
        "Usage: mutator -fasta <in.fa> -rate <u> -generations <G> "
        "-mode <0|1|2|3> -threads <T> -seed <S> -out_prefix <P> "
        "[-bed <regions.bed>] [-TsTv <ratio>]\n"
        "  -TsTv: transition/transversion ratio (default 1.0)\n";
}

bool parse_args(int argc, char* argv[],
                std::string &fasta, double &mu,
                int &G, int &mode, int &threads,
                int &seed, std::string &prefix,
                std::string &bedfile, double &tsTv)
{
    if (argc < 15 || argc % 2 == 0) {
        print_usage();
        return false;
    }
    tsTv = 1.0;  // default
    for (int i = 1; i + 1 < argc; i += 2) {
        std::string a = argv[i], v = argv[i+1];
        if      (a == "-fasta")       fasta   = v;
        else if (a == "-rate")        mu      = std::stod(v);
        else if (a == "-generations") G       = std::stoi(v);
        else if (a == "-mode")        mode    = std::stoi(v);
        else if (a == "-threads")     threads = std::stoi(v);
        else if (a == "-seed")        seed    = std::stoi(v);
        else if (a == "-out_prefix")  prefix  = v;
        else if (a == "-bed")         bedfile = v;
        else if (a == "-TsTv")        tsTv    = std::stod(v);
        else {
            std::cerr << "Unknown argument: " << a << '\n';
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

// ────────────────────────────────────────────────────────────────────
// BED loader
// ────────────────────────────────────────────────────────────────────
bool load_bed(const std::string &bedfile,
              std::unordered_map<std::string, std::vector<BedRegion>> &bed_map)
{
    std::ifstream in(bedfile);
    if (!in) { std::cerr << "Cannot open BED file: " << bedfile << '\n'; return false; }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0]=='#') continue;
        std::istringstream ss(line);
        BedRegion r; std::string na, strand;
        if (!(ss >> r.chrom >> r.start >> r.end >> na >> na >> strand)) {
            std::cerr << "Skipping malformed BED line: " << line << '\n';
            continue;
        }
        std::getline(ss, r.extra);
        if (r.end <= r.start) continue;
        bed_map[r.chrom].push_back(r);
    }
    for (auto &kv : bed_map) {
        auto &v = kv.second;
        std::sort(v.begin(), v.end(), [](auto &a, auto &b){ return a.start < b.start; });
    }
    return true;
}

// ────────────────────────────────────────────────────────────────────
// FASTA helpers
// ────────────────────────────────────────────────────────────────────
bool load_fasta(const std::string &fasta,
                std::vector<std::pair<std::string,std::string>> &seqs)
{
    std::ifstream in(fasta);
    if (!in) { std::cerr << "Cannot open FASTA file: " << fasta << '\n'; return false; }
    std::string line, hdr, seq;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0]=='>') {
            if (!hdr.empty()) seqs.emplace_back(hdr, seq);
            hdr = line.substr(1); seq.clear();
        } else {
            seq += line;
        }
    }
    if (!hdr.empty()) seqs.emplace_back(hdr, seq);
    return true;
}

bool write_fasta(const std::string &out,
                 const std::vector<std::pair<std::string,std::string>> &seqs)
{
    std::ofstream o(out);
    if (!o) { std::cerr << "Cannot write FASTA file: " << out << '\n'; return false; }
    for (auto &p : seqs) {
        o << '>' << p.first << '\n';
        auto &s = p.second;
        for (size_t i = 0; i < s.size(); i += 60)
            o << s.substr(i, 60) << '\n';
    }
    return true;
}

// ────────────────────────────────────────────────────────────────────
// Substitution map
// ────────────────────────────────────────────────────────────────────
const std::unordered_map<char,std::vector<char>>& get_subs() {
    static const std::unordered_map<char,std::vector<char>> m = {
        {'A',{'C','G','T'}}, {'C',{'A','G','T'}},
        {'G',{'A','C','T'}}, {'T',{'A','C','G'}},
        {'a',{'c','g','t'}}, {'c',{'a','g','t'}},
        {'g',{'a','c','t'}}, {'t',{'a','c','g'}}
    };
    return m;
}

// Transition helper (case-insensitive compare, but preserve case in output)
inline bool is_transition(char a, char b) {
    a = std::toupper(static_cast<unsigned char>(a));
    b = std::toupper(static_cast<unsigned char>(b));
    return (a=='A'&&b=='G')||(a=='G'&&b=='A')||(a=='C'&&b=='T')||(a=='T'&&b=='C');
}

// For deterministic transition choice preserving case
inline char transition_of(char base) {
    switch (base) {
        case 'A': return 'G';
        case 'G': return 'A';
        case 'C': return 'T';
        case 'T': return 'C';
        case 'a': return 'g';
        case 'g': return 'a';
        case 'c': return 't';
        case 't': return 'c';
        default:  return base;
    }
}

// ────────────────────────────────────────────────────────────────────
// Decay-distributed generation sampler
// ────────────────────────────────────────────────────────────────────
inline int sample_decayG(int G_user, std::mt19937 &rng) {
    static constexpr double LAMBDA = 1.4e-6;
    static constexpr int MAX_DECAY = 3'000'000;
    int G = std::clamp(G_user, 0, MAX_DECAY);
    if (G==0) return 0;
    std::uniform_real_distribution<> U(0.0,1.0);
    double u = U(rng);
    double denom = 1.0 - std::exp(-LAMBDA*G);
    double x = -std::log(1.0 - u*denom)/LAMBDA;
    return static_cast<int>(std::round(x));
}

// ────────────────────────────────────────────────────────────────────
// Main
// ────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    std::string fasta, prefix, bedfile;
    double mu = 0.0, tsTv = 1.0;
    int G=0, mode=0, threads=1, seed=0;

    if (!parse_args(argc, argv, fasta, mu, G, mode, threads, seed, prefix, bedfile, tsTv))
        return EXIT_FAILURE;

    std::cout << "Using Ts/Tv ratio: " << tsTv << "\n";

    if (threads<1) threads=1;
    int max_threads = omp_get_max_threads();
    if (threads>max_threads) {
        std::cerr<<"Requested "<<threads<<" threads, but only "
                 <<max_threads<<" available. Using "<<max_threads<<".\n";
        threads = max_threads;
    }
    omp_set_dynamic(0);
    omp_set_num_threads(threads);

    std::unordered_map<std::string,std::vector<BedRegion>> bed_map;
    if (mode>=2) {
        std::cout<<"Loading BED "<<bedfile<<" …\n";
        if (!load_bed(bedfile, bed_map)) return EXIT_FAILURE;
        std::cout<<"  "<<bed_map.size()<<" chromosomes indexed\n";
    }

    std::cout<<"Loading FASTA "<<fasta<<" …\n";
    auto t0 = std::chrono::high_resolution_clock::now();
    std::vector<std::pair<std::string,std::string>> seqs;
    if (!load_fasta(fasta, seqs)) return EXIT_FAILURE;
    size_t N = seqs.size();
    if (N==0) { std::cerr<<"No sequences loaded.\n"; return EXIT_FAILURE; }
    std::cout<<"  "<<N<<" sequences\n";

    unsigned long long genome_size = 0;
    for (auto &p: seqs) genome_size += p.second.size();
    const auto &subs = get_subs();

    // Outputs and counters
    std::vector<std::pair<std::string,std::string>> out_seqs(N);
    unsigned long long total_muts = 0, total_ts = 0, total_tv = 0;

    // New exact counters for recurrence
    unsigned long long total_unique_sites_hit = 0; // sites with >=1 successful mutation
    unsigned long long total_recurrent_hits   = 0; // sum over sites of (hits - 1) when hits>=1

    // ───────── parallel mutator ─────────
    #pragma omp parallel reduction(+:total_muts,total_ts,total_tv,total_unique_sites_hit,total_recurrent_hits)
    {
        std::mt19937 rng(seed + 7*omp_get_thread_num());
        std::uniform_real_distribution<> U(0.0,1.0);

        #pragma omp for schedule(static)
        for (size_t i=0; i<N; ++i) {
            auto hdr = seqs[i].first;
            auto seq = seqs[i].second;
            size_t L  = seq.size();

            // Per-sequence hit counter to track recurrence (only for successful mutations)
            std::vector<uint16_t> pos_hits(L, 0);

            // build segments
            std::vector<Segment> segs;
            if (mode<2) {
                int g = (mode==0)?G:sample_decayG(G, rng);
                if (L) segs.emplace_back(0,L,g);
            } else {
                size_t prev = 0;
                auto it = bed_map.find(hdr);
                if (it!=bed_map.end()) {
                    for (auto &r: it->second) {
                        if (r.start>prev && prev<L)
                            segs.emplace_back(prev, std::min(r.start,L), G);
                        int g = (mode==2)
                              ? ((r.extra.rfind("gene",0)==0)?G:sample_decayG(G,rng))
                              : ((r.extra.find("new_TE")!=std::string::npos)
                                 ? sample_decayG(G,rng) : G);
                        if (r.start<L && r.end>r.start)
                            segs.emplace_back(r.start, std::min(r.end,L), g);
                        prev = std::min(r.end,L);
                    }
                }
                if (prev<L) segs.emplace_back(prev,L,G);
            }

            unsigned long long seq_muts = 0;

            for (auto &S: segs) {
                if (S.end<=S.start) continue;
                double lambda_segment = double(S.generations)*mu*double(S.end-S.start);
                std::poisson_distribution<> pd(lambda_segment);
                int mcount = pd(rng);

                for (int m=0; m<mcount; ++m) {
                    std::uniform_int_distribution<size_t> pickPos(S.start, S.end-1);
                    size_t pos = pickPos(rng);
                    char orig = seq[pos];
                    auto it2 = subs.find(orig);
                    if (it2==subs.end()) continue; // skip non-ACGT; don't count as a mutation

                    // True Ts vs Tv decision:
                    double p_ts = tsTv/(tsTv+1.0);
                    char nb;
                    if (U(rng) < p_ts) {
                        // transition
                        nb = transition_of(orig);
                        ++total_ts;
                    } else {
                        // transversion: pick uniformly among the two non-transition cands
                        const auto &cands = it2->second;
                        char tvs[2]; int idx=0;
                        for (char c: cands)
                            if (!is_transition(orig,c))
                                tvs[idx++] = c;
                        std::uniform_int_distribution<int> pickTv(0,1);
                        nb = tvs[pickTv(rng)];
                        ++total_tv;
                    }

                    // apply mutation
                    seq[pos] = nb;
                    ++seq_muts;
                    // track recurrence (successful hit only)
                    if (pos_hits[pos] < std::numeric_limits<uint16_t>::max())
                        ++pos_hits[pos];
                }
            }

            // tally unique and recurrent hits for this sequence
            unsigned long long unique_sites = 0;
            unsigned long long recurrent_hits = 0;
            for (size_t p = 0; p < L; ++p) {
                uint16_t h = pos_hits[p];
                if (h > 0) {
                    ++unique_sites;           // site hit at least once
                    if (h > 1) recurrent_hits += (h - 1); // extra hits beyond the first
                }
            }

            total_unique_sites_hit += unique_sites;
            total_recurrent_hits   += recurrent_hits;
            total_muts             += seq_muts;
            out_seqs[i] = { std::move(hdr), std::move(seq) };

            if ((i+1)%10000==0 && omp_get_thread_num()==0)
                std::cout<<"  processed "<<(i+1)<<"/"<<N<<" sequences\n";
        }
    }

    // ───────── write outputs ─────────
    std::string outfa = prefix + ".fa";
    std::cout<<"Writing "<<outfa<<" …\n";
    if (!write_fasta(outfa, out_seqs)) return EXIT_FAILURE;

    std::string outtxt = prefix + ".txt";
    std::ofstream ot(outtxt);
    if (!ot) { std::cerr<<"Cannot write metrics file: "<<outtxt<<"\n"; return EXIT_FAILURE; }
    
    // Exact (non-recurrent) and total rates
    // unique_sites_hit == number of sites with >=1 successful mutation
    // Non-recurrent mutations/site counts only the first hit to each site
    const double muts_per_site_exact_nonrec = (genome_size > 0)
        ? static_cast<double>(total_unique_sites_hit) / static_cast<double>(genome_size)
        : 0.0;
    // Total mutations per site (requested behavior for "Mutions / site")
    const double muts_per_site_total = (genome_size > 0)
        ? static_cast<double>(total_muts) / static_cast<double>(genome_size)
        : 0.0;

    // Output summary (with the requested field names/spelling)
    ot << "Genome size:  " << genome_size << "\n"
       << "Total mutations: " << total_muts << "\n"
       << "Recurrent mutations: " << total_recurrent_hits << "\n"
       << "Mutions / site:    " << std::setprecision(10) << muts_per_site_total << "\n"
       << "Mutions / site * 2:    " << std::setprecision(10) << (muts_per_site_total * 2.0) << "\n"
       << "Non-recurrent Mutions / site:    " << std::setprecision(10) << muts_per_site_exact_nonrec << "\n"
       << "Non-recurrent Mutions / site * 2:    " << std::setprecision(10) << (muts_per_site_exact_nonrec * 2.0) << "\n";

    // Debug prints
    std::cout<<"Transition count:      "<<total_ts<<"\n";
    std::cout<<"Transversion count:    "<<total_tv<<"\n";
    std::cout<<"Ts + Tv = "<<(total_ts+total_tv)
             <<" (should equal Total mutations: "<<total_muts<<")\n";
    std::cout<<"Unique sites hit:      "<<total_unique_sites_hit<<"\n";
    std::cout<<"Recurrent mutations:   "<<total_recurrent_hits<<"\n";

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt = t1 - t0;
    std::cout<<"Done in "<<dt.count()<<" s\n";
    return EXIT_SUCCESS;
}
