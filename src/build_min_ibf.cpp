#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <iostream>
#include <minimizer.hpp>


struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

template <class IBFType>
void fill_ibf(IBFType & ibf,
              std::vector<std::filesystem::path> const & dir_path,
              uint8_t const n_k,
              uint64_t const n_w,
              uint64_t const n_bins,
              std::string const & extension)
{
    minimizer mini{window{n_w}, kmer{n_k}};

    for (uint64_t cur_bin = 0; cur_bin < n_bins; ++cur_bin)
    {
        std::filesystem::path bin_path{dir_path[0]};
        bin_path /= ("bin_" + std::to_string(cur_bin) + extension);

        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{bin_path};

        for (auto & [seq] : fin)
        {
            mini.compute(seq);
            for (auto && hash : mini.minimizer_hash)
                ibf.emplace(hash, seqan3::bin_index{cur_bin});
        }
    }
}

// Use, when minimisers were precalculated and can be found in a binary file for each experiment
template <class IBFType>
void fill_ibf(IBFType & ibf,
              std::vector<std::filesystem::path> const & in_path)
{
    std::ifstream infile;

    for (uint64_t cur_bin = 0; cur_bin < in_path.size(); ++cur_bin)
    {
        infile.open(in_path[cur_bin], std::ios::binary);
        if (!infile.is_open()) {
               std::cerr << "Error in open file.\n";
               break;
        }
        uint64_t num;
        while(infile.read((char*)&num, sizeof(num)))
            ibf.emplace(num, seqan3::bin_index{cur_bin});
        infile.close();
    }
}

void run_program(std::vector<std::filesystem::path> const & dir_path,
                 std::filesystem::path const & out_path,
                 uint8_t const n_k,
                 uint64_t const n_w,
                 uint64_t const n_bins,
                 uint64_t const n_bits,
                 uint64_t const n_hash,
                 bool const gz,
                 bool const bz2,
                 bool precalculated)
{
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>  ibf{seqan3::bin_count{n_bins},
                                         seqan3::bin_size{n_bits},
                                         seqan3::hash_function_count{n_hash}};

    std::string extension{".fasta"};
    if (gz)
        extension += ".gz";
    if (bz2)
        extension += ".bz2";

    if (precalculated)
        fill_ibf(ibf, dir_path);
    else
        fill_ibf(ibf, dir_path, n_k, n_w, n_bins, extension);

    std::ofstream os{out_path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(ibf);
}

void get_minimizer_files(std::vector<std::filesystem::path> const & in_paths,
                         std::filesystem::path const & out_path,
                         uint8_t const n_k,
                         uint64_t const n_w)
{
    std::unordered_map<uint64_t, uint8_t> hash_table{}; // storage for minimizers
    std::ofstream outfile;
    std::ofstream headerfile;
    uint64_t count{0};
    uint64_t filesize{0};
    uint16_t cutoff{50};
    // Cutoffs and bounds from Mantis
    std::vector<int> cutoffs{1,3,10,20};
    std::vector<uint64_t> cutoff_bounds{314572800, 524288000, 1073741824, 3221225472};

    minimizer mini{window{n_w}, kmer{n_k}};

    for(auto & file: in_paths)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq>> fin{file};

        for (auto & [seq] : fin)
        {
            mini.compute(seq);
            for (auto && hash : mini.minimizer_hash)
                hash_table[hash] = std::min<uint8_t>(126u, hash_table[hash]+1);
        }

        filesize = std::filesystem::file_size(file);
        filesize = filesize + filesize; // Filesize multiplied by two, because Mantis based on fastq files, we use fasta files

        for(unsigned k = 0; k < cutoff_bounds.size(); ++k)
        {
            if (filesize <= cutoff_bounds[k])
            {
                cutoff = cutoffs[k];
                break;
            }
        }

        outfile.open(std::string(out_path) + std::string(file.stem()) + ".minimizer", std::ios::binary);
        {
            for (auto & hash : hash_table)
            {
                if (hash.second > cutoff)
                {
                    outfile.write((char*) & hash.first, sizeof(hash.first));
                    ++count;
                }

            }
        }
        outfile.close();

        headerfile.open(std::string(out_path) + std::string(file.stem()) + ".header");
        headerfile << (uint64_t) n_k << "\t" << n_w << "\t" << cutoff << "\t" << count << "\n";
        headerfile.close();
        count = 0;
        cutoff = 50;
        hash_table.clear();
    }
}

struct cmd_arguments
{
    std::vector<std::filesystem::path> bin_path{};
    std::filesystem::path out_path{"./"};
    uint64_t w{23};
    uint8_t k{20};
    uint64_t bins{64};
    uint64_t bits{4096};
    uint64_t hash{2};
    bool gz{false};
    bool bz2{false};
    bool binary{false};
    bool precalculated{false};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Build an IBF using minimizers.";
    parser.info.version = "1.0.0";
    parser.add_positional_option(args.bin_path, "Please provide one path to a directory containing one FASTA file for "
                                                "each bin. Or provide a list of binary minimizer files as outputed by "
                                                "build_min_ibf --binary.");
    parser.add_option(args.out_path, 'o', "out", "Please provide a valid output path.", seqan3::option_spec::DEFAULT);
    parser.add_option(args.w, '\0', "window", "Choose the window size.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 1000});
    parser.add_option(args.k, '\0', "kmer", "Choose the kmer size.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 32});
    parser.add_option(args.bins, '\0', "bins", "Choose the number of bins.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 65536});
    parser.add_option(args.bits, '\0', "bits", "Choose the size in bits of one bin.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 35184372088832});
    parser.add_option(args.hash, '\0', "hash", "Choose the number of hashes.", seqan3::option_spec::DEFAULT,
                      seqan3::arithmetic_range_validator{1, 4});
    parser.add_flag(args.gz, '\0', "gz", "Expect FASTA files to be gz compressed.");
    parser.add_flag(args.bz2, '\0', "bz2", "Expect FASTA files to be bz2 compressed.");
    parser.add_flag(args.binary, '\0', "binary", "Set to true when not IBF should be build but only binary minimizer files.");
    parser.add_flag(args.precalculated, '\0', "pre", "Set to true when precalculated binary minimizer files are used.");
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"build_min_ibf", argc, argv, false};
    cmd_arguments args{};
    initialize_argument_parser(myparser, args);
    try
    {
         myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cout << "[Error] " << ext.what() << "\n";
        return -1;
    }

    if (args.binary)
    {
        get_minimizer_files(args.bin_path, args.out_path, args.k, args.w);
        return 0;
    }

    if (args.gz && args.bz2)
        throw seqan3::argument_parser_error{"Files cannot be both gz and bz2 compressed."};

    if (args.k > args.w)
        throw seqan3::argument_parser_error{"The kmer size cannot be bigger than the window size."};

    run_program(args.bin_path, args.out_path, args.k, args.w, args.bins, args.bits, args.hash, args.gz, args.bz2,
                args.precalculated);
    return 0;
}
