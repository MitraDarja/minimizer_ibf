// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/eseiler/minimizer_ibf/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides minimizer.
 */

#pragma once

#include <seqan/seq_io.h>

#include <strong_types.hpp>

//!\brief Whether to use xor for computing the hash value.
enum use_xor : bool
{
    //!\brief Do not use xor.
    no,
    //!\brief Use xor.
    yes
};

template<use_xor do_xor = use_xor::yes>
struct minimizer
{
private:
    //!\brief The alphabet type.
    using alphabet_t = seqan::Dna;
    //!\brief The text type.
    using text_t = seqan::String<alphabet_t>;
    //!\brief The type of the complemented text.
    using complement_t = seqan::ModifiedString<text_t const, seqan::ModComplementDna>;
    //!\brief The type of the reverse complemented text.
    using reverse_complement_t = seqan::ModifiedString<complement_t, seqan::ModReverse>;
    //!\brief The seqan::Shape type to compute rolling hashes with.
    using shape_t = seqan::Shape<alphabet_t, seqan::SimpleShape>;

    //!\brief The window size of the minimizer.
    uint64_t w{};
    //!\brief The size of the k-mers.
    uint64_t k{};
    //!\brief Random but fixed value to xor k-mers with. Counteracts consecutive minimizers.
    uint64_t seed{};

    //!\brief Shape for computing the forward strand k-mers.
    shape_t forward_shape{};
    //!\brief Shape for computing the reverse strand k-mers.
    shape_t reverse_shape{};

    //!\brief Stores the k-mer hashes of the forward strand.
    std::vector<uint64_t> forward_hashes;
    //!\brief Stores the k-mer hashes of the reverse complement strand.
    std::vector<uint64_t> reverse_hashes;


public:

    //!\brief Stores the hashes of the minimizers.
    std::vector<uint64_t> minimizer_hash;
    //!\brief Stores the begin positions of the minimizers.
    std::vector<uint64_t> minimizer_begin;
    //!\brief Stores the end positions of the minimizers.
    std::vector<uint64_t> minimizer_end;

    minimizer() = default;                                        //!< Defaulted
    minimizer(minimizer const &) = default;                       //!< Defaulted
    minimizer(minimizer &&) = default;                            //!< Defaulted
    minimizer & operator=(minimizer const &) = default;           //!< Defaulted
    minimizer & operator=(minimizer &&) = default;                //!< Defaulted
    ~minimizer() = default;                                       //!< Defaulted

    /*!\brief Constructs a minimizer from given k-mer, window size and a seed.
     * \param[in] w_    The window size.
     * \param[in] k_    The k-mer size.
     * \param[in] seed_ The seed to use. Default: 0x8F3F73B5CF1C9ADE.
     */
    minimizer(window const w_, kmer const k_, uint64_t const seed_ = 0x8F3F73B5CF1C9ADE) :
        w{w_.v}, k{k_.v}, seed{seed_}
    {
        seqan::resize(forward_shape, k);
        seqan::resize(reverse_shape, k);
    }

    /*!\brief Resize the minimizer.
     * \param[in] w_    The new window size.
     * \param[in] k_    The new k-mer size.
     * \param[in] seed_ The new seed to use. Default: 0x8F3F73B5CF1C9ADE.
     */
    void resize(window const w_, kmer const k_, uint64_t const seed_ = 0x8F3F73B5CF1C9ADE)
    {
        w = w_.v;
        k = k_.v;
        seed = seed_;
        seqan::resize(forward_shape, k);
        seqan::resize(reverse_shape, k);
    }

    void compute(text_t const & text)
    {
        uint64_t text_length = seqan::length(text);

        forward_hashes.clear();
        reverse_hashes.clear();
        minimizer_hash.clear();
        minimizer_begin.clear();
        minimizer_end.clear();

        // Return empty vector if text is shorter than k.
        if (k > text_length)
            return;

        reverse_complement_t rc_text{text}; // This does not copy.
        uint64_t possible_minimizers = text_length > w ? text_length - w + 1u : 1u;
        uint64_t possible_kmers = text_length - k + 1;
        uint64_t kmers_per_window = w - k + 1u;

        // Compute all k-mer hashes for both forward and reverse strand.

        // Helper lambda for xor'ing values depending on `do_xor`.
        auto hash_impl = [this] (uint64_t const val)
        {
            if constexpr(do_xor)
                return val ^ seed;
            else
                return val;
        };

        forward_hashes.reserve(possible_kmers);
        reverse_hashes.reserve(possible_kmers);
        seqan::hashInit(forward_shape, seqan::begin(text));
        seqan::hashInit(reverse_shape, seqan::begin(rc_text));

        for (uint64_t i = 0; i < possible_kmers; ++i)
        {
            forward_hashes.push_back(hash_impl(seqan::hashNext(forward_shape, seqan::begin(text) + i)));
            reverse_hashes.push_back(hash_impl(seqan::hashNext(reverse_shape, seqan::begin(rc_text) + i)));
        }

        // Choose the minimizers.
        minimizer_hash.reserve(possible_minimizers);
        minimizer_begin.reserve(possible_minimizers);
        minimizer_end.reserve(possible_minimizers);

        // Stores hash, begin and end for all k-mers in the window
        std::deque<std::tuple<uint64_t, uint64_t, uint64_t>> window_values;

        // Initialisation. We need to compute all hashes for the first window.
        for (uint64_t i = 0; i < kmers_per_window; ++i)
        {
            // Get smallest canonical k-mer.
            uint64_t forward_hash = forward_hashes[i];
            uint64_t reverse_hash = reverse_hashes[possible_kmers - i - 1];
            window_values.emplace_back(std::min(forward_hash, reverse_hash), i, i + k - 1);
        }

        auto min = std::min_element(std::begin(window_values), std::end(window_values));
        minimizer_hash.push_back(std::get<0>(*min));
        minimizer_begin.push_back(std::get<1>(*min));
        minimizer_end.push_back(std::get<2>(*min));


        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        bool minimizer_changed{false};
        for (uint64_t i = 1; i < possible_minimizers; ++i)
        {
            // Shift the window.
            // If current minimizer leaves the window, we need to decide on a new one.
            if (min == std::begin(window_values))
            {
                window_values.pop_front();
                min = std::min_element(std::begin(window_values), std::end(window_values));
                minimizer_changed = true;
            }
            else
            {
                window_values.pop_front();
            }

            uint64_t forward_hash = forward_hashes[kmers_per_window - 1 + i];
            uint64_t reverse_hash = reverse_hashes[possible_kmers - kmers_per_window - i];
            window_values.emplace_back(std::min(forward_hash, reverse_hash),
                                       kmers_per_window + i - 1,
                                       kmers_per_window + i + k - 2);

            if (std::get<0>(window_values.back()) < std::get<0>(*min))
            {
                min = std::prev(std::end(window_values));
                minimizer_changed = true;
            }

            if (minimizer_changed)
            {
                minimizer_hash.push_back(std::get<0>(*min));
                minimizer_begin.push_back(std::get<1>(*min));
                minimizer_end.push_back(std::get<2>(*min));
                minimizer_changed = false;
            }
        }
        return;
    }

};
