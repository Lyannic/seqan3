// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/kmer_delete_mask_hash.hpp>
#include <seqan3/range/views/kmer_nthash.hpp>
#include <seqan3/range/views/kmer_gapped_hash.hpp>
#include <seqan3/range/views/kmer_gapped_delete_mask_hash.hpp>
#include <seqan3/range/views/kmer_gapped_perfect_hash.hpp>
#include <seqan3/range/views/kmer_gapped_delete_mask_perfect_hash.hpp>
#include <seqan3/range/views/kmer_issh_hash.hpp>
#include <seqan3/range/views/kmer_issh_precomp_hash.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/naive_kmer_hash.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/contrib/issh/SpacedQmer.hpp>

#ifdef SEQAN3_HAS_SEQAN2
    #include <seqan/index.h>
#endif // SEQAN3_HAS_SEQAN2

inline benchmark::Counter bp_per_second(size_t const basepairs)
{
    return benchmark::Counter(basepairs,
                              benchmark::Counter::kIsIterationInvariantRate,
                              benchmark::Counter::OneK::kIs1000);
}

inline seqan3::shape make_gapped_shape(size_t const k)
{
    seqan3::shape shape_{};

    for (size_t i{0}; i < k - 1; ++i)
        shape_.push_back((i + 1) % 2);

    shape_.push_back(1u);
    shape_.push_back(0u);
    return shape_;
}


static void arguments(benchmark::internal::Benchmark* b)
{
    for (int32_t sequence_length : {1'000, 50'000, /*1'000'000*/})
    {
        for (int32_t k : {8, /*16, 24,*/ 30})
        {
            b->Args({sequence_length, k});
        }
    }
}

static void seqan_kmer_hash_ungapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | seqan3::views::kmer_hash(seqan3::ungapped{static_cast<uint8_t>(k)}))
            benchmark::DoNotOptimize(sum += h);
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void seqan_kmer_delete_mask_hash_ungapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = test::generate_sequence<dna4>(sequence_length, 0, 0);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | views::kmer_delete_mask_hash(ungapped{static_cast<uint8_t>(k)}))
            benchmark::DoNotOptimize(sum += h);
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void seqan_kmer_nthash_ungapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = test::generate_sequence<dna4>(sequence_length, 0, 0);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | views::kmer_nthash(ungapped{static_cast<uint8_t>(k)}))
            benchmark::DoNotOptimize(sum += h);
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void nthash_ungapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto generatedSeq = test::generate_sequence<dna4>(sequence_length, 0, 0);
    auto genereatedSeqVec = generatedSeq | views::to_char;
    std::string seq(genereatedSeqVec.begin(), genereatedSeqVec.end());

    volatile size_t sum{0};

    for (auto _ : state)
    {
        uint64_t hVal;
        hVal = NTF64(seq.c_str(), k);
        benchmark::DoNotOptimize(sum += hVal);
        for (size_t i = 1; i < seq.length() - k + 1; i++) {
            hVal = NTF64(hVal, k, seq[i-1], seq[i-1+k]);
            benchmark::DoNotOptimize(sum += hVal);
        }
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void seqan_kmer_hash_gapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | seqan3::views::kmer_hash(make_gapped_shape(k)))
            benchmark::DoNotOptimize(sum += h);
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void seqan_kmer_gapped_hash_gapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = test::generate_sequence<dna4>(sequence_length, 0, 0);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | views::kmer_gapped_hash(make_gapped_shape(k)))
            benchmark::DoNotOptimize(sum += h);
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void seqan_kmer_gapped_delete_mask_hash_gapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = test::generate_sequence<dna4>(sequence_length, 0, 0);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | views::kmer_gapped_delete_mask_hash(make_gapped_shape(k)))
            benchmark::DoNotOptimize(sum += h);
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void seqan_kmer_gapped_perfect_hash_gapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = test::generate_sequence<dna4>(sequence_length, 0, 0);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | views::kmer_gapped_perfect_hash(make_gapped_shape(k)))
            benchmark::DoNotOptimize(sum += h);
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void seqan_kmer_gapped_delete_mask_perfect_hash_gapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = test::generate_sequence<dna4>(sequence_length, 0, 0);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | views::kmer_gapped_delete_mask_perfect_hash(make_gapped_shape(k)))
            benchmark::DoNotOptimize(sum += h);
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void seqan_kmer_issh_hash_gapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = test::generate_sequence<dna4>(sequence_length, 0, 0);

    auto shape = make_gapped_shape(k);

    std::string shape_string = "";
    for (size_t i{0}; i < shape.size() - 1u; ++i)
    {
        if (shape[i])
        {
            shape_string += "1";
        } else
        {
            shape_string += "0";
        }
    }
    if (shape[shape.size() - 1u])
    {
        shape_string += "1";
    }
    else 
    {
        shape_string += "0";
    }

    auto spaced_qmer = SpacedQmer(shape_string, 0);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | views::kmer_issh_hash(shape, spaced_qmer.GetMultipleShifts(), spaced_qmer.GetPosOne()))
            benchmark::DoNotOptimize(sum += h);
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void seqan_kmer_issh_hash_gapped_with_precomputations(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = test::generate_sequence<dna4>(sequence_length, 0, 0);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | views::kmer_issh_precomp_hash(make_gapped_shape(k)))
            benchmark::DoNotOptimize(sum += h);
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void naive_kmer_hash(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | seqan3::views::naive_kmer_hash(k))
            benchmark::DoNotOptimize(sum += h);
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

#ifdef SEQAN3_HAS_SEQAN2
inline auto make_gapped_shape_seqan2(size_t const k)
{
    seqan::String<char> bitmap;

    for (size_t i{0}; i < k - 1; ++i)
        seqan::append(bitmap, seqan::CharString(std::to_string((i + 1) % 2)));

    seqan::append(bitmap, seqan::CharString("1"));

    seqan::Shape<seqan::Dna, seqan::GenericShape> s(bitmap);
    return s;
}

static void seqan2_kmer_hash_ungapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = seqan3::test::generate_sequence_seqan2<seqan::Dna>(sequence_length, 0, 0);
    seqan::Shape<seqan::Dna, seqan::SimpleShape> s;
    seqan::resize(s, k);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        auto it = seqan::begin(seq);
        seqan::hashInit(s, it);
        for (size_t i{0}; i < seqan::length(seq) - k + 1; ++i, ++it)
        {
            benchmark::DoNotOptimize(sum += seqan::hashNext(s, it));
        }
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void seqan2_kmer_hash_gapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = seqan3::test::generate_sequence_seqan2<seqan::Dna>(sequence_length, 0, 0);
    seqan::Shape<seqan::Dna, seqan::GenericShape> s = make_gapped_shape_seqan2(k);

    volatile size_t sum{0};

    for (auto _ : state)
    {
        auto it = seqan::begin(seq);
        seqan::hashInit(s, it);
        for (size_t i{0}; i < seqan::length(seq) - k + 1; ++i, ++it)
        {
            benchmark::DoNotOptimize(sum += seqan::hashNext(s, it));
        }
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

BENCHMARK(seqan2_kmer_hash_ungapped)->Apply(arguments);
BENCHMARK(seqan2_kmer_hash_gapped)->Apply(arguments);
#endif // SEQAN3_HAS_SEQAN2

BENCHMARK(seqan_kmer_hash_ungapped)->Apply(arguments);
BENCHMARK(seqan_kmer_delete_mask_hash_ungapped)->Apply(arguments);
BENCHMARK(seqan_kmer_nthash_ungapped)->Apply(arguments);
BENCHMARK(nthash_ungapped)->Apply(arguments);
BENCHMARK(seqan_kmer_hash_gapped)->Apply(arguments);
BENCHMARK(seqan_kmer_gapped_hash_gapped)->Apply(arguments);
BENCHMARK(seqan_kmer_gapped_delete_mask_hash_gapped)->Apply(arguments);
BENCHMARK(seqan_kmer_gapped_perfect_hash_gapped)->Apply(arguments);
BENCHMARK(seqan_kmer_gapped_delete_mask_perfect_hash_gapped)->Apply(arguments);
BENCHMARK(seqan_kmer_issh_hash_gapped)->Apply(arguments);
BENCHMARK(seqan_kmer_issh_hash_gapped_with_precomputations)->Apply(arguments);
BENCHMARK(naive_kmer_hash)->Apply(arguments);

BENCHMARK_MAIN();
