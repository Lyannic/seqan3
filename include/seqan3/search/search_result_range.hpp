// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::search_result_range.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/range/views/single_pass_input.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief An input range over the search results generated by the given search algorithm.
 * \ingroup search
 * \implements std::ranges::input_range
 *
 * \tparam search_algorithm_t The type of the search algorithm; must model std::invocable with the value type of
 *                            the given range query over queries.
 * \tparam query_range_t The type of the range over the queries; must model std::ranges::view and
 *                       std::ranges::forward_range.
 *
 * \details
 *
 * Provides a lazy input-range interface over the search results generated by the underlying search algorithm.
 */
template <typename search_algorithm_t, std::ranges::view query_range_t>
//!\cond
    requires std::ranges::forward_range<query_range_t> &&
             std::invocable<search_algorithm_t, std::ranges::range_reference_t<query_range_t>>
//!\endcond
class search_result_range
{
private:
    class search_result_range_iterator;

    //!\brief The type of a single query.
    using query_t = std::ranges::range_value_t<query_range_t>;
    //!\brief The type of the buffer returned by invoking the search on a single query.
    using hit_range_t = std::invoke_result_t<search_algorithm_t, query_t &>;
    //!\brief The value type of the buffer returned by invoking the search on a single query.
    using result_buffer_value_t = std::pair<size_t, std::ranges::range_value_t<hit_range_t>>;
    //!\brief The type of the buffer used to store the hits.
    using result_buffer_t = std::vector<result_buffer_value_t>;
    //!\brief The wrapped query_range_t to add single pass behaviour.
    using single_pass_query_range_t = decltype(views::zip(std::views::iota(0), std::declval<query_range_t &&>())
                                               | views::single_pass_input);
    //!\brief The iterator over the wrapped query range.
    using single_pass_query_range_iterator = std::ranges::iterator_t<single_pass_query_range_t>;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    search_result_range() = default; //!< Defaulted.
    search_result_range(search_result_range const &) = delete; //!< This is a move-only type.
    search_result_range(search_result_range &&) = default; //!< Defaulted.
    search_result_range & operator=(search_result_range const &) = delete; //!< This is a move-only type.
    search_result_range & operator=(search_result_range &&) = default; //!< Defaulted.
    ~search_result_range() = default; //!< Defaulted.

    /*!\brief Constructs a search result range.
     *
     * \param[in] search_algorithm The search algorithm to invoke on a single query.
     * \param[in] query_range The range over the queries to be searched.
     */
    search_result_range(search_algorithm_t search_algorithm, query_range_t query_range) :
        search_algorithm{std::move(search_algorithm)},
        single_pass_query_range{views::zip(std::views::iota(0), std::move(query_range))},
        current_query_it{single_pass_query_range.begin()}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */

    /*!\brief Returns an iterator to the first element of the search result range.
     * \return An iterator to the first element.
     *
     * \details
     *
     * The fist invocation of this function will invoke the search for the first query in the underlying range.
     */
    constexpr search_result_range_iterator begin()
    {
        search_result_range_iterator tmp_it{*this, first_call_of_begin};
        first_call_of_begin = false;
        return tmp_it;
    }

    search_result_range_iterator begin() const = delete;
    search_result_range_iterator cbegin() const = delete;

    /*!\brief Returns a sentinel signaling the end of the search result range.
     * \return a sentinel.
     *
     * \details
     *
     * The search result range is an input range and the end is reached when all queries have been processed.
     */
    constexpr std::ranges::default_sentinel_t end() noexcept
    {
        return {};
    }

    constexpr std::ranges::default_sentinel_t end() const = delete;
    constexpr std::ranges::default_sentinel_t cend() const = delete;
    //!\}

private:
    /*!\brief Receives the search results from the next query.
     *
     * \returns `true` if hits could be found for the current query, otherwise `false`.
     *
     * \details
     *
     * Invokes the search on the current query and increments the query iterator.
     */
    bool next()
    {
        result_buffer.clear();

        if (at_end())
            return false;

        auto && [query_id, query] = *current_query_it;
        for (auto res : search_algorithm(query))
            result_buffer.emplace_back(query_id, std::move(res));

        ++current_query_it;
        return !std::ranges::empty(result_buffer);
    }

    //!\brief Returns `true` if all queries in the underlying query range have been processed, otherwise false.
    bool at_end() noexcept
    {
        return current_query_it == single_pass_query_range.end();
    }

    //!\brief The search algorithm called on every query in the underlying query_range.
    search_algorithm_t search_algorithm{};
    //!\brief The underlying query range wrapped in an single pass input range.
    single_pass_query_range_t single_pass_query_range{};
    //!\brief The current iterator over the single pass query range wrapper.
    single_pass_query_range_iterator current_query_it{};
    //!\brief Stores the current search results.
    result_buffer_t result_buffer{};
    //!\brief A flag that indicates whether begin was already called.
    bool first_call_of_begin{true};
};

/*!\brief The iterator of seqan3::detail::search_result_range.
 * \implements std::input_iterator
 */
template <typename search_algorithm_t, std::ranges::view query_range_t>
//!\cond
    requires std::ranges::forward_range<query_range_t> &&
             std::invocable<search_algorithm_t, std::ranges::range_reference_t<query_range_t>>
//!\endcond
class search_result_range<search_algorithm_t, query_range_t>::search_result_range_iterator
{
private:
    //!\brief The internal result buffer iterator.
    using result_buffer_iterator = std::ranges::iterator_t<result_buffer_t>;
public:
    /*!\name Associated types
     * \{
     */
    //!\brief Type for distances between iterators.
    using difference_type = std::ptrdiff_t;
    //!\brief Value type of container elements.
    using value_type = result_buffer_value_t;
    //!\brief Use reference type defined by container.
    using reference = result_buffer_value_t &;
    //!\brief Pointer type is pointer of container element type.
    using pointer = std::add_pointer_t<value_type>;
    //!\brief Sets iterator category as input iterator.
    using iterator_category = std::input_iterator_tag;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr search_result_range_iterator() noexcept = default; //!< Defaulted.
    constexpr search_result_range_iterator(search_result_range_iterator const &) noexcept = default; //!< Defaulted.
    constexpr search_result_range_iterator(search_result_range_iterator &&) noexcept = default; //!< Defaulted.
    constexpr search_result_range_iterator & operator=(search_result_range_iterator const &) noexcept = default; //!< Defaulted.
    constexpr search_result_range_iterator & operator=(search_result_range_iterator &&) noexcept = default; //!< Defaulted.
    ~search_result_range_iterator() = default; //!< Defaulted.

    /*!\brief Construct from the associated search result range.
     * \param[in] range The associated search result range.
     * \param[in] first_call_of_begin A flag indicating whether this is the first call of begin.
     *
     * \details
     *
     * Initialises the iterator with the associated search result range. If `first_call_of_begin` is `true`,
     * the search algorithm is invoked on the first query.
     */
    constexpr search_result_range_iterator(search_result_range & range, bool const first_call_of_begin) :
        range_ptr{& range},
        at_end(range.at_end())
    {
        if (first_call_of_begin)
            fetch_next_query_results();
    }
    //!\}

    /*!\name Access operators
     * \{
     */

    /*!\brief Access the pointed-to element.
     * \returns A reference to the current element.
     */
    //!\brief Returns the current search result.
    reference operator*() const noexcept
    {
        return *current_buffer_it;
    }

    //!\brief Returns a pointer to the current search result.
    pointer operator->() const noexcept
    {
        return std::addressof(*current_buffer_it);
    }
    //!\}

    /*!\name Increment operators
     * \{
     */
    //!\brief Increments the iterator by one.
    search_result_range_iterator & operator++(/*pre*/)
    {
        assert(range_ptr != nullptr);

        if (++current_buffer_it == last_buffer_it)
            fetch_next_query_results();

        return *this;
    }

    //!\brief Increments the iterator by one.
    void operator++(int /*post*/)
    {
        ++(*this);
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Checks whether `lhs` is equal to the sentinel.
    friend constexpr bool operator==(search_result_range_iterator const & lhs,
                                     std::ranges::default_sentinel_t const &) noexcept
    {
        return lhs.at_end && (lhs.current_buffer_it == lhs.last_buffer_it);
    }

    //!\brief Checks whether `lhs` is equal to `rhs`.
    friend constexpr bool operator==(std::ranges::default_sentinel_t const & lhs,
                                     search_result_range_iterator const & rhs) noexcept
    {
        return rhs == lhs;
    }

    //!\brief Checks whether `lhs` is not equal to the sentinel.
    friend constexpr bool operator!=(search_result_range_iterator const & lhs,
                                     std::ranges::default_sentinel_t const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Checks whether `lhs` is not equal to `rhs`.
    friend constexpr bool operator!=(std::ranges::default_sentinel_t const & lhs,
                                     search_result_range_iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}

private:
    /*!\brief Fetches the next query results and updates the corresponding iterators to the new search result buffer.
     *
     * \details
     *
     * If the last search did not generate any hits, the next query is searched. This process is repeated
     * until either a query produced at least one hit or all queries have been searched.
     */
    void fetch_next_query_results()
    {
        while (!(range_ptr->next() || range_ptr->at_end()));

        current_buffer_it = std::ranges::begin(range_ptr->result_buffer);
        last_buffer_it = std::ranges::end(range_ptr->result_buffer);

        at_end = range_ptr->at_end();
    }

    //!\brief The current iterator over the buffered search results.
    result_buffer_iterator current_buffer_it{};
    //!\brief The end iterator over the buffered search results.
    result_buffer_iterator last_buffer_it{};
    //!\brief Pointer to the underlying range.
    search_result_range * range_ptr{};
    //!\brief Indicates the end of the underlying range.
    bool at_end{true};
};

} // namespace seqan3
