// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::views::kmer_gapped_delete_mask_hash.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/math.hpp>
#include <seqan3/range/hash.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

namespace seqan3::detail
{
// ---------------------------------------------------------------------------------------------------------------------
// kmer_gapped_delete_mask_hash_view class
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief The type returned by seqan3::views::kmer_gapped_delete_mask_hash.
 * \tparam urng_t The type of the underlying ranges, must model std::forward_range, the reference type must model
 *                seqan3::semialphabet.
 * \implements std::ranges::view
 * \implements std::ranges::random_access_range
 * \implements std::ranges::sized_range
 * \ingroup views
 *
 * \details
 *
 * Note that most members of this class are generated by ranges::view_interface which is not yet documented here.
 */
template <std::ranges::view urng_t>
class kmer_gapped_delete_mask_hash_view : public std::ranges::view_interface<kmer_gapped_delete_mask_hash_view<urng_t>>
{
private:
    static_assert(std::ranges::forward_range<urng_t const>, "The kmer_gapped_delete_mask_hash_view only works on forward_ranges");
    static_assert(semialphabet<std::ranges::range_reference_t<urng_t>>,
                  "The reference type of the underlying range must model seqan3::semialphabet.");

    //!\brief The underlying range.
    urng_t urange;

    //!\brief The shape to use.
    shape shape_;

    template <typename rng_t>
    class shape_iterator;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    kmer_gapped_delete_mask_hash_view()                                       = default; //!< Defaulted.
    kmer_gapped_delete_mask_hash_view(kmer_gapped_delete_mask_hash_view const & rhs)             = default; //!< Defaulted.
    kmer_gapped_delete_mask_hash_view(kmer_gapped_delete_mask_hash_view && rhs)                  = default; //!< Defaulted.
    kmer_gapped_delete_mask_hash_view & operator=(kmer_gapped_delete_mask_hash_view const & rhs) = default; //!< Defaulted.
    kmer_gapped_delete_mask_hash_view & operator=(kmer_gapped_delete_mask_hash_view && rhs)      = default; //!< Defaulted.
    ~kmer_gapped_delete_mask_hash_view()                                      = default; //!< Defaulted.

    /*!\brief Construct from a view and a given shape.
     * \throws std::invalid_argument if hashes resulting from the shape/alphabet combination cannot be represented in
     *         `uint64_t`, i.e. \f$s>\frac{64}{\log_2\sigma}\f$ with shape size \f$s\f$ and alphabet size \f$\sigma\f$.
     */
    kmer_gapped_delete_mask_hash_view(urng_t urange_, shape const & s_) : urange{std::move(urange_)}, shape_{s_}
    {
        if (shape_.count() > (64 / std::log2(alphabet_size<std::ranges::range_reference_t<urng_t>>)))
        {
            throw std::invalid_argument{"The chosen shape/alphabet combination is not valid. "
                                        "The alphabet or shape size must be reduced."};
        }
    }

    /*!\brief Construct from a non-view that can be view-wrapped and a given shape.
     * \throws std::invalid_argument if hashes resulting from the shape/alphabet combination cannot be represented in
     *         `uint64_t`, i.e. \f$s>\frac{64}{\log_2\sigma}\f$ with shape size \f$s\f$ and alphabet size \f$\sigma\f$.
     */
    template <typename rng_t>
    //!\cond
     requires !std::same_as<remove_cvref_t<rng_t>, kmer_gapped_delete_mask_hash_view> &&
              std::ranges::viewable_range<rng_t> &&
              std::constructible_from<urng_t, ranges::ref_view<std::remove_reference_t<rng_t>>>
    //!\endcond
    kmer_gapped_delete_mask_hash_view(rng_t && urange_, shape const & s_) :
        urange{std::views::all(std::forward<rng_t>(urange_))}, shape_{s_}
    {
        if (shape_.count() > (64 / std::log2(alphabet_size<std::ranges::range_reference_t<urng_t>>)))
        {
            throw std::invalid_argument{"The chosen shape/alphabet combination is not valid. "
                                        "The alphabet or shape size must be reduced."};
        }
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the range.
     * \returns Iterator to the first element.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    auto begin() noexcept
    {
        return shape_iterator<urng_t>{std::ranges::begin(urange), std::ranges::end(urange), shape_};
    }

    //!\copydoc begin()
    auto begin() const noexcept
    //!\cond
        requires const_iterable_range<urng_t>
    //!\endcond
    {
        return shape_iterator<urng_t const>{std::ranges::begin(urange), std::ranges::end(urange), shape_};
    }

    //!\copydoc begin()
    auto cbegin() const noexcept
    //!\cond
        requires const_iterable_range<urng_t>
    //!\endcond
    {
        return begin();
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator to the end.
     *
     * \details
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    auto end() noexcept
    {
        return std::ranges::end(urange);
    }

    //!\copydoc end()
    auto end() const noexcept
    //!\cond
        requires const_iterable_range<urng_t>
    //!\endcond
    {
        return std::ranges::end(urange);
    }

    //!\copydoc end()
    auto cend() const noexcept
    //!\cond
        requires const_iterable_range<urng_t>
    //!\endcond
    {
        return end();
    }
    //!\}

    /*!\brief Returns the size of the range, if the underlying range is a std::ranges::sized_range.
     * \returns Size of range.
     */
    auto size() const
    //!\cond
        requires std::ranges::sized_range<urng_t>
    //!\endcond
    {
        using size_type = decltype(std::ranges::size(urange));
        return std::max<size_type>(std::ranges::size(urange) + 1, shape_.size()) - shape_.size();
    }
};

/*!\brief Iterator for calculating hash values via a given seqan3::shape.
 * \tparam urng_t Type of the text. Must model std::forward_range. Reference type must model seqan3::semialphabet.
 *
 * \details
 *
 * The shape_iterator can be used to iterate over the hash values of a text. A shape_iterator needs an iterator of
 * the text and a seqan3::shape that defines how to hash the text.
 *
 * Depending on the type of the iterator passed to the shape_iterator, different functionality is available:
 *
 * | Concept modelled by passed text iterator | Available functions             |
 * |------------------------------------------|---------------------------------|
 * | std::forward_iterator                    | \ref shape_iterator_comparison "Comparison operators"<br>\ref operator++ "Pre-increment (++it)"<br>\ref operator++(int) "Post-increment (it++)"<br>\ref operator* "Indirection operator (*it)" |
 * | std::bidirectional_iterator              | \ref operator-- "Pre-decrement (--it)"<br>\ref operator--(int) "Post-decrement (it--)" |
 * | std::random_access_iterator              | \ref operator+= "Forward (it +=)"<br>\ref operator+ "Forward copy (it +)"<br>\ref operator-= "Decrement(it -=)"<br>\ref shape_iterator_operator-decrement "Decrement copy (it -)"<br>\ref shape_iterator_operator-difference "Difference (it1 - it2)"<br>\ref operator[] "Subscript (it[])" |
 *
 * When using a gapped seqan3::shape, the `0`s of the seqan3::shape are virtually removed from the hashed k-mer.
 * Note that any shape is expected to start with a `1` and end with a `1`.
 *
 * ### Implementation detail
 *
 * To avoid dereferencing the sentinel when iterating, the shape_iterator computes the hash value up until
 * the second to last position and performs the addition of the last position upon
 * access (\ref operator* and \ref operator[]).
 */
template <std::ranges::view urng_t>
template <typename rng_t>
class kmer_gapped_delete_mask_hash_view<urng_t>::shape_iterator
{
private:
    //!\brief The iterator type of the underlying range.
    using it_t = std::ranges::iterator_t<rng_t>;
    //!\brief The sentinel type of the underlying range.
    using sentinel_t = std::ranges::sentinel_t<rng_t>;

    template <typename urng2_t>
    friend class shape_iterator;

public:
    /*!\name Associated types
     * \{
     */
    //!\brief Type for distances between iterators.
    using difference_type = typename std::iter_difference_t<it_t>;
    //!\brief Value type of this iterator.
    using value_type = size_t;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief Reference to `value_type`.
    using reference = value_type;
    //!\brief Tag this class as input iterator.
    using iterator_category = typename std::iterator_traits<it_t>::iterator_category;
    //!\brief Tag this class depending on which concept `it_t` models.
    using iterator_concept = std::conditional_t<std::contiguous_iterator<it_t>,
                                                 typename std::random_access_iterator_tag,
                                                 iterator_tag_t<it_t>>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr shape_iterator()                                   = default; //!< Defaulted.
    constexpr shape_iterator(shape_iterator const &)             = default; //!< Defaulted.
    constexpr shape_iterator(shape_iterator &&)                  = default; //!< Defaulted.
    constexpr shape_iterator & operator=(shape_iterator const &) = default; //!< Defaulted.
    constexpr shape_iterator & operator=(shape_iterator &&)      = default; //!< Defaulted.
    ~shape_iterator()                                            = default; //!< Defaulted.

    //!\brief Allow iterator on a const range to be constructible from an iterator over a non-const range.
    template <typename urng2_t>
    //!\cond
        requires std::same_as<std::remove_const_t<urng_t>, urng2_t>
    //!\endcond
    shape_iterator(shape_iterator<urng2_t> it) :
        hash_value{std::move(it.hash_value)},
        roll_factor{std::move(it.roll_factor)},
        shape_{std::move(it.shape_)},
        text_left{std::move(it.text_left)},
        text_right{std::move(it.text_right)}
    {}

    /*!\brief Construct from a given iterator on the text and a seqan3::shape.
    * /param[in] it_start Iterator pointing to the first position of the text.
    * /param[in] s_       The seqan3::shape that determines which positions participate in hashing.
    *
    * \details
    *
    * ### Complexity
    *
    * Linear in size of shape.
    */
    shape_iterator(it_t it_start, sentinel_t it_end, shape s_) :
        shape_{s_}, text_left{it_start}, text_right{std::ranges::next(it_start, s_.size(), it_end)}
    {
        assert(std::ranges::size(shape_) > 0);

        if (shape_.size() <= std::ranges::distance(text_left, text_right))
        {
            roll_factor = pow(sigma, static_cast<size_t>(std::ranges::size(shape_) - 1));

            for (size_t i{0}; i < shape_.size() - 1u; ++i)
            {
                delete_mask <<= 2;
                delete_mask += 3;
            }

            if(!shape_.all()) 
            {
                for (size_t i{0}; i < shape_.size() - 1u; ++i)
                {
                    shape_mask += shape_[i];
                    shape_mask *= 2;
                    shape_mask += shape_[i];
                    shape_mask *= 2;
                }
                shape_mask += shape_[shape_.size() - 1u];
                shape_mask *= 2;
                shape_mask += shape_[shape_.size() - 1u];
            }

            hash_full();
        }
    }
    //!\}

    //!\anchor shape_iterator_comparison
    //!\name Comparison operators
    //!\{

    //!\brief Compare to iterator on text.
    friend bool operator==(shape_iterator const & lhs, sentinel_t const & rhs) noexcept
    {
        return lhs.text_right == rhs;
    }

    //!\brief Compare to iterator on text.
    friend bool operator==(sentinel_t const & lhs, shape_iterator const & rhs) noexcept
    {
        return lhs == rhs.text_right;
    }

    //!\brief Compare to another shape_iterator.
    friend bool operator==(shape_iterator const & lhs, shape_iterator const & rhs) noexcept
    {
        return std::tie(lhs.text_right, lhs.shape_) == std::tie(rhs.text_right, rhs.shape_);
    }

    //!\brief Compare to iterator on text.
    friend bool operator!=(shape_iterator const & lhs, sentinel_t const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to iterator on text.
    friend bool operator!=(sentinel_t const & lhs, shape_iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to another shape_iterator.
    friend bool operator!=(shape_iterator const & lhs, shape_iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to another shape_iterator.
    friend bool operator<(shape_iterator const & lhs, shape_iterator const & rhs) noexcept
    {
        return (lhs.shape_ <= rhs.shape_) && (lhs.text_right < rhs.text_right);
    }

    //!\brief Compare to another shape_iterator.
    friend bool operator>(shape_iterator const & lhs, shape_iterator const & rhs) noexcept
    {
        return (lhs.shape_ >= rhs.shape_) && (lhs.text_right > rhs.text_right);
    }

    //!\brief Compare to another shape_iterator.
    friend bool operator<=(shape_iterator const & lhs, shape_iterator const & rhs) noexcept
    {
        return (lhs.shape_ <= rhs.shape_) && (lhs.text_right <= rhs.text_right);
    }

    //!\brief Compare to another shape_iterator.
    friend bool operator>=(shape_iterator const & lhs, shape_iterator const & rhs) noexcept
    {
        return (lhs.shape_ >= rhs.shape_) && (lhs.text_right >= rhs.text_right);
    }

    //!\}

    //!\brief Pre-increment.
    shape_iterator & operator++() noexcept
    {
        hash_forward();
        return *this;
    }

    //!\brief Post-increment.
    shape_iterator operator++(int) noexcept
    {
        shape_iterator tmp{*this};
        hash_forward();
        return tmp;
    }

    /*!\brief Pre-decrement.
     * \attention This function is only avaible if `it_t` models std::bidirectional_iterator.
     */
    shape_iterator & operator--() noexcept
    //!\cond
        requires std::bidirectional_iterator<it_t>
    //!\endcond
    {
        hash_backward();
        return *this;
    }

    /*!\brief Post-decrement.
     * \attention This function is only avaible if `it_t` models std::bidirectional_iterator.
     */
    shape_iterator operator--(int) noexcept
    //!\cond
        requires std::bidirectional_iterator<it_t>
    //!\endcond
    {
        shape_iterator tmp{*this};
        hash_backward();
        return tmp;
    }

    /*!\brief Forward this iterator.
     * \attention This function is only avaible if `it_t` models std::random_access_iterator.
     */
    shape_iterator & operator+=(difference_type const skip) noexcept
    //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        hash_forward(skip);
        return *this;
    }

    /*!\brief Forward copy of this iterator.
     * \attention This function is only avaible if `it_t` models std::random_access_iterator.
     */
    shape_iterator operator+(difference_type const skip) const noexcept
    //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        shape_iterator tmp{*this};
        return tmp += skip;
    }

    /*!\brief Non-member operator+ delegates to non-friend operator+.
     * \attention This function is only avaible if `it_t` models std::random_access_iterator.
     */
    friend shape_iterator operator+(difference_type const skip, shape_iterator const & it) noexcept
    //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        return it + skip;
    }

    /*!\brief Decrement iterator by `skip`.
     * \attention This function is only avaible if `it_t` models std::random_access_iterator.
     */
    shape_iterator & operator-=(difference_type const skip) noexcept
    //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        hash_backward(skip);
        return *this;
    }

    /*!\anchor shape_iterator_operator-decrement
     * \brief Return decremented copy of this iterator.
     * \attention This function is only avaible if `it_t` models std::random_access_iterator.
     */
    shape_iterator operator-(difference_type const skip) const noexcept
    //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        shape_iterator tmp{*this};
        return tmp -= skip;
    }

    /*!\brief Non-member operator- delegates to non-friend operator-.
     * \attention This function is only avaible if `it_t` models std::random_access_iterator.
     */
    friend shape_iterator operator-(difference_type const skip, shape_iterator const & it) noexcept
    //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        return it - skip;
    }

    /*!\anchor shape_iterator_operator-difference
     * \brief Return offset between this and remote iterator's position.
     * \attention This function is only avaible if `it_t` models std::random_access_iterator.
     */
    difference_type operator-(shape_iterator const & lhs) const noexcept
    //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        return static_cast<difference_type>(text_right - lhs.text_right);
    }

    /*!\brief Return offset between remote sentinel's position and this.
     * \attention This function is only avaible if sentinel_t and it_t model std::sized_sentinel_for.
     */
    friend difference_type operator-(sentinel_t const & lhs, shape_iterator const & rhs) noexcept
    //!\cond
        requires std::sized_sentinel_for<sentinel_t, it_t>
    //!\endcond
    {
        return static_cast<difference_type>(lhs - rhs.text_right);
    }

    /*!\brief Return offset this and remote sentinel's position.
     * \attention This function is only avaible if it_t and sentinel_t model std::sized_sentinel_for.
     */
    friend difference_type operator-(shape_iterator const & lhs, sentinel_t const & rhs) noexcept
    //!\cond
        requires std::sized_sentinel_for<it_t, sentinel_t>
    //!\endcond
    {
        return static_cast<difference_type>(lhs.text_right - rhs);
    }

    /*!\brief Move the iterator by a given offset and return the corresponding hash value.
     * \attention This function is only avaible if `it_t` models std::random_access_iterator.
     */
    reference operator[](difference_type const n) const
    //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        return *(*this + n);
    }

    //!\brief Return the hash value.
    value_type operator*() const noexcept
    {
        if(shape_.all()) 
        {
            return hash_value + to_rank(*text_right);
        }
        return (hash_value + to_rank(*text_right)) & shape_mask;
    }

private:
    //!\brief The alphabet type of the passed iterator.
    using alphabet_t = std::iter_value_t<it_t>;

    //!\brief The alphabet size.
    static constexpr auto const sigma{alphabet_size<alphabet_t>};

    //!\brief The hash value.
    size_t hash_value{0};

    //!\brief The factor for the left most position of the hash value.
    size_t roll_factor{0};

    //!\brief The shape to use.
    shape shape_;

    size_t shape_mask{0};

    size_t delete_mask{0};

    //!\brief Iterator to the leftmost position of the k-mer.
    it_t text_left;

    //!\brief Iterator to the rightmost position of the k-mer.
    it_t text_right;

    //!\brief Increments iterator by 1.
    void hash_forward()
    {
        hash_roll_forward();
    }

    /*!\brief Increments iterator by `skip`.
     * \param skip Amount to increment.
     * \attention This function is only avaible if `it_t` models std::random_access_iterator.
     */
    void hash_forward(difference_type const skip)
    //!\cond
        requires std::random_access_iterator<it_t>
    //!\endcond
    {
        std::ranges::advance(text_left, skip);
        hash_full();
    }

    /*!\brief Decrements iterator by 1.
     * \attention This function is only avaible if `it_t` models std::bidirectional_iterator.
     */
    void hash_backward()
    //!\cond
        requires std::bidirectional_iterator<it_t>
    //!\endcond
    {
        if (shape_.all())
        {
            hash_roll_backward();
        }
        else
        {
            std::ranges::advance(text_left,  -1);
            hash_full();
        }
    }

    /*!\brief Decrements iterator by `skip`.
     * \param skip Amount to decrement.
     * \attention This function is only avaible if `it_t` models std::bidirectional_iterator.
     */
    void hash_backward(difference_type const skip)
    {
        std::ranges::advance(text_left, -skip);
        hash_full();
    }

    //!\brief Calculates a hash value by explicitly looking at each position.
    void hash_full()
    {
        text_right = text_left;
        hash_value = 0;

        for (size_t i{0}; i < shape_.size() - 1u; ++i)
        {
            hash_value += shape_[i] * to_rank(*text_right);
            hash_value *= shape_[i] ? sigma : 1;
            std::ranges::advance(text_right, 1);
        }
    }

    //!\brief Calculates the next hash value via rolling hash.
    void hash_roll_forward()
    {
        hash_value &= delete_mask;
        hash_value += to_rank(*(text_right));
        hash_value *= sigma;

        std::ranges::advance(text_left,  1);
        std::ranges::advance(text_right, 1);
    }

    /*!\brief Calculates the previous hash value via rolling hash.
     * \attention This function is only avaible if `it_t` models std::bidirectional_iterator.
     */
    void hash_roll_backward()
        //!\cond
        requires std::bidirectional_iterator<it_t>
        //!\endcond
    {
        std::ranges::advance(text_left,  -1);
        std::ranges::advance(text_right, -1);

        hash_value /= sigma;
        hash_value -= to_rank(*(text_right));
        hash_value += to_rank(*(text_left)) * roll_factor;
    }
};

//!\brief A deduction guide for the view class template.
template <std::ranges::viewable_range rng_t>
kmer_gapped_delete_mask_hash_view(rng_t &&, shape const & shape_) -> kmer_gapped_delete_mask_hash_view<std::ranges::all_view<rng_t>>;

// ---------------------------------------------------------------------------------------------------------------------
// kmer_gapped_delete_mask_hash_fn (adaptor definition)
// ---------------------------------------------------------------------------------------------------------------------

//![adaptor_def]
//!\brief views::kmer_gapped_delete_mask_hash's range adaptor object type (non-closure).
struct kmer_gapped_delete_mask_hash_fn
{
    //!\brief Store the shape and return a range adaptor closure object.
    constexpr auto operator()(shape const & shape_) const
    {
        return adaptor_from_functor{*this, shape_};
    }

    /*!\brief            Call the view's constructor with the underlying view and a seqan3::shape as argument.
     * \param[in] urange The input range to process. Must model std::ranges::viewable_range and the reference type
     *                   of the range must model seqan3::semialphabet.
     * \param[in] shape_ The seqan3::shape to use for hashing.
     * \throws std::invalid_argument if resulting hash values would be too big for a 64 bit integer.
     * \returns          A range of converted elements.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, shape const & shape_) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::kmer_gapped_delete_mask_hash cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::kmer_gapped_delete_mask_hash must model std::ranges::forward_range.");
        static_assert(semialphabet<std::ranges::range_reference_t<urng_t>>,
            "The range parameter to views::kmer_gapped_delete_mask_hash must be over elements of seqan3::semialphabet.");

        return kmer_gapped_delete_mask_hash_view{std::forward<urng_t>(urange), shape_};
    }
};
//![adaptor_def]

} // namespace seqan3::detail

namespace seqan3::views
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief               Computes hash values for each position of a range via a given shape.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] shape     The seqan3::shape that determines how to compute the hash value.
 * \returns             A range of std::size_t where each value is the hash of the resp. k-mer.
 *                      See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * \attention
 * For the alphabet size \f$\sigma\f$ of the alphabet of `urange` and the number of 1s \f$s\f$ of `shape` it must hold
 * that \f$s>\frac{64}{\log_2\sigma}\f$, i.e. hashes resulting from the shape/alphabet combination can be represented
 * in an `uint64_t`.
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)   | `rrng_t` (returned range type)   |
 * |----------------------------------|:----------------------------------:|:--------------------------------:|
 * | std::ranges::input_range         | *required*                         | *preserved*                      |
 * | std::ranges::forward_range       | *required*                         | *preserved*                      |
 * | std::ranges::bidirectional_range |                                    | *preserved*                      |
 * | std::ranges::random_access_range |                                    | *preserved*                      |
 * | std::ranges::contiguous_range    |                                    | *lost*                           |
 * |                                  |                                    |                                  |
 * | std::ranges::viewable_range      | *required*                         | *guaranteed*                     |
 * | std::ranges::view                |                                    | *guaranteed*                     |
 * | std::ranges::sized_range         |                                    | *preserved*                      |
 * | std::ranges::common_range        |                                    | *lost*                           |
 * | std::ranges::output_range        |                                    | *lost*                           |
 * | seqan3::const_iterable_range     |                                    | *preserved*                      |
 * |                                  |                                    |                                  |
 * | std::ranges::range_reference_t   | seqan3::semialphabet               | std::size_t                      |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \include test/snippet/range/views/kmer_gapped_delete_mask_hash.cpp
 *
 * \hideinitializer
 */
inline constexpr auto kmer_gapped_delete_mask_hash = detail::kmer_gapped_delete_mask_hash_fn{};

//!\}

} // namespace seqan3::views