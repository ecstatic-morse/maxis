// We need access to boost::dynamic_bitset's private fields to efficiently
// implement a hash function.
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS

#include "boost/dynamic_bitset.hpp"
#include "boost/functional/hash.hpp"
#include "boost/iterator/iterator_facade.hpp"

#ifndef MAXIS_BIT_VECTOR_H
#define MAXIS_BIT_VECTOR_H

// Boost insists on not defining itertors for dynamic_bitset. This will do for
// now
#define FOR_EACH_BV(v, bv)                                                     \
  for (auto v = (bv).find_first(); v != boost::dynamic_bitset<>::npos;         \
       v = (bv).find_next(v))

namespace maxis {

// The type used to store sets of vertices.
using BitVec = boost::dynamic_bitset<>;

namespace detail {

// Define an iterator for boost::dynamic_bitset which iterates over the indices
// of the set bits.
template <typename Container, typename Reference>
class bitset_iterator
    : public boost::iterator_facade<bitset_iterator<Container, Reference>, bool,
                                    boost::bidirectional_traversal_tag,
                                    Reference, ssize_t> {
public:
  bitset_iterator() : b{nullptr}, pos{0} {}
  bitset_iterator(Container *b, typename Container::size_type i = 0)
      : b{b}, pos{i} {}

  template <typename OtherContainer, typename OtherReference>
  bitset_iterator(const bitset_iterator<OtherContainer, OtherReference> &other)
      : b{other.b}, pos{other.pos} {}

private:
  friend class boost::iterator_core_access;

  Container *const b;
  typename Container::size_type pos;

  Reference dereference() const { return (*b)[pos]; }
  void increment() { ++pos; }
  void decrement() { --pos; }
  void advance(ssize_t n) { pos += n; }

  template <typename OtherContainer, typename OtherReference>
  bool
  equal(const bitset_iterator<OtherContainer, OtherReference> &other) const {
    return this->b == other.b && this->pos == other.pos;
  }

  template <typename OtherContainer, typename OtherReference>
  long int distance_to(
      const bitset_iterator<OtherContainer, OtherReference> &other) const {
    return this->pos - other.pos;
  }
};

} // namespace detail

using bitset_iterator =
    detail::bitset_iterator<boost::dynamic_bitset<>,
                            boost::dynamic_bitset<>::reference>;

using const_bitset_iterator =
    detail::bitset_iterator<const boost::dynamic_bitset<>,
                            boost::dynamic_bitset<>::const_reference>;

auto begin(boost::dynamic_bitset<> &b) -> bitset_iterator;
auto begin(const boost::dynamic_bitset<> &b) -> const_bitset_iterator;
auto end(boost::dynamic_bitset<> &b) -> bitset_iterator;
auto end(const boost::dynamic_bitset<> &b) -> const_bitset_iterator;
auto cbegin(const boost::dynamic_bitset<> &b) -> const_bitset_iterator;
auto cend(const boost::dynamic_bitset<> &b) -> const_bitset_iterator;

} // namespace maxis

namespace std {

// Define a hash function for boost::dynamic_bitset
template <typename B, typename A> struct hash<boost::dynamic_bitset<B, A>> {
  size_t operator()(const boost::dynamic_bitset<B, A> &s) const {
    std::size_t res = boost::hash_value(s.m_num_bits);
    boost::hash_combine(res, boost::hash_value(s.m_bits));
    return res;
  }
};

} // namespace std

#endif // MAXIS_BIT_VECTOR_H
