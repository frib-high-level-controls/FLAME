#ifndef UTIL_H
#define UTIL_H

#include <stdexcept>
#include <iterator>

struct key_error : public std::runtime_error
{
    key_error(const std::string& s) : std::runtime_error(s) {}
};

//! An iterator which adapts the second member of a std::pair
//! Allows iteration over a values of a std::map<K,V> to look like
//! iteration over a std::list<V> or std::vector<V>.
template<typename I>
class value_proxy_iterator
{
    I orig;
public:
    typedef typename I::iterator_category iterator_category;
    typedef typename I::value_type::second_type value_type;
    typedef std::ptrdiff_t                   difference_type;
    typedef value_type*                  pointer;
    typedef value_type&                  reference;

    value_proxy_iterator() {}
    explicit value_proxy_iterator(const I& o) :orig(o) {}
    value_proxy_iterator(const value_proxy_iterator& o) :orig(o.orig) {}

    reference operator*() const { return orig->second; }
    pointer operator->() const { return &orig->second; }
    reference operator[](size_t i) const { return orig[i]->second; }

    bool operator==(const value_proxy_iterator& o) const { return orig==o.orig; }
    bool operator!=(const value_proxy_iterator& o) const { return orig!=o.orig; }

    value_proxy_iterator& operator++() { ++orig; return *this; }
    value_proxy_iterator  operator++(int) { return value_proxy_iterator(orig++); }
    value_proxy_iterator& operator--() { --orig; return *this; }
    value_proxy_iterator  operator--(int) { return value_proxy_iterator(orig--); }

    value_proxy_iterator& operator+=(size_t i) { orig+=i; return *this; }
    value_proxy_iterator& operator-=(size_t i) { orig-=i; return *this; }
    value_proxy_iterator  operator+(size_t i) const { return value_proxy_iterator(orig+i); }
    value_proxy_iterator  operator-(size_t i) const { return value_proxy_iterator(orig-i); }
};

#endif // UTIL_H
