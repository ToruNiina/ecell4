#ifndef ECELL4_REAL3_HPP
#define ECELL4_REAL3_HPP

#include <ostream>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <cmath>
#include <boost/array.hpp>

#include <ecell4/core/config.h>

#include "types.hpp"
#include "functions.hpp"

#include "hash.hpp"

namespace ecell4
{

struct Real3
{
    typedef boost::array<Real, 3> base_type;
    typedef base_type::value_type             value_type;
    typedef base_type::size_type              size_type;
    typedef base_type::iterator               iterator;
    typedef base_type::const_iterator         const_iterator;
    typedef base_type::reverse_iterator       reverse_iterator;
    typedef base_type::const_reverse_iterator const_reverse_iterator;
    typedef base_type::reference              reference;
    typedef base_type::const_reference        const_reference;
    typedef base_type::difference_type        difference_type;

    Real3& operator+=(const Real3& rhs);
    Real3& operator-=(const Real3& rhs);
    Real3& operator*=(const Real3::value_type& rhs);
    Real3& operator/=(const Real3::value_type& rhs);

    Real3()
    {
        this->base_[0] = 0;
        this->base_[1] = 0;
        this->base_[2] = 0;
    }

    Real3(value_type p0, value_type p1, value_type p2)
    {
        this->base_[0] = p0;
        this->base_[1] = p1;
        this->base_[2] = p2;
    }

    Real3(const Real3 &rhs)
    {
        this->base_[0] = rhs[0];
        this->base_[1] = rhs[1];
        this->base_[2] = rhs[2];
    }

    // Real3(const Real (&a)[3])
    //     : base_type(*reinterpret_cast<const base_type*>(&a))
    // {
    //     ;
    // }

    // Real3(const Real a[3])
    //     : base_type(*reinterpret_cast<const base_type*>(a))
    // {
    //     ;
    // }

    // Real3(const base_type& a)
    //     : base_type(a)
    // {
    //     ;
    // }

    Real3& operator=(const Real3& other)
    {
        this->base_ = other.base_;
        return *this;
    }
    template<typename T>
    Real3& operator=(const boost::array<T, 3>& other)
    {
        this->base_ = other;
        return *this;
    }

    size_type     size() const throw() {return base_.size();}
    bool         empty() const throw() {return base_.empty();}
    size_type max_size() const throw() {return base_.max_size();}

    Real& operator[](size_type i)       throw() {return base_[i];}
    Real  operator[](size_type i) const throw() {return base_[i];}
    Real&         at(size_type i)       throw() {return base_.at(i);}
    Real          at(size_type i) const throw() {return base_.at(i);}

    Real& front()       throw() {return base_.front();}
    Real  front() const throw() {return base_.front();}
    Real& back()        throw() {return base_.back();}
    Real  back()  const throw() {return base_.back();}

    iterator       begin()        throw() {return base_.begin();}
    iterator       end()          throw() {return base_.end();}
    const_iterator begin()  const throw() {return base_.begin();}
    const_iterator end()    const throw() {return base_.end();}
    const_iterator cbegin() const throw() {return base_.cbegin();}
    const_iterator cend()   const throw() {return base_.cend();}

    reverse_iterator       rbegin()        throw() {return base_.rbegin();}
    reverse_iterator       rend()          throw() {return base_.rend();}
    const_reverse_iterator rbegin()  const throw() {return base_.rbegin();}
    const_reverse_iterator rend()    const throw() {return base_.rend();}
    const_reverse_iterator crbegin() const throw() {return base_.crbegin();}
    const_reverse_iterator crend()   const throw() {return base_.crend();}

    const Real* data() const {return base_.data();}
    Real*    c_array()       {return base_.c_array();}

    void swap  (Real3& other) {this->base_.swap(other.base_);}
    void assign(const Real& v){this->base_.fill(v);}

    base_type&       base()       throw() {return base_;}
    base_type const& base() const throw() {return base_;}

  private:
    base_type base_;
};

// derived from boost::array

inline void swap(Real3& lhs, Real3& rhs) {lhs.swap(rhs);}
template<std::size_t N>
inline Real  get(const Real3& v) {return boost::get<N>(v.base());}
template<std::size_t N>
inline Real& get(Real3& v)       {return boost::get<N>(v.base());}

inline bool operator==(const Real3& lhs, const Real3& rhs)
{
    return lhs.base() == rhs.base();
}
inline bool operator!=(const Real3& lhs, const Real3& rhs)
{
    return lhs.base() != rhs.base();
}
inline bool operator< (const Real3& lhs, const Real3& rhs)
{
    return lhs.base() <  rhs.base();
}
inline bool operator<=(const Real3& lhs, const Real3& rhs)
{
    return lhs.base() <= rhs.base();
}
inline bool operator> (const Real3& lhs, const Real3& rhs)
{
    return lhs.base() >  rhs.base();
}
inline bool operator>=(const Real3& lhs, const Real3& rhs)
{
    return lhs.base() >= rhs.base();
}

inline Real3 add(const Real3& p1, const Real3& p2)
{
    Real3 retval;
    retval[0] = p1[0] + p2[0];
    retval[1] = p1[1] + p2[1];
    retval[2] = p1[2] + p2[2];
    return retval;
}

inline Real3 subtract(const Real3& p1, const Real3& p2)
{
    Real3 retval;
    retval[0] = p1[0] - p2[0];
    retval[1] = p1[1] - p2[1];
    retval[2] = p1[2] - p2[2];
    return retval;
}

inline Real3 divide(const Real3& p1, const Real3::value_type& p2)
{
    Real3 retval;
    retval[0] = p1[0] / p2;
    retval[1] = p1[1] / p2;
    retval[2] = p1[2] / p2;
    return retval;
}

inline Real3 multiply(const Real3& p1, const Real3::value_type& p2)
{
    Real3 retval;
    retval[0] = p1[0] * p2;
    retval[1] = p1[1] * p2;
    retval[2] = p1[2] * p2;
    return retval;
}

inline Real3 modulo(const Real3& p1, const Real3::value_type& p2)
{
    Real3 retval;
    retval[0] = modulo(p1[0], p2);
    retval[1] = modulo(p1[1], p2);
    retval[2] = modulo(p1[2], p2);
    return retval;
}

inline Real3 modulo(const Real3& p1, const Real3& p2)
{
    Real3 retval;
    retval[0] = modulo(p1[0], p2[0]);
    retval[1] = modulo(p1[1], p2[1]);
    retval[2] = modulo(p1[2], p2[2]);
    return retval;
}

inline Real3 abs(const Real3& v)
{
    Real3 retval;
    retval[0] = abs(v[0]);
    retval[1] = abs(v[1]);
    retval[2] = abs(v[2]);
    return retval;
}

inline Real3::value_type dot_product(
    const Real3& p1, const Real3& p2)
{
    return p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
}

inline Real3 cross_product(const Real3& p1, const Real3& p2)
{
    Real3 retval;
    retval[0] = p1[1] * p2[2] - p1[2] * p2[1];
    retval[1] = p1[2] * p2[0] - p1[0] * p2[2];
    retval[2] = p1[0] * p2[1] - p1[1] * p2[0];
    return retval;
}

inline Real3::value_type length_sq(const Real3& r)
{
    return pow_2(r[0]) + pow_2(r[1]) + pow_2(r[2]);
}

inline Real3::value_type length(const Real3& r)
{
    return std::sqrt(length_sq(r));
}

inline Real3 operator+(const Real3& lhs, const Real3& rhs)
{
    return add(lhs, rhs);
}

inline Real3 operator-(const Real3& lhs, const Real3& rhs)
{
    return subtract(lhs, rhs);
}

inline Real3 operator/(
    const Real3& lhs, const Real3::value_type& rhs)
{
    return divide(lhs, rhs);
}

inline Real3 operator*(
    const Real3& lhs, const Real3::value_type& rhs)
{
    return multiply(lhs, rhs);
}

inline Real3& Real3::operator+=(const Real3& rhs)
{
    *this = add(*this, rhs);
    return *this;
}

inline Real3& Real3::operator-=(const Real3& rhs)
{
    *this = subtract(*this, rhs);
    return *this;
}

inline Real3& Real3::operator*=(const Real3::value_type& rhs)
{
    *this = multiply(*this, rhs);
    return *this;
}

inline Real3& Real3::operator/=(const Real3::value_type& rhs)
{
    *this = divide(*this, rhs);
    return *this;
}

template<typename Tstrm_, typename Ttraits_>
inline std::basic_ostream<Tstrm_, Ttraits_>& operator<<(
    std::basic_ostream<Tstrm_, Ttraits_>& strm, const Real3& v)
{
    strm << std::setprecision(12)
         << "(" << v[0] <<  ", " << v[1] <<  ", " << v[2] << ")";
    return strm;
}

inline Real3 ones()
{
    return Real3(1.0, 1.0, 1.0);
}

inline Real3 unitx()
{
    return Real3(1.0, 0.0, 0.0);
}

inline Real3 unity()
{
    return Real3(0.0, 1.0, 0.0);
}

inline Real3 unitz()
{
    return Real3(0.0, 0.0, 1.0);
}

} // ecell4

ECELL4_DEFINE_HASH_BEGIN()

template<>
struct hash<ecell4::Real3>
{
    typedef ecell4::Real3 argument_type;

    std::size_t operator()(argument_type const& val)
    {
        return hash<argument_type::value_type>()(val[0]) ^
            hash<argument_type::value_type>()(val[1]) ^
            hash<argument_type::value_type>()(val[2]);
    }
};

ECELL4_DEFINE_HASH_END()

#endif /* ECELL4_REAL3_HPP */
