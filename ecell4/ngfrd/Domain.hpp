#ifndef ECELL4_NGFRD_DOMAIN_HPP
#define ECELL4_NGFRD_DOMAIN_HPP
#include <ecell4/ngfrd/SingleSphericalDomain.hpp>
#include <ecell4/ngfrd/SingleCircularDomain.hpp>
#include <ecell4/ngfrd/SingleConicalDomain.hpp>
#include <ecell4/ngfrd/MultiDomain.hpp>

namespace ecell4
{
namespace ngfrd
{

struct Domain
{
public:
    enum class DomainKind : int // boost::variant::which returns an int.
    {
        Uninitialized   = 0,
        SingleSpherical = 1,
        SingleCircular  = 2,
        SingleConical   = 3,
        Multi           = 4,
    };

    using storage_type = boost::variant<
            boost::blank,
            SingleSphericalDomain,
            SingleCircularDomain,
            SingleConicalDomain,
            MultiDomain
        >;

private:

    struct multiplicity_visitor : boost::static_visitor<std::size_t>
    {
        std::size_t operator()(const boost::blank&) const
        {
            throw_exception<IllegalState>("Domain is not initialized");
        }
        std::size_t operator()(const SingleSphericalDomain& dom) const noexcept
        {
            return dom.multiplicity();
        }
        std::size_t operator()(const SingleCircularDomain& dom) const noexcept
        {
            return dom.multiplicity();
        }
        std::size_t operator()(const SingleConicalDomain& dom) const noexcept
        {
            return dom.multiplicity();
        }
        std::size_t operator()(const MultiDomain& dom) const noexcept
        {
            return dom.multiplicity();
        }
    };

public:

    Domain() noexcept: storage_(boost::blank{}) {}

    template<typename D>
    explicit Domain(D&& d): storage_(std::forward<D>(d)) {}

    DomainKind kind() const noexcept {return DomainKind(storage_.which());}

    bool is_single_spherical() const noexcept
    {
        return this->kind() == DomainKind::SingleSpherical;
    }
    bool is_single_circular() const noexcept
    {
        return this->kind() == DomainKind::SingleCircular;
    }
    bool is_single_conical() const noexcept
    {
        return this->kind() == DomainKind::SingleConical;
    }
    bool is_multi() const noexcept
    {
        return this->kind() == DomainKind::Multi;
    }

    bool is_2D() const noexcept {return is_single_circular() || is_single_conical();}
    bool is_3D() const noexcept {return is_single_spherical();}

    SingleSphericalDomain const& as_single_spherical() const
    {
        return boost::get<SingleSphericalDomain>(storage_);
    }
    SingleSphericalDomain&       as_single_spherical()
    {
        return boost::get<SingleSphericalDomain>(storage_);
    }

    SingleCircularDomain const& as_single_circular() const
    {
        return boost::get<SingleCircularDomain>(storage_);
    }
    SingleCircularDomain&       as_single_circular()
    {
        return boost::get<SingleCircularDomain>(storage_);
    }

    SingleConicalDomain const& as_single_conical() const
    {
        return boost::get<SingleConicalDomain>(storage_);
    }
    SingleConicalDomain&       as_single_conical()
    {
        return boost::get<SingleConicalDomain>(storage_);
    }

    MultiDomain const& as_multi() const
    {
        return boost::get<MultiDomain>(storage_);
    }
    MultiDomain&       as_multi()
    {
        return boost::get<MultiDomain>(storage_);
    }

    storage_type const& as_variant() const noexcept {return storage_;}
    storage_type&       as_variant()       noexcept {return storage_;}

    std::size_t multiplicity() const
    {
        return boost::apply_visitor(multiplicity_visitor(), storage_);
    }

private:

    storage_type storage_;
};

} //ngfrd
} //ecell4
#endif// ECELL4_NGFRD_DOMAIN_HPP
