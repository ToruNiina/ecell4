#ifndef __ECELL4_ROD_HPP
#define __ECELL4_ROD_HPP

#include "Shape.hpp"

namespace ecell4
{

struct RodSurface;

// A Rod is parallel with the x-axis.
// The center of a Rod is the origin.
struct Rod
    : public Shape
{
public:

    Rod();
    Rod(const Real& length, const Real& radius);
    Rod(const Real& length, const Real& radius, const Real3& origin);
    Rod(const Rod& rhs);
    const Real& lengthX() const;
    const Real& radius() const;
    const Real3& origin() const;
    void shift(const Real3& vec);
    Real is_inside(const Real3& pos) const;
    Real distance(const Real3& pos) const;
    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const;
    RodSurface surface() const;

protected:

    Real length_; // LengthX
    Real radius_; // LengthY/2
    Real3 origin_;
};

struct RodSurface
    : public Shape
{
public:

    RodSurface();
    RodSurface(const Real& length, const Real& radius);
    RodSurface(const Real& length, const Real& radius, const Real3& origin);
    RodSurface(const RodSurface& rhs);
    const Real& lengthX() const;
    const Real& radius() const;
    const Real3& origin() const;
    void shift(const Real3& vec);
    Real is_inside(const Real3& pos) const;
    Real distance(const Real3& pos) const;
    Real3 draw_position(boost::shared_ptr<RandomNumberGenerator>& rng) const;
    Rod inside() const;

    dimension_kind dimension() const
    {
        return TWO;
    }

protected:

    Real length_; // LengthX
    Real radius_; // LengthY/2
    Real3 origin_; // origin
};

} // ecell4

#endif /* __ECELL4_ROD_HPP */
