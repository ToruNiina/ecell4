#ifndef ECELL4_NGFRD_BD_MATH_HPP
#define ECELL4_NGFRD_BD_MATH_HPP
#include <ecell4/core/types.hpp>

namespace ecell4
{
namespace ngfrd
{
// determine the distance between particles that are dissociated from each other
namespace bd_math
{
Real drawR_gbd_3D(const Real sigma, const Real t, const Real D, const Real rnd) noexcept;
} // bd_math

} // ngfrd
} // ecell4
#endif /* ECELL4_NGFRD_BD_MATH_HPP */
