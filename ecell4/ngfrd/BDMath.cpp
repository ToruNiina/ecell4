#include <ecell4/core/exceptions.hpp>
#include <ecell4/ngfrd/BDMath.hpp>
#include <greens_functions/freeFunctions.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/tools/roots.hpp>
#include <ciso646>

namespace ecell4
{
namespace ngfrd
{
// determine the distance between particles that are dissociated from each other
namespace bd_math
{
Real Igbd_3d(const Real sigma, const Real t, const Real D) noexcept
{
    constexpr Real sqrtPi = boost::math::constants::root_pi<Real>();

    const Real Dt      = D * t;
    const Real Dt2     = Dt + Dt;
    const Real sqrtDt  = std::sqrt(Dt);
    const Real sigmasq = sigma * sigma;

    constexpr Real term1 = 1 / (3 * sqrtPi);
    const     Real term2 = sigmasq - Dt2;
    const     Real term3 = Dt2 - 3 * sigmasq;
    const     Real term4 = sqrtPi * sigmasq * sigma *
                           boost::math::erfc<Real>(sigma / sqrtDt);

    return term1 * (-sqrtDt * (term2 * std::exp(-sigmasq / Dt) + term3) + term4);
}

Real Igbd_r_3d(const Real r, const Real sigma, const Real t, Real D) noexcept
{
    constexpr Real one_div_root_pi = boost::math::constants::one_div_root_pi<Real>();

    const Real Dt  = D * t;
    const Real Dt2 = Dt + Dt;
    const Real Dt4 = Dt2 + Dt2;
    const Real sqrtDt  = std::sqrt(Dt);
    const Real sqrtDt4 = 2 * sqrtDt;
    const Real sigmasq = sigma * sigma;
    const Real sigmacb = sigmasq * sigma;
    const Real rcb     = r * r * r;

    const Real rsigma =  r * sigma;
    const Real rps_sq = (r + sigma) * (r + sigma);
    const Real rms_sq = (r - sigma) * (r - sigma);

    const Real term1 = -2 * sqrtDt * one_div_root_pi;
    const Real term2 =  std::exp(-sigmasq / Dt) * (sigmasq - Dt2);
    const Real term3 = -std::exp(-rps_sq / Dt4) * (rms_sq + rsigma - Dt2);
    const Real term4 =  std::exp(-rms_sq / Dt4) * (rps_sq - rsigma - Dt2);
    const Real term5 = -sigmasq * 3 + Dt2;

    const Real term6 =  (sigmacb - rcb)     * boost::math::erf((r - sigma) / sqrtDt4);
    const Real term7 = -(sigmacb + sigmacb) * boost::math::erf( sigma      / sqrtDt );
    const Real term8 =  (sigmacb + rcb)     * boost::math::erf((r + sigma) / sqrtDt4);

    return (term1 * (term2 + term3 + term4 + term5) + term6 + term7 + term8) / 6;
}

Real drawR_gbd_3D(const Real sigma, const Real t, const Real D, const Real rnd) noexcept
{
    constexpr Real abs_tol = 1e-18;
    constexpr Real rel_tol = 1e-12;
    const Real target = Igbd_3d(sigma, t, D) * rnd; // rnd is in [0.0, 1.0)

    const Real low  = sigma;
    const Real high = sigma + 10 * std::sqrt(6 * D * t);

    const auto tol = [=](const Real a, const Real b) noexcept {
        return std::abs(a - b) < abs_tol || std::abs(a / b - 1.0) < rel_tol;
    };
    const auto Reqn = [=](const Real x) noexcept {
        return Igbd_r_3d(x, sigma, t, D) - target;
    };

    constexpr std::size_t max_iteration = 100;
    std::size_t iterated = max_iteration;

    const std::pair<Real, Real> t_range =
        boost::math::tools::toms748_solve(Reqn, low, high, tol, iterated);

    if(iterated == max_iteration)
    {
        throw_exception<std::runtime_error>(
            "ngfrd::BDPropagator::bd_math::drawR_gbd_3D: did not converge. "
            "(rnd=", rnd, ", r12=", sigma, ", dt=", t, ", D12=", D, ")");
    }
    return t_range.first;
}
} // bd_math
} // ngfrd
} // ecell4
