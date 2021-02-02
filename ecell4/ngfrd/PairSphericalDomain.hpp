#ifndef ECELL4_NGFRD_PAIR_SPHERICAL_DOMAIN_HPP
#define ECELL4_NGFRD_PAIR_SPHERICAL_DOMAIN_HPP
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/Particle.hpp>
#include <ecell4/ngfrd/ShellID.hpp>
#include <ecell4/ngfrd/Shell.hpp>
#include <greens_functions/GreensFunction3DRadAbs.hpp>
#include <greens_functions/GreensFunction3DAbsSym.hpp>
#include <cstdint>

namespace ecell4
{
namespace ngfrd
{

// 3D, spherical shell
class PairSphericalDomain
{
  public:
    enum class EventKind : std::uint8_t
    {
        ComEscape,
        IpvEscape,
        SingleReaction1,
        SingleReaction2,
        PairReaction,
        PairEvent, // not yet determined
        Unknown,
    };

  public:

    PairSphericalDomain()
        : kind_(EventKind::Unknown), dt_(0.0), begin_time_(0.0),
          com_radius_(0.0), ipv_radius_(0.0),
          shell_id_(ShellID{}), shell_(Sphere(Real3(0.0, 0.0, 0.0), 0.0)),
          particle1_id_(ParticleID{}), particle2_id_(ParticleID{}),
          gf_ipv_(0.0, 0.0, 0.0, 0.0, 0.0), gf_com_(0.0, 0.0)
    {}
    PairSphericalDomain(
            const EventKind kind, const Real dt, const Real begin,
            const ShellID&  shid, const SphericalShell& sh,
            const Real com_radius, const Real ipv_radius,
            const ParticleID& pid1, const ParticleID& pid2, const Real3& ipv,
            greens_functions::GreensFunction3DAbsSym gf_com,
            greens_functions::GreensFunction3DRadAbs gf_ipv)
        : kind_(kind), dt_(dt), begin_time_(begin),
          com_radius_(com_radius), ipv_radius_(ipv_radius),
          shell_id_(shid), shell_(sh),
          particle1_id_(pid1), particle2_id_(pid2), ipv_(ipv),
          gf_com_(std::move(gf_com)), gf_ipv_(std::move(gf_ipv))
    {}
    ~PairSphericalDomain() = default;

    ShellID&              shell_id()          noexcept {return shell_id_;}
    ShellID const&        shell_id()    const noexcept {return shell_id_;}
    SphericalShell&       shell()             noexcept {return shell_;}
    SphericalShell const& shell()       const noexcept {return shell_;}

    ParticleID&           particle1_id()       noexcept {return particle1_id_;}
    ParticleID const&     particle1_id() const noexcept {return particle1_id_;}
    ParticleID&           particle2_id()       noexcept {return particle2_id_;}
    ParticleID const&     particle2_id() const noexcept {return particle2_id_;}

    Real& dt()       noexcept {return dt_;}
    Real  dt() const noexcept {return dt_;}
    Real& begin_time()       noexcept {return begin_time_;}
    Real  begin_time() const noexcept {return begin_time_;}

    EventKind  eventkind() const {return kind_;}
    EventKind& eventkind()       {return kind_;}

    Real  com_radius() const noexcept {return com_radius_;}
    Real& com_radius()       noexcept {return com_radius_;}
    Real  ipv_radius() const noexcept {return ipv_radius_;}
    Real& ipv_radius()       noexcept {return ipv_radius_;}

    Real3 const& ipv() const noexcept {return ipv_;}

    greens_functions::GreensFunction3DRadAbs const& gf_ipv() const noexcept {return gf_ipv_;}
    greens_functions::GreensFunction3DRadAbs&       gf_ipv()       noexcept {return gf_ipv_;}

    greens_functions::GreensFunction3DAbsSym const& gf_com() const noexcept {return gf_com_;}
    greens_functions::GreensFunction3DAbsSym&       gf_com()       noexcept {return gf_com_;}

    constexpr std::size_t num_shells()   const noexcept {return 1;}
    constexpr std::size_t multiplicity() const noexcept {return 2;}

  private:

    EventKind      kind_;
    Real           dt_;
    Real           begin_time_;
    Real           com_radius_; // effective radius of com domain
    Real           ipv_radius_; // effective radius of ipv domain
    ShellID        shell_id_;
    SphericalShell shell_;
    ParticleID     particle1_id_;
    ParticleID     particle2_id_;
    Real3          ipv_;        // initial inter particle vector
    greens_functions::GreensFunction3DAbsSym gf_com_;
    greens_functions::GreensFunction3DRadAbs gf_ipv_;
};

} // ngfrd
} // ecell4
#endif /* ECELL4_NGFRD_PAIR_SPHERICAL_DOMAIN_HPP */
