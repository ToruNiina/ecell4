#ifndef ECELL4_NGFRD_SINGLE_CONICAL_DOMAIN_HPP
#define ECELL4_NGFRD_SINGLE_CONICAL_DOMAIN_HPP
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/Particle.hpp>
#include <ecell4/ngfrd/ShellID.hpp>
#include <ecell4/ngfrd/Shell.hpp>
#include <greens_functions/GreensFunction2DRefWedgeAbs.hpp>
#include <cstdint>

namespace ecell4
{
namespace ngfrd
{

// 2D, conical shell
class SingleConicalDomain
{
  public:

    enum class EventKind : std::uint8_t
    {
        Escape,
        Reaction,
        Unknown,
    };

  public:

    SingleConicalDomain()
        : kind_(EventKind::Unknown), dt_(0.0), begin_time_(0.0), effective_radius_(0.0),
          shell_(0.0, Circle(0.0, Real3(0.0, 0.0, 0.0), Real3(0.0, 0.0, 0.0)), FaceID{}),
          gf_(0.0, 0.0)
    {}
    SingleConicalDomain(
            const EventKind kind, const Real dt, const Real begin,
            const ShellID&  shid, const ConicalShell& sh,
            const ParticleID& pid, const Real D, const Real effective_radius,
            greens_functions::GreensFunction2DAbsSym gf)
        : kind_(kind), dt_(dt), begin_time_(begin), effective_radius_(effective_radius),
          shell_id_(shid), shell_(sh), particle_id_(pid), gf_(std::move(gf))
    {}
    ~SingleConicalDomain() = default;

    ShellID&              shell_id()          noexcept {return shell_id_;}
    ShellID const&        shell_id()    const noexcept {return shell_id_;}
    ConicalShell&         shell()             noexcept {return shell_;}
    ConicalShell const&   shell()       const noexcept {return shell_;}
    ParticleID&           particle_id()       noexcept {return particle_id_;}
    ParticleID const&     particle_id() const noexcept {return particle_id_;}
    VertexID&             vertex_id()         noexcept {return shell_.vid();}
    VertexID const&       vertex_id()   const noexcept {return shell_.vid();}

    Real& dt()       noexcept {return dt_;}
    Real  dt() const noexcept {return dt_;}
    Real& begin_time()       noexcept {return begin_time_;}
    Real  begin_time() const noexcept {return begin_time_;}

    EventKind  eventkind() const {return kind_;}
    EventKind& eventkind()       {return kind_;}

    Real effective_radius() const noexcept {return effective_radius_;}

    greens_functions::GreensFunction2DRefWedgeAbs const& gf() const noexcept {return gf_;}
    greens_functions::GreensFunction2DRefWedgeAbs&       gf()       noexcept {return gf_;}

    constexpr std::size_t num_shells()   const noexcept {return 1;}
    constexpr std::size_t multiplicity() const noexcept {return 1;}

  private:

    EventKind      kind_;
    Real           dt_;
    Real           begin_time_;
    Real           effective_radius_;
    ShellID        shell_id_;
    ConicalShell   shell_;
    ParticleID     particle_id_;
    greens_functions::GreensFunction2DRefWedgeAbs gf_;
};

} // ngfrd
} // ecell4
#endif /* ECELL4_NGFRD_SINGLE_CONICAL_DOMAIN_HPP */
