#ifndef ECELL4_NGFRD_SINGLE_CIRCULAR_DOMAIN_HPP
#define ECELL4_NGFRD_SINGLE_CIRCULAR_DOMAIN_HPP
#include <ecell4/core/Identifier.hpp>
#include <ecell4/core/Particle.hpp>
#include <ecell4/ngfrd/ShellID.hpp>
#include <ecell4/ngfrd/Shell.hpp>
#include <greens_functions/GreensFunction2DAbsSym.hpp>
#include <cstdint>

namespace ecell4
{
namespace ngfrd
{

// 2D, circular shell
class SingleCircularDomain
{
  public:

    enum class EventKind : std::uint8_t
    {
        Escape,
        Reaction,
        Unknown,
    };

  public:

    SingleCircularDomain()
        : kind_(EventKind::Unknown), dt_(0.0), begin_time_(0.0),
          shell_(Circle(0.0, Real3(0.0, 0.0, 0.0), Real3(0.0, 0.0, 0.0)), FaceID{}),
          gf_(0.0, 0.0)
    {}
    SingleCircularDomain(
            const EventKind kind, const Real dt, const Real begin,
            const ShellID&  shid, const CircularShell& sh,
            const ParticleID& pid, const Real D, const Real effective_radius,
            greens_functions::GreensFunction2DAbsSym gf)
        : kind_(kind), dt_(dt), begin_time_(begin), shell_id_(shid), shell_(sh),
          particle_id_(pid), gf_(std::move(gf))
    {}
    ~SingleCircularDomain() = default;

    ShellID&              shell_id()          noexcept {return shell_id_;}
    ShellID const&        shell_id()    const noexcept {return shell_id_;}
    CircularShell&        shell()             noexcept {return shell_;}
    CircularShell const&  shell()       const noexcept {return shell_;}
    ParticleID&           particle_id()       noexcept {return particle_id_;}
    ParticleID const&     particle_id() const noexcept {return particle_id_;}
    FaceID&               face_id()           noexcept {return shell_.fid();}
    FaceID const&         face_id()     const noexcept {return shell_.fid();}

    Real& dt()       noexcept {return dt_;}
    Real  dt() const noexcept {return dt_;}
    Real& begin_time()       noexcept {return begin_time_;}
    Real  begin_time() const noexcept {return begin_time_;}

    EventKind  eventkind() const {return kind_;}
    EventKind& eventkind()       {return kind_;}

    greens_functions::GreensFunction2DAbsSym const& gf() const noexcept {return gf_;}
    greens_functions::GreensFunction2DAbsSym&       gf()       noexcept {return gf_;}

    constexpr std::size_t num_shells()   const noexcept {return 1;}
    constexpr std::size_t multiplicity() const noexcept {return 1;}

  private:

    EventKind      kind_;
    Real           dt_;
    Real           begin_time_;
    ShellID        shell_id_;
    CircularShell  shell_;
    ParticleID     particle_id_;
    greens_functions::GreensFunction2DAbsSym gf_;
};

} // ngfrd
} // ecell4
#endif /* ECELL4_NGFRD_SINGLE_SPHERICAL_DOMAIN_HPP */
