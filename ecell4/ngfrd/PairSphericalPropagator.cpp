#include <ecell4/core/exceptions.hpp>
#include <ecell4/ngfrd/BDMath.hpp>
#include <ecell4/ngfrd/PairSphericalPropagator.hpp>
#include <ecell4/ngfrd/NGFRDSimulator.hpp>
#include <ecell4/ngfrd/Logger.hpp>

#include <greens_functions/freeFunctions.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/tools/roots.hpp>
#include <numeric>

namespace ecell4
{
namespace ngfrd
{

boost::container::static_vector<ParticleID, 3>
PairSphericalPropagator::com_escape(const PairSphericalDomain& dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    constexpr Real pi = boost::math::constants::pi<Real>();

    const auto pid1 = dom.particle1_id();
    const auto pid2 = dom.particle2_id();
    assert( ! world_.on_which_face(pid1));
    assert( ! world_.on_which_face(pid2));

    Particle p1 = world_.get_particle(pid1).second;
    Particle p2 = world_.get_particle(pid2).second;

    const Real D1 = p1.D();
    const Real D2 = p2.D();
    const Real D12 = D1 + D2;

    const Real R_com = dom.com_radius();

    // after adding ipv, we will apply boundary. (this optimization is correct
    // as long as the domain is smaller than the whole world, and this
    // assumption is always correct.)
    const Real3 com = dom.shell().shape().center() + rng_.direction3d(R_com);

    const auto& gf_ipv   = dom.gf_ipv();
    const Real R_ipv     = gf_ipv.drawR(rng_.uniform(0.0, 1.0), dom.dt());
    const Real theta_ipv = gf_ipv.drawTheta(rng_.uniform(0.0, 1.0), R_ipv, dom.dt());
    const Real phi_ipv   = rng_.uniform(0.0, 2 * pi);

    assert(p1.radius() + p2.radius() <= R_ipv);

    // from p1 -> p2
    const Real3 ipv = this->generate_ipv(dom.ipv(), R_ipv, theta_ipv, phi_ipv);

    p1.position() = world_.boundary().apply_boundary(com - ipv * (D1 / D12));
    p2.position() = world_.boundary().apply_boundary(com + ipv * (D2 / D12));

    world_.update_particle_3D(pid1, p1);
    world_.update_particle_3D(pid2, p2);

    return boost::container::static_vector<ParticleID, 2>{{pid1, pid2}};
}

boost::container::static_vector<ParticleID, 3>
PairSphericalPropagator::ipv_escape(const PairSphericalDomain& dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    constexpr Real pi = boost::math::constants::pi<Real>();

    const auto pid1 = dom.particle1_id();
    const auto pid2 = dom.particle2_id();
    assert( ! world_.on_which_face(pid1));
    assert( ! world_.on_which_face(pid2));

    Particle p1 = world_.get_particle(pid1).second;
    Particle p2 = world_.get_particle(pid2).second;

    const Real D1 = p1.D();
    const Real D2 = p2.D();
    const Real D12 = D1 + D2;

    const auto& gf_com = dom.gf_com();
    const Real R_com   = gf_com.drawR(rng_.uniform(0.0, 1.0), dom.dt());

    // later we will do apply_boundary. (this optimization is correct as long as
    // the domain is smaller than the whole world, and this assumption is always
    // correct.)
    const Real3 com = dom.shell().shape().center() + rng_.direction3d(R_com);

    const auto& gf_ipv   = dom.gf_ipv();
    const Real R_ipv     = dom.ipv_radius();
    const Real theta_ipv = gf_ipv.drawTheta(rng_.uniform(0.0, 1.0), R_ipv, dom.dt());
    const Real phi_ipv   = rng_.uniform(0.0, 2 * pi);

    assert(p1.radius() + p2.radius() <= R_ipv);

    // from p1 -> p2
    const Real3 ipv = this->generate_ipv(dom.ipv(), R_ipv, theta_ipv, phi_ipv);

    p1.position() = world_.boundary().apply_boundary(com - ipv * (D1 / D12));
    p2.position() = world_.boundary().apply_boundary(com + ipv * (D2 / D12));

    world_.update_particle_3D(pid1, p1);
    world_.update_particle_3D(pid2, p2);

    return boost::container::static_vector<ParticleID, 2>{{pid1, pid2}};
}

boost::container::static_vector<ParticleID, 3>
PairSphericalPropagator::single_reaction(const PairSphericalDomain& dom,
        const ParticleID& pid, const ParticleID& other)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    // move particle and update world
    this->propagate(dom, dom.dt());

    // fetch propagated particle positions
    const auto p     = world_.get_particle(pid).second;
    const auto& rule = this->determine_single_reaction_rule(
            p.species(), rng_.uniform(0.0, 1.0));

    switch(rule.products().size())
    {
        case 0:
        {
            world_.remove_particle(pid);
            last_reactions_.emplace_back(rule,
                    make_degradation_reaction_info(world_.t(), pid, p));
            return boost::container::static_vector<ParticleID, 3>{other};
        }
        case 1:
        {
            return this->attempt_1to1_reaction(dom, pid, p, rule, other);
        }
        case 2:
        {
            return this->attempt_1to2_reaction(dom, pid, p, rule, other);
        }
        default:
        {
            throw_exception<NotSupported>("invalid number of products");
        }
    }
}

boost::container::static_vector<ParticleID, 3>
PairSphericalPropagator::attempt_1to1_reaction(const PairSphericalDomain& dom,
        const ParticleID& pid, const Particle& p, const ReactionRule& rule,
        const ParticleID& other)
{
    assert(rule.products().size() == 1);

    const auto species_new = rule.products().front();
    const auto molinfo     = world_.get_molecule_info(species_new);
    const Real radius_new  = molinfo.radius;
    const Real D_new       = molinfo.D;

    Particle newp(species_new, p.position(), radius_new, D_new);

    // first, check if the new particle collides with the paired particle.
    const auto b       = this->world_.boundary();
    const auto partner = this->world_.get_particle(other).second;

    // if it collides, the reaction should be rejected.
    if(length_sq(newp.position() - b.periodic_transpose(partner.position(), newp.position()))
        <= (radius_new + partner.radius()) * (radius_new + partner.radius()))
    {
        this->rejected_move_count_ += 1;
        return boost::container::static_vector<ParticleID, 3>{pid, other};
    }

    // next, check if the new (different-sized) particle sticks out of the shell.
    if( ! this->is_inside_of_shell(dom, newp.position(), newp.radius()))
    {
        if(world_.has_overlapping_faces(newp.position(), radius_new))
        {
            this->rejected_move_count_ += 1;
            return boost::container::static_vector<ParticleID, 3>{pid, other};
        }

        // determine positions of particles in overlapping domains
        sim_.determine_positions_XD(newp.position(), radius_new, self_id_);

        // then check if there is any particle that overlaps
        if(world_.has_overlapping_particles_XD(
                    newp.position(), radius_new, /*ignore = */ pid))
        {
            this->rejected_move_count_ += 1;
            return boost::container::static_vector<ParticleID, 3>{pid, other};
        }
    }

    // it passes all the checks. update the particle species.
    // (the `other` particle is already updated via `propagate` function.)
    world_.update_particle_3D(pid, newp);

    last_reactions_.emplace_back(rule, make_unimolecular_reaction_info(
                world_.t(), pid, p, pid, newp));
    return boost::container::static_vector<ParticleID, 3>{pid, other};
}

boost::container::static_vector<ParticleID, 3>
PairSphericalPropagator::attempt_1to2_reaction(const PairSphericalDomain& dom,
        const ParticleID& pid, const Particle& p, const ReactionRule& rule,
        const ParticleID& other)
{
    assert(rule.products().size() == 2);

    const auto sp1      = rule.products().at(0);
    const auto sp2      = rule.products().at(1);
    const auto molinfo1 = world_.get_molecule_info(sp1);
    const auto molinfo2 = world_.get_molecule_info(sp2);

    const Real D1  = molinfo1.D;
    const Real D2  = molinfo2.D;
    const Real r1  = molinfo1.radius;
    const Real r2  = molinfo2.radius;
    const Real D12 = D1 + D2;
    const Real r12 = r1 + r2;

    if(D1 == 0 && D2 == 0)
    {
        throw_exception<NotSupported>("ngfrd::PairSphericalPropagator attempts"
            " unbinding reaction but both particles are immovable => ",
            rule.as_string());
    }

    const auto boundary = world_.boundary();
    const auto partner = this->world_.get_particle(other).second;

    const Real partner_distance_threshold_1_sq = (partner.radius() + r1) * (partner.radius() + r1);
    const Real partner_distance_threshold_2_sq = (partner.radius() + r2) * (partner.radius() + r2);

    Real3 pos1_new(p.position());
    Real3 pos2_new(p.position());

    const Real separation_length = r12 * NGFRDSimulator::SAFETY_EXPAND;
    std::size_t separation_count = 1 + max_retry_count_;
    while(separation_count != 0)
    {
        --separation_count;

        // Here, the time when the reaciton happens is precisely determined.
        // So, unlike BD, we don't need to consider the inter-particle distance.
        // We can place the two particles in contact with each other.

        const Real  R   = separation_length;
        const Real3 ipv = rng_.direction3d(R);

        Real3 disp1 = ipv * (D1 / D12);
        Real3 disp2 = disp1 - ipv; // disp1 + (-disp2) = ipv

        pos1_new = boundary.apply_boundary(p.position() + disp1);
        pos2_new = boundary.apply_boundary(p.position() + disp2);

        // numerical-error-phobia
        if(length_sq(boundary.periodic_transpose(pos1_new, pos2_new) - pos2_new)
                < r12 * r12)
        {
            continue;
        }

        // check if any of the products collide with the paired (non-reactant) particle.
        if(length_sq(boundary.periodic_transpose(pos1_new, partner.position()) - partner.position())
                < partner_distance_threshold_1_sq ||
           length_sq(boundary.periodic_transpose(pos2_new, partner.position()) - partner.position())
                < partner_distance_threshold_2_sq)
        {
            continue;
        }

        // burst domains around the reactants, if needed.
        if( ! this->is_inside_of_shell(dom, pos1_new, r1))
        {
            sim_.determine_positions_XD(pos1_new, r1, self_id_);
        }
        if( ! this->is_inside_of_shell(dom, pos2_new, r2))
        {
            sim_.determine_positions_XD(pos2_new, r2, self_id_);
        }

        if(world_.has_overlapping_faces(pos1_new, r1) ||
           world_.has_overlapping_faces(pos2_new, r2))
        {
            continue;
        }
        // since we already check the paired (non-reactant) particle, we can skip it.
        if(world_.has_overlapping_particles_XD(pos1_new, r1, pid, other) ||
           world_.has_overlapping_particles_XD(pos2_new, r2, pid, other))
        {
            continue;
        }
    }
    if(separation_count == 0)
    {
        // could not find an appropreate configuration in max_retry_count_.
        this->rejected_move_count_ += 1;
        return boost::container::static_vector<ParticleID, 3>{pid, other};
    }

    Particle p1_new(sp1, pos1_new, r1, D1);
    Particle p2_new(sp2, pos2_new, r2, D2);

    const auto result1 = world_.update_particle_3D(pid, p1_new);
    const auto result2 = world_.new_particle_3D   (     p2_new);

    assert(not result1); // should be already exist (it returns true if it's new)
    assert(result2.second); // should succeed

    const auto pid2 = result2.first.first;

    last_reactions_.emplace_back(rule, make_unbinding_reaction_info(world_.t(),
                pid, p, pid, p1_new, pid2, p2_new));

    return boost::container::static_vector<ParticleID, 3>{pid, pid2, other};
}

boost::container::static_vector<ParticleID, 3>
PairSphericalPropagator::pair_reaction(const PairSphericalDomain& dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    constexpr Real pi = boost::math::constants::pi<Real>();

    const auto pid1 = dom.particle1_id();
    const auto pid2 = dom.particle2_id();
    assert( ! world_.on_which_face(pid1));
    assert( ! world_.on_which_face(pid2));

    auto p1 = world_.get_particle(pid1).second;
    auto p2 = world_.get_particle(pid2).second;

    const auto D1  = p1.D();
    const auto D2  = p2.D();
    const auto D12 = D1 + D2;

    // determine com and ipv
    const auto& gf_com = dom.gf_com();
    const Real  R_com  = gf_com.drawR(rng_.uniform(0.0, 1.0), dom.dt());
    const Real3 com    = dom.shell().shape().center() + rng_.direction3d(R_com);

    const auto& gf_ipv    = dom.gf_ipv();
    const Real  R_ipv     = (p1.radius() + p2.radius()) * NGFRDSimulator::SAFETY_EXPAND;
    const Real  theta_ipv = gf_ipv.drawTheta(rng_.uniform(0.0, 1.0), p1.radius() + p2.radius(), dom.dt());
    const Real  phi_ipv   = rng_.uniform(0.0, 2 * pi);
    const Real3 ipv       = this->generate_ipv(dom.ipv(), R_ipv, theta_ipv, phi_ipv);

    // determine the positions of particles immediately before the reaction

    const auto b = world_.boundary();
    p1.position() = b.apply_boundary(com - ipv * (D1 / D12));
    p2.position() = b.apply_boundary(com + ipv * (D2 / D12));

    // XXX since rejection of pair reaction is unlikely (the product particle
    // should be implausibly large, we here do not update the positions of those
    // particles in the `world_`. If rejection happens, then we need to update.

    // determine reaction
    const auto& rule = this->determine_pair_reaction_rule(
            p1.species(), p2.species(), rng_.uniform(0.0, 1.0));

    switch(rule.products().size())
    {
        case 0:
        {
            world_.remove_particle(pid1);
            world_.remove_particle(pid2);
            last_reactions_.emplace_back(rule, make_degradation_reaction_info(
                        world_.t(), pid1, p1, pid2, p2));
            return boost::container::static_vector<ParticleID, 3>{/*empty*/};
        }
        case 1:
        {
            return attempt_2to1_reaction(dom, pid1, p1, pid2, p2, rule, com);
        }
        default:
        {
            throw_exception<NotSupported>("ngfrd::PairSphericalPropagator: "
                    "invalid number of products, 2 -> 2 reaction");
        }
    }
}

boost::container::static_vector<ParticleID, 3>
PairSphericalPropagator::attempt_2to1_reaction(const PairSphericalDomain& dom,
        const ParticleID& pid1, const Particle& p1,
        const ParticleID& pid2, const Particle& p2,
        const ReactionRule& rule, const Real3& com)
{
    ECELL4_NGFRD_LOG_FUNCTION();

    assert(rule.products().size() == 1);

    const auto species_new = rule.products().front();
    const auto molinfo     = world_.get_molecule_info(species_new);
    const Real radius_new  = molinfo.radius;
    const Real D_new       = molinfo.D;

    Particle newp(species_new, com, radius_new, D_new);

    // Although it's very unlikely, since a product can be larger than the reactant,
    // there is a small possibility that the product sticks out of the shell
    if( ! this->is_inside_of_shell(dom, newp.position(), newp.radius()))
    {
        if(world_.has_overlapping_faces(newp.position(), radius_new))
        {
            // since it does not update the positions of particles, we need to
            // move particles to the state immediately before reaction.
            world_.update_particle_3D(pid1, p1);
            world_.update_particle_3D(pid2, p2);

            this->rejected_move_count_ += 1;
            return boost::container::static_vector<ParticleID, 3>{pid1, pid2};
        }

        // determine positions of particles in overlapping domains
        sim_.determine_positions_XD(newp.position(), radius_new, self_id_);

        // then check if there is any particle that overlaps
        if(world_.has_overlapping_particles_XD(
                    newp.position(), radius_new, /*ignore = */ pid1, pid2))
        {
            // since it does not update the positions of particles, we need to
            // move particles to the state immediately before reaction.
            world_.update_particle_3D(pid1, p1);
            world_.update_particle_3D(pid2, p2);

            this->rejected_move_count_ += 1;
            return boost::container::static_vector<ParticleID, 3>{pid1, pid2};
        }
    }

    // looks okay.
    world_.update_particle_3D(pid1, newp);

    last_reactions_.emplace_back(rule, make_binding_reaction_info(
                world_.t(), pid1, p1, pid2, p2, pid1, newp));
    return boost::container::static_vector<ParticleID, 3>{pid1};
}

std::array<ParticleID, 2>
PairSphericalPropagator::propagate(const PairSphericalDomain& dom, const Real dt)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    constexpr Real pi = boost::math::constants::pi<Real>();

    const auto pid1 = dom.particle1_id();
    const auto pid2 = dom.particle2_id();
    assert( ! world_.on_which_face(pid1));
    assert( ! world_.on_which_face(pid2));

    Particle p1 = world_.get_particle(pid1).second;
    Particle p2 = world_.get_particle(pid2).second;

    const Real D1 = p1.D();
    const Real D2 = p2.D();
    const Real D12 = D1 + D2;

    const auto& gf_com = dom.gf_com();
    const Real R_com   = gf_com.drawR(rng_.uniform(0.0, 1.0), dt);

    const Real3 com = world_.boundary().apply_boundary(
        dom.shell().shape().center() + rng_.direction3d(R_com));

    const auto& gf_ipv   = dom.gf_ipv();
    const Real R_ipv     = gf_ipv.drawR(rng_.uniform(0.0, 1.0), dt);
    const Real theta_ipv = gf_ipv.drawTheta(rng_.uniform(0.0, 1.0), R_ipv, dt);
    const Real phi_ipv   = rng_.uniform(0.0, 2 * pi);

    assert(p1.radius() + p2.radius() <= R_ipv);

    // from p1 -> p2
    const Real3 ipv = this->generate_ipv(dom.ipv(), R_ipv, theta_ipv, phi_ipv);

    p1.position() = world_.boundary().apply_boundary(com - ipv * (D1 / D12));
    p2.position() = world_.boundary().apply_boundary(com + ipv * (D2 / D12));

    world_.update_particle_3D(pid1, p1);
    world_.update_particle_3D(pid2, p2);

    return std::array<ParticleID, 2>{{pid1, pid2}};
}

ReactionRule const& PairSphericalPropagator::determine_single_reaction_rule(
        const Species& sp, const Real rnd)
{
    const auto rules = this->model_.query_reaction_rules(sp);
    assert(!rules.empty());

    if(rules.size() == 1)
    {
        return rules.front();
    }

    const Real k_tot = std::accumulate(rules.begin(), rules.end(), Real(0.0),
            [](const Real k, const ReactionRule& rule) -> Real {
                return k + rule.k();
            });
    const Real threshold = k_tot * rnd;

    Real k_cumm = 0.0;
    for(const auto& rule : rules)
    {
        k_cumm += rule.k();
        if(threshold < k_cumm)
        {
            return rule;
        }
    }
    // workaround for numerical error ...
    return rules.back();
}

ReactionRule const& PairSphericalPropagator::determine_pair_reaction_rule(
        const Species& sp1, const Species& sp2, const Real rnd)
{
    const auto rules = this->model_.query_reaction_rules(sp1, sp2);
    assert(!rules.empty());

    if(rules.size() == 1)
    {
        return rules.front();
    }

    const Real k_tot = std::accumulate(rules.begin(), rules.end(), Real(0.0),
            [](const Real k, const ReactionRule& rule) -> Real {
                return k + rule.k();
            });
    const Real threshold = k_tot * rnd;

    Real k_cumm = 0.0;
    for(const auto& rule : rules)
    {
        k_cumm += rule.k();
        if(threshold < k_cumm)
        {
            return rule;
        }
    }
    // workaround for numerical error ...
    return rules.back();
}

bool PairSphericalPropagator::is_inside_of_shell(
        const PairSphericalDomain& dom, const Real3& center, const Real radius)
{
    const auto b  = this->world_.boundary();
    const Real threshold = dom.shell().shape().radius() - radius;
    const auto dr = b.periodic_transpose(dom.shell().shape().center(), center) - center;
    return length_sq(dr) < (threshold * threshold);
}

Real3 PairSphericalPropagator::generate_ipv(
        const Real3& ipv, const Real r, const Real theta, const Real phi)
{
    constexpr Real pi = boost::math::constants::pi<Real>();

    const auto sin_theta = std::sin(theta);
    const auto cos_theta = std::cos(theta);
    const auto sin_phi   = std::sin(phi);
    const auto cos_phi   = std::cos(phi);

    // an angle between `ipv` and z axis
    const auto angle = std::acos(boost::algorithm::clamp(ipv[2] / length(ipv), -1.0, 1.0));
    if(angle == 0)
    {
        return Real3(r * sin_theta * cos_phi,
                     r * sin_theta * sin_phi,
                     r * cos_theta);
    }
    else if(angle == pi)
    {
        return Real3(r * sin_theta * cos_phi,
                     r * sin_theta * sin_phi,
                    -r * cos_theta);
    }
    const Real3 new_ipv(r * sin_theta * cos_phi,
                        r * sin_theta * sin_phi,
                        r * cos_theta);
    const Real3 cross_z(-ipv[1], ipv[0], 0);

    return rotate(angle, cross_z, new_ipv);
}

} // ngfrd
} // ecell4
