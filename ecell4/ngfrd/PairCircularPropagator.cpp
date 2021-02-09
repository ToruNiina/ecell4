#include <ecell4/core/exceptions.hpp>
#include <ecell4/ngfrd/BDMath.hpp>
#include <ecell4/ngfrd/PairCircularPropagator.hpp>
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
PairCircularPropagator::com_escape(const PairCircularDomain& dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    constexpr Real pi = boost::math::constants::pi<Real>();

    const auto pid1 = dom.particle1_id();
    const auto pid2 = dom.particle2_id();
    assert(world_.on_which_face(pid1));
    assert(world_.on_which_face(pid2));

    Particle p1 = world_.get_particle(pid1).second;
    Particle p2 = world_.get_particle(pid2).second;

    const Real D1 = p1.D();
    const Real D2 = p2.D();
    const Real D12 = D1 + D2;

    // calculate CoM position

    const Real R_com = dom.com_radius();
    const auto com = ecell4::polygon::travel(world_.polygon(),
            std::make_pair(dom.shell().shape().center(), dom.shell().fid()),
            this->draw_2D_displacement(dom.shell().fid(), R_com));

    // calculate inter-particle vector (p1 -> p2)

    const auto& gf_ipv   = dom.gf_ipv();
    const Real R_ipv     = gf_ipv.drawR(rng_.uniform(0.0, 1.0), dom.dt());
    const Real theta_ipv = gf_ipv.drawTheta(rng_.uniform(0.0, 1.0), R_ipv, dom.dt());
    assert(p1.radius() + p2.radius() <= R_ipv);

    const Real3 ipv = this->generate_ipv(dom.shell().fid(), dom.ipv(), R_ipv, theta_ipv);

    // update particles in the world

    const auto p1_newpos = ecell4::polygon::travel(world_.polygon(), com, ipv * (-D1 / D12));
    const auto p2_newpos = ecell4::polygon::travel(world_.polygon(), com, ipv * ( D2 / D12));

    p1.position() = p1_newpos.first;
    p2.position() = p2_newpos.first;

    world_.update_particle_2D(pid1, p1, p1_newpos.second);
    world_.update_particle_2D(pid2, p2, p2_newpos.second);

    return boost::container::static_vector<ParticleID, 2>{{pid1, pid2}};
}

boost::container::static_vector<ParticleID, 3>
PairCircularPropagator::ipv_escape(const PairCircularDomain& dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    constexpr Real pi = boost::math::constants::pi<Real>();

    const auto pid1 = dom.particle1_id();
    const auto pid2 = dom.particle2_id();
    assert(world_.on_which_face(pid1));
    assert(world_.on_which_face(pid2));

    Particle p1 = world_.get_particle(pid1).second;
    Particle p2 = world_.get_particle(pid2).second;

    const Real D1 = p1.D();
    const Real D2 = p2.D();
    const Real D12 = D1 + D2;

    // calculate CoM position

    const auto& gf_com = dom.gf_com();
    const Real R_com   = gf_com.drawR(rng_.uniform(0.0, 1.0), dom.dt());
    const auto com = ecell4::polygon::travel(world_.polygon(),
            std::make_pair(dom.shell().shape().center(), dom.shell().fid()),
            this->draw_2D_displacement(dom.shell().fid(), R_com));

    // calculate inter-particle vector (p1 -> p2)

    const auto& gf_ipv   = dom.gf_ipv();
    const Real R_ipv     = dom.ipv_radius();
    const Real theta_ipv = gf_ipv.drawTheta(rng_.uniform(0.0, 1.0), R_ipv, dom.dt());
    assert(p1.radius() + p2.radius() <= R_ipv);

    const Real3 ipv = this->generate_ipv(dom.shell().fid(), dom.ipv(), R_ipv, theta_ipv);

    // update particles in the world

    const auto p1_newpos = ecell4::polygon::travel(world_.polygon(), com, ipv * (-D1 / D12));
    const auto p2_newpos = ecell4::polygon::travel(world_.polygon(), com, ipv * ( D2 / D12));

    p1.position() = p1_newpos.first;
    p2.position() = p2_newpos.first;

    world_.update_particle_2D(pid1, p1, p1_newpos.second);
    world_.update_particle_2D(pid2, p2, p2_newpos.second);

    return boost::container::static_vector<ParticleID, 2>{{pid1, pid2}};
}

boost::container::static_vector<ParticleID, 3>
PairCircularPropagator::single_reaction(const PairCircularDomain& dom,
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
PairCircularPropagator::attempt_1to1_reaction(const PairCircularDomain& dom,
        const ParticleID& pid, const Particle& p, const ReactionRule& rule,
        const ParticleID& other)
{
    assert(rule.products().size() == 1);

    const auto species_new = rule.products().front();
    const auto molinfo     = world_.get_molecule_info(species_new);
    const Real radius_new  = molinfo.radius;
    const Real D_new       = molinfo.D;

    FaceID   newfid = this->world_.on_which_face(pid).value();
    Particle newp(species_new, p.position(), radius_new, D_new);

    // first, check if the new particle collides with the paired particle.
    const auto partner = this->world_.get_particle(other).second;
    const auto partner_fid = this->world_.on_which_face(other).value();

    // if it collides, the reaction should be rejected.
    if(ecell4::polygon::distance(this->world_.polygon(),
            std::make_pair(newp.position(), newfid),
            std::make_pair(partner.position(), partner_fid))
            <= radius_new + partner.radius())
    {
        this->rejected_move_count_ += 1;
        return boost::container::static_vector<ParticleID, 3>{pid, other};
    }

    // next, check if the new (different-sized) particle sticks out of the shell.
    if( ! this->is_inside_of_shell(dom, std::make_pair(newp.position(), newfid), newp.radius()))
    {
        // determine positions of particles in overlapping domains.
        // check 3D domain as a spherical particle, and 2D domain as a 2D patch
        sim_.determine_positions_3D(newp.position(), radius_new, self_id_);
        sim_.determine_positions_2D(std::make_pair(newp.position(), newfid),
                                    radius_new, self_id_);

        // then check if there is any particle that overlaps.
        if(world_.has_overlapping_particles_3D(newp.position(), radius_new, pid) ||
           world_.has_overlapping_particles_2D(std::make_pair(newp.position(), newfid), radius_new, pid))
        {
            this->rejected_move_count_ += 1;
            return boost::container::static_vector<ParticleID, 3>{pid, other};
        }
    }

    // it passes all the checks. update the particle species.
    // (the `other` particle is already updated via `propagate` function.)
    world_.update_particle_2D(pid, newp, newfid);

    last_reactions_.emplace_back(rule, make_unimolecular_reaction_info(
                world_.t(), pid, p, pid, newp));
    return boost::container::static_vector<ParticleID, 3>{pid, other};
}

boost::container::static_vector<ParticleID, 3>
PairCircularPropagator::attempt_1to2_reaction(const PairCircularDomain& dom,
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
        throw_exception<NotSupported>("ngfrd::PairCircularPropagator attempts"
            " unbinding reaction but both particles are immovable => ",
            rule.as_string());
    }

    const auto fid         = this->world_.on_which_face(pid).value();
    const auto partner     = this->world_.get_particle(other).second;
    const auto partner_fid = this->world_.on_which_face(other).value();

    const Real partner_distance_threshold_1 = partner.radius() + r1;
    const Real partner_distance_threshold_2 = partner.radius() + r2;

    auto pos1_new = std::make_pair(p.position(), fid);
    auto pos2_new = std::make_pair(p.position(), fid);

    const Real separation_length = r12 * NGFRDSimulator::SAFETY_EXPAND;
    std::size_t separation_count = 1 + max_retry_count_;
    while(separation_count != 0)
    {
        --separation_count;

        const Real3 ipv = this->draw_2D_displacement(fid, separation_length);

        Real3 disp1 = ipv * (D1 / D12);
        Real3 disp2 = disp1 - ipv; // disp1 + (-disp2) = ipv

        pos1_new = ecell4::polygon::travel(world_.polygon(), pos1_new, disp1);
        pos2_new = ecell4::polygon::travel(world_.polygon(), pos2_new, disp2);

        // numerical-error-phobia
        if(ecell4::polygon::distance(world_.polygon(), pos1_new, pos2_new) < r12)
        {
            continue;
        }
        if(ecell4::polygon::distance(world_.polygon(),
                pos1_new, std::make_pair(partner.position(), partner_fid)) <
                partner_distance_threshold_1)
        {
            continue;
        }
        if(ecell4::polygon::distance(world_.polygon(),
                pos1_new, std::make_pair(partner.position(), partner_fid)) <
                partner_distance_threshold_2)
        {
            continue;
        }

        if(this->is_inside_of_shell(dom, pos1_new, r1))
        {
            sim_.determine_positions_2D(pos1_new,       r1, self_id_);
            sim_.determine_positions_3D(pos1_new.first, r1, self_id_);
        }
        if(this->is_inside_of_shell(dom, pos2_new, r2))
        {
            sim_.determine_positions_2D(pos2_new,       r2, self_id_);
            sim_.determine_positions_3D(pos2_new.first, r2, self_id_);
        }

        if(world_.has_overlapping_particles_2D(pos1_new, r1, pid) ||
           world_.has_overlapping_particles_2D(pos2_new, r2, pid))
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

    // burst domains around the reactants.
    if(this->is_inside_of_shell(dom, pos1_new, r1))
    {
        sim_.determine_positions_2D(pos1_new,       r1, self_id_);
        sim_.determine_positions_3D(pos1_new.first, r1, self_id_);
    }
    if(this->is_inside_of_shell(dom, pos2_new, r2))
    {
        sim_.determine_positions_2D(pos2_new,       r2, self_id_);
        sim_.determine_positions_3D(pos2_new.first, r2, self_id_);
    }
    if(world_.has_overlapping_particles_2D(pos1_new, r1, pid) ||
       world_.has_overlapping_particles_2D(pos2_new, r2, pid) ||
       world_.has_overlapping_particles_3D(pos1_new.first, r1, pid) ||
       world_.has_overlapping_particles_3D(pos2_new.first, r2, pid))
    {
        this->rejected_move_count_ += 1;
        return boost::container::static_vector<ParticleID, 2>{pid};
    }

    Particle p1_new(sp1, pos1_new.first, r1, D1);
    Particle p2_new(sp2, pos2_new.first, r2, D2);

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
PairCircularPropagator::pair_reaction(const PairCircularDomain& dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    constexpr Real pi = boost::math::constants::pi<Real>();

    // The logic is different from this->propagate() because the ipv length is
    // at sigma. Since the value of R is at the boundary, the probability
    // distribution of ipv theta differs from normal cases.

    const auto pid1 = dom.particle1_id();
    const auto pid2 = dom.particle2_id();
    assert(world_.on_which_face(pid1).has_value());
    assert(world_.on_which_face(pid2).has_value());

    auto p1 = world_.get_particle(pid1).second;
    auto p2 = world_.get_particle(pid2).second;

    const auto D1  = p1.D();
    const auto D2  = p2.D();
    const auto D12 = D1 + D2;

    // determine com and ipv
    const auto& gf_com = dom.gf_com();
    const Real  R_com  = gf_com.drawR(rng_.uniform(0.0, 1.0), dom.dt());

    const auto com = ecell4::polygon::travel(world_.polygon(),
            std::make_pair(dom.shell().shape().center(), dom.shell().fid()),
            this->draw_2D_displacement(dom.shell().fid(), R_com));

    const auto& gf_ipv    = dom.gf_ipv();
    const Real  R_ipv     = (p1.radius() + p2.radius()) * NGFRDSimulator::SAFETY_EXPAND;
    const Real  theta_ipv = gf_ipv.drawTheta(rng_.uniform(0.0, 1.0), p1.radius() + p2.radius(), dom.dt());
    const Real3 ipv       = this->generate_ipv(dom.shell().fid(), dom.ipv(), R_ipv, theta_ipv);

    // determine the positions of particles immediately before the reaction

    const auto p1_newpos = ecell4::polygon::travel(world_.polygon(), com, ipv * (-D1 / D12));
    const auto p2_newpos = ecell4::polygon::travel(world_.polygon(), com, ipv * ( D2 / D12));

    p1.position() = p1_newpos.first;
    p2.position() = p2_newpos.first;

    // XXX Since rejection of pair reaction is unlikely (the product particle
    // should be implausibly large), we here do not update the positions of those
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
            return attempt_2to1_reaction(dom, pid1, p1, p1_newpos.second,
                                              pid2, p2, p2_newpos.second, rule, com);
        }
        default:
        {
            throw_exception<NotSupported>("ngfrd::PairCircularPropagator: "
                    "invalid number of products, 2 -> 2 reaction");
        }
    }
}

boost::container::static_vector<ParticleID, 3>
PairCircularPropagator::attempt_2to1_reaction(const PairCircularDomain& dom,
        const ParticleID& pid1, const Particle& p1, const FaceID& fid1,
        const ParticleID& pid2, const Particle& p2, const FaceID& fid2,
        const ReactionRule& rule, const std::pair<Real3, FaceID>& com)
{
    ECELL4_NGFRD_LOG_FUNCTION();

    assert(rule.products().size() == 1);

    const auto species_new = rule.products().front();
    const auto molinfo     = world_.get_molecule_info(species_new);
    const Real radius_new  = molinfo.radius;
    const Real D_new       = molinfo.D;

    Particle newp(species_new, com.first, radius_new, D_new);

    // Although it's very unlikely, since a product can be larger than the reactant,
    // there is a small possibility that the product sticks out of the shell.
    if( ! this->is_inside_of_shell(dom, com, newp.radius()))
    {
        // determine positions of particles in overlapping domains
        sim_.determine_positions_3D(newp.position(), radius_new, self_id_);
        sim_.determine_positions_2D(com, radius_new, self_id_);

        // then check if there is any particle that overlaps
        if(world_.has_overlapping_particles_3D(newp.position(), radius_new, pid1, pid2) ||
           world_.has_overlapping_particles_2D(com, radius_new, pid1, pid2))
        {
            // since it does not update the positions of particles, we need to
            // move particles to the state immediately before reaction.
            world_.update_particle_2D(pid1, p1, fid1);
            world_.update_particle_2D(pid2, p2, fid2);

            this->rejected_move_count_ += 1;
            return boost::container::static_vector<ParticleID, 3>{pid1, pid2};
        }
    }

    // looks okay.
    world_.update_particle_2D(pid1, newp, com.second);

    last_reactions_.emplace_back(rule, make_binding_reaction_info(
                world_.t(), pid1, p1, pid2, p2, pid1, newp));
    return boost::container::static_vector<ParticleID, 3>{pid1};
}

std::array<ParticleID, 2>
PairCircularPropagator::propagate(const PairCircularDomain& dom, const Real dt)
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

    // calculate CoM position

    const auto& gf_com = dom.gf_com();
    const Real R_com   = gf_com.drawR(rng_.uniform(0.0, 1.0), dt);

    const auto com = ecell4::polygon::travel(world_.polygon(),
            std::make_pair(dom.shell().shape().center(), dom.shell().fid()),
            this->draw_2D_displacement(dom.shell().fid(), R_com));

    // calculate inter-particle vector (p1 -> p2)

    const auto& gf_ipv   = dom.gf_ipv();
    const Real R_ipv     = gf_ipv.drawR(rng_.uniform(0.0, 1.0), dt);
    const Real theta_ipv = gf_ipv.drawTheta(rng_.uniform(0.0, 1.0), R_ipv, dt);

    assert(p1.radius() + p2.radius() <= R_ipv);

    const Real3 ipv = this->generate_ipv(dom.shell().fid(), dom.ipv(), R_ipv, theta_ipv);

    // update particles in world

    const auto p1_newpos = ecell4::polygon::travel(world_.polygon(), com, ipv * (-D1 / D12));
    const auto p2_newpos = ecell4::polygon::travel(world_.polygon(), com, ipv * ( D2 / D12));

    p1.position() = p1_newpos.first;
    p2.position() = p2_newpos.first;

    world_.update_particle_2D(pid1, p1, p1_newpos.second);
    world_.update_particle_2D(pid2, p2, p2_newpos.second);

    return std::array<ParticleID, 2>{{pid1, pid2}};
}

ReactionRule const& PairCircularPropagator::determine_single_reaction_rule(
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

ReactionRule const& PairCircularPropagator::determine_pair_reaction_rule(
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

bool PairCircularPropagator::is_inside_of_shell(const PairCircularDomain& dom,
        const std::pair<Real3, FaceID>& pos, const Real radius) const
{
    if(dom.shell().thickness() < radius)
    {
        return false;
    }
    const Real dist = ecell4::polygon::distance(world_.polygon(), pos,
            std::make_pair(dom.shell().position(), dom.shell().fid()));

    return dist + radius < dom.shell().shape().radius();
}

Real3 PairCircularPropagator::generate_ipv(const FaceID& fid,
        const Real3& ipv, const Real r, const Real theta) const
{
    const auto& tri  = world_.polygon().triangle_at(fid);
    return rotate(theta, tri.normal(), ipv * (r / length(ipv)));
}

Real3 PairCircularPropagator::draw_2D_displacement(const FaceID& fid, const Real len)
{
    constexpr Real pi = boost::math::constants::pi<Real>();

    const auto& tri  = world_.polygon().triangle_at(fid);
    const auto& disp = (tri.represent() / length(tri.represent())) * len;

    return rotate(rng_.uniform(0.0, 2 * pi), tri.normal(), disp);
}

} // ngfrd
} // ecell4
