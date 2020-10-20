#include <ecell4/core/exceptions.hpp>
#include <ecell4/ngfrd/SingleSphericalPropagator.hpp>
#include <ecell4/ngfrd/NGFRDSimulator.hpp>
#include <ecell4/ngfrd/Logger.hpp>

#include <greens_functions/freeFunctions.hpp>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/tools/roots.hpp>

namespace ecell4
{
namespace ngfrd
{

boost::container::static_vector<ParticleID, 2>
SingleSphericalPropagator::escape(const SingleSphericalDomain& dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();

    const auto pid = dom.particle_id();
    assert( ! world_.on_which_face(pid));
    Particle p = world_.get_particle(pid).second;

    const Real R = shell_.shape().radius() - p.radius();
    p.position() = world_.boundary().apply_boundary(
        sh.shape().center() + rng_.direction3d(R));

    world_.update_particle(pid, p);

    return boost::container::static_vector<ParticleID, 2>{pid};
}

boost::container::static_vector<ParticleID, 2>
SingleSphericalPropagator::reaction(const SingleSphericalDomain& dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    // move particle and update world
    this->propagate(dom, this->domain_.dt());

    // fetch propagated particle positions
    const auto pid   = dom.particle_id();
    const auto p     = world_.get_particle(pid).second;
    const auto& rule = this->determine_reaction_rule(
            p.species(), rng_.uniform(0.0, 1.0));

    switch(rule.products().size())
    {
        case 0:
        {
            world_.remove_particle(pid);
            last_reactions_.emplace_back(rule,
                    make_degradation_reaction_info(world_.t(), pid, p));
            return boost::container::static_vector<ParticleID, 2>{};
        }
        case 1:
        {
            return this->attempt_1to1_reaction(dom, pid, p, rule);
        }
        case 2:
        {
            return this->attempt_1to2_reaction(dom, pid, p, rule);
        }
        default:
        {
            format_throw<NotSupported>("invalid number of products");
        }
    }
}

boost::container::static_vector<ParticleID, 2>
SingleSphericalPropagator::attempt_1to1_reaction(const SingleSphericalDomain& dom,
        const ParticleID& pid, const Particle& p, const ReactionRule& rule)
{
    assert(rule.products().size() == 1);
    const auto species_new = rule.products().front();
    const auto molinfo     = world_.get_molecule_info(species_new);
    const Real radius_new  = molinfo.radius;
    const Real D_new       = molinfo.D;

    Particle newp(species_new, p.position(), radius_new, D_new);

    // since a product could be larger than the reactant,
    // it may sticks out of the shell.
    if( ! this->is_inside_of_shell(dom, newp.position(), newp.radius()))
    {
        if(world_.has_overlapping_faces(newp.position(), radius_new))
        {
            return false;
        }
        // determine positions of particles in overlapping domains
        sim_.determine_positions_3D(newp.position(), radius_new);

        // then check if there is any particle that overlaps
        if(world_.has_overlapping_particles_3D(
                    newp.position(), radius_new, /*ignore = */ pid))
        {
            return false;
        }
    }

    // looks okay.
    world_.update_particle_3D(pid, newp);

    last_reactions_.emplace_back(rule, make_unimolecular_reaction_info(
                world_.t(), pid, p, pid, newp));
    return boost::container::static_vector<ParticleID, 2>{pid};
}

boost::container::static_vector<ParticleID, 2>
SingleSphericalPropagator::attempt_1to2_reaction(const SingleSphericalDomain& dom,
        const ParticleID& pid, const Particle& p, const ReactionRule& rule)
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
        throw_exception<NotSupported>("ngfrd::SingleSphericalPropagator attempts"
            " unbinding reaction but both particles are immovable => ",
            rule.as_string());
    }

    const auto boundary = world_.boundary();

    Real3 pos1_new(p.position());
    Real3 pos2_new(p.position());

    const Real separation_length = r12 * NGFRDSimulator::SAFETY; // = 1+epsilon

    // XXX Since single shell considers only the radius of reactant, not the
    //     products. The radius of products could be larger than that of
    //     reactant, they may stick out of the shell and overlaps with other
    //     particles or faces.
    //         To make checking easier, we first determine particles around the
    //     reacting volume which is a sphere centering at p.position(), and its
    //     radius is separation_length.
    //         It rarely bursts other domains.
    if(this->is_inside_of_shell(dom, p.position(), separation_length))
    {
        // if reacting volume sticks out of the shell, the products may also
        // stick out of the shell.
        sim_.determine_positions_3D(p.position(), separation_length);
    }

    std::size_t separation_count = 1 + max_retry_count_;
    while(separation_count != 0)
    {
        --separation_count;

        const Real3 ipv = draw_ipv_3D(separation_length, dt_, D12);
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

        // we already bursted other domains that might overlaps with them.

        if(world_.has_overlapping_faces(pos1_new, r1) ||
           world_.has_overlapping_faces(pos2_new, r2))
        {
            continue;
        }
        if(world_.has_overlapping_particles_3D(pos1_new, r1, pid) ||
           world_.has_overlapping_particles_3D(pos2_new, r2, pid))
        {
            continue;
        }
    }
    if(separation_count == 0)
    {
        // could not find an appropreate configuration in max_retry_count_.
        this->rejected_move_count_ += 1;
        return false;
    }

    Particle p1_new(sp1, pos1_new, r1, D1);
    Particle p2_new(sp2, pos2_new, r2, D2);

    const auto result1 = world_.update_particle_3D(pid, p1_new);
    const auto result2 = world_.new_particle_3D   (     p2_new);

    assert(not result1); // should be already exist (it returns true if it's new)
    assert(result2.second); // should succeed

    const auto pid2 = result2.first.first;
    this->particles_.push_back(pid2);

    last_reactions_.emplace_back(rule, make_unbinding_reaction_info(world_.t(),
                pid, p, pid, p1_new, pid2, p2_new));

    return boost::container::static_vector<ParticleID, 2>{pid, pid2};
}

//
// Since proapgate is called only when burst is called, no reaction happens.
// Also, SingleSphericalDomain is formed avoiding collision, no collision happens.
//
ParticleID SingleSphericalPropagator::propagate(const SingleSphericalDomain& dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();

    const auto pid = dom.particle_id();
    assert( ! world_.on_which_face(pid));

    Particle p = world_.get_particle(pid).second;

    const auto& gf = dom.gf();
    const Real R = gf.drawR(rng_.uniform(0.0, 1.0));
    p.position() = world_.boundary().apply_boundary(
        sh.shape().center() + rng_.direction3d(R));

    world_.update_particle(pid, p);
    return pid;
}

ReactionRule const& SingleSphericalPropagator::determine_reaction_rule(
        const Species& sp, const Real rnd)
{
    const auto rules = this->model_.query_reaction_rules(p.species());
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

bool SingleSphericalPropagator::is_inside_of_shell(
        const SingleSphericalDomain& dom, const Real3& center, const Real radius)
{
    const auto b  = this->world_.boundary();
    const Real threshold = dom.shape().radius() - radius;
    const auto dr = b.periodic_transpose(dom.shape().center(), center) - center;
    return length_sq(dr) < (threshold * threshold);
}

} // ngfrd
} // ecell4
