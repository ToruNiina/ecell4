#include <ecell4/core/exceptions.hpp>
#include <ecell4/ngfrd/BDMath.hpp>
#include <ecell4/ngfrd/SingleCircularPropagator.hpp>
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

boost::container::static_vector<ParticleID, 2>
SingleCircularPropagator::escape(const SingleCircularDomain& dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();

    const auto pid = dom.particle_id();

    FaceID fid = world_.on_which_face(pid).value();
    Particle p = world_.get_particle(pid).second;

    const Real R = dom.shell().shape().radius() - p.radius();

    const auto disp    = this->draw_2D_displacement(fid, R);
    const auto old_pos = std::make_pair(p.position(), fid);
    const auto new_pos = ecell4::polygon::travel(world_.polygon(), old_pos, disp);

    std::tie(p.position(), fid) = new_pos;

    world_.update_particle_2D(pid, p, fid);

    return boost::container::static_vector<ParticleID, 2>{pid};
}

boost::container::static_vector<ParticleID, 2>
SingleCircularPropagator::reaction(const SingleCircularDomain& dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    // move particle and update world
    this->propagate(dom, dom.dt());

    // fetch propagated particle positions
    const auto pid   = dom.particle_id();
    const auto p     = world_.get_particle(pid).second;
    const auto fid   = world_.on_which_face(pid).value();

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
            return this->attempt_1to1_reaction(dom, pid, p, fid, rule);
        }
        case 2:
        {
            return this->attempt_1to2_reaction(dom, pid, p, fid, rule);
        }
        default:
        {
            throw_exception<NotSupported>("invalid number of products");
        }
    }
}

boost::container::static_vector<ParticleID, 2>
SingleCircularPropagator::attempt_1to1_reaction(const SingleCircularDomain& dom,
        const ParticleID& pid, const Particle& p, const FaceID& fid,
        const ReactionRule& rule)
{
    assert(rule.products().size() == 1);

    const auto species_new = rule.products().front();
    const auto molinfo     = world_.get_molecule_info(species_new);
    const Real radius_new  = molinfo.radius;
    const Real D_new       = molinfo.D;

    Particle newp(species_new, p.position(), radius_new, D_new);
    const auto pos = std::make_pair(p.position(), fid);

    // since a product could be larger than the reactant,
    // it may sticks out of the shell.
    if( ! this->is_inside_of_shell(dom, pos, newp.radius()))
    {
        // determine positions of particles in overlapping domains
        sim_.determine_positions_2D(pos,       radius_new, self_id_);
        sim_.determine_positions_3D(pos.first, radius_new, self_id_);

        // then check if there is any particle that overlaps
        if(world_.has_overlapping_particles_2D(pos,       radius_new, pid) ||
           world_.has_overlapping_particles_3D(pos.first, radius_new, pid))
        {
            this->rejected_move_count_ += 1;
            return boost::container::static_vector<ParticleID, 2>{pid};
        }
    }

    // looks okay.
    world_.update_particle_2D(pid, newp, fid);

    last_reactions_.emplace_back(rule, make_unimolecular_reaction_info(
                world_.t(), pid, p, pid, newp));
    return boost::container::static_vector<ParticleID, 2>{pid};
}

boost::container::static_vector<ParticleID, 2>
SingleCircularPropagator::attempt_1to2_reaction(const SingleCircularDomain& dom,
        const ParticleID& pid, const Particle& p, const FaceID& fid,
        const ReactionRule& rule)
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
        throw_exception<NotSupported>("ngfrd::SingleCircularPropagator attempts"
            " unbinding reaction but both particles are immovable => ",
            rule.as_string());
    }

    const auto boundary = world_.boundary();

    std::pair<Real3, FaceID> pos1_new(p.position(), fid);
    std::pair<Real3, FaceID> pos2_new(p.position(), fid);

    const Real separation_length = r12 * NGFRDSimulator::SAFETY_EXPAND;
    std::size_t separation_count = 1 + max_retry_count_;
    while(separation_count != 0)
    {
        --separation_count;

        // Here, the time when the reaciton happens is precisely determined.
        // So, unlike BD, we don't need to consider the reaction volume.
        // We can place the two particles in contact with each other.

        const Real3 ipv(random_circular_uniform(separation_length, fid));

        Real3 disp1(ipv * ( r1 / r12));
        Real3 disp2(ipv * (-r2 / r12));

        pos1_new = ecell4::polygon::travel(world_.polygon(), pos1_new, disp1);
        pos2_new = ecell4::polygon::travel(world_.polygon(), pos2_new, disp2);

        // numerical-error-phobia
        if(ecell4::polygon::distance(world_.polygon(), pos1_new, pos2_new) < r12)
        {
            continue;
        }

        // burst domains around the reactants.
        if(this->is_inside_of_shell(dom, pos1_new, r1))
        {
            sim_.determine_positions_2D(pos1_new, r1, self_id_);
            sim_.determine_positions_3D(pos1_new, r1, self_id_);
        }
        if(this->is_inside_of_shell(dom, pos2_new, r2))
        {
            sim_.determine_positions_2D(pos2_new, r2, self_id_);
            sim_.determine_positions_3D(pos2_new, r2, self_id_);
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
        return boost::container::static_vector<ParticleID, 2>{pid};
    }

    // burst domains around the reactants.
    if(this->is_inside_of_shell(dom, pos1_new, r1))
    {
        sim_.determine_positions_2D(pos1_new, r1, self_id_);
        sim_.determine_positions_3D(pos1_new, r1, self_id_);
    }
    if(this->is_inside_of_shell(dom, pos2_new, r2))
    {
        sim_.determine_positions_2D(pos2_new, r2, self_id_);
        sim_.determine_positions_3D(pos2_new, r2, self_id_);
    }

    // The reason why we use 2D and 3D, not XD is the following.
    // - 2D particles only interacts along the polygon surface.
    // - 3D distance between 2D particles does not affects.
    // - XD function checks overlap considering particles are spherical, so it
    //   counts 2D-2D 3D overlap.
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

    const auto result1 = world_.update_particle_2D(pid, p1_new, pos1_new.second);
    const auto result2 = world_.new_particle_2D   (     p2_new, pos2_new.second);

    assert(not result1); // should be already exist (it returns true if it's new)
    assert(result2.second); // should succeed

    const auto pid2 = result2.first.first;

    last_reactions_.emplace_back(rule, make_unbinding_reaction_info(world_.t(),
                pid, p, pid, p1_new, pid2, p2_new));

    return boost::container::static_vector<ParticleID, 2>{pid, pid2};
}

//
// Since proapgate is called only when burst is called, no reaction happens.
// Also, SingleCircularDomain is formed avoiding collision, no collision happens.
//
ParticleID SingleCircularPropagator::propagate(
        const SingleCircularDomain& dom, const Real dt)
{
    ECELL4_NGFRD_LOG_FUNCTION();

    const auto pid = dom.particle_id();

    FaceID fid = world_.on_which_face(pid).value();
    Particle p = world_.get_particle(pid).second;

    const auto& gf = dom.gf();
    const Real R = gf.drawR(rng_.uniform(0.0, 1.0), dt);

    const auto disp    = this->draw_2D_displacement(fid, R);
    const auto old_pos = std::make_pair(p.position(), fid);
    const auto new_pos = ecell4::polygon::travel(world_.polygon(), old_pos, disp);

    assert(this->is_inside_of_shell(dom, new_pos, p.radius()));

    std::tie(p.position(), fid) = new_pos;

    world_.update_particle_2D(pid, p, fid);
    return pid;
}

ReactionRule const& SingleCircularPropagator::determine_reaction_rule(
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

bool SingleCircularPropagator::is_inside_of_shell(const SingleCircularDomain& dom,
        const std::pair<Real3, FaceID>& pos, const Real radius)
{
    if(dom.shell().thickness() < radius)
    {
        return false;
    }
    const Real dist = ecell4::polygon::distance(world_.polygon(), pos,
            std::make_pair(dom.shell().position(), dom.shell().fid()));

    return dist + radius < dom.shell().shape().radius();
}

Real3 SingleCircularPropagator::draw_2D_displacement(const FaceID& fid, const Real len)
{
    constexpr Real pi = boost::math::constants::pi<Real>();

    const auto& tri  = world_.polygon().triangle_at(fid);
    const auto& disp = (tri.represent() / length(tri.represent())) * len;

    return rotate(rng_.uniform(0.0, 2 * pi), tri.normal(), disp);
}

} // ngfrd
} // ecell4
