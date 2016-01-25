#include <iterator>

#include <ecell4/core/exceptions.hpp>
#include <ecell4/core/Species.hpp>

#include "BDPropagator.hpp"
#include <ecell4/core/PlanarSurface.hpp>


#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>


namespace ecell4
{

namespace bd
{

template <typename T>
inline
signed int extract_sign(const T &number)
{
    return T(0) < number ? 1 : -1;
} 

boost::tuple<bool, Real3, Real3> refrection(const boost::shared_ptr<PlanarSurface> surface, const Real3& from, const Real3& displacement) 
{
    // return value:
    //  tuple(is_cross, intrusion_point, remaining_displacement_from_intrusion_point)
    Real3 temporary_destination(from + displacement);
    Real is_inside_from( surface->is_inside(from) );
    Real is_inside_dest( surface->is_inside(temporary_destination) );
    bool sign_of_dest_is_inside = 0. < is_inside_dest ? true : false;
    if (0 < is_inside_from * is_inside_dest) {
        return boost::make_tuple(false, from , displacement);
    }
    else if (0 < is_inside_from && is_inside_dest < 0) {
        // Inside -> refrection -> Inside
        Real distance_from_surface(std::abs(is_inside_dest));
        Real3 new_pos = temporary_destination + multiply(surface->normal(), (-2.0) * distance_from_surface);
        Real ratio( std::abs(is_inside_from) / (std::abs(is_inside_from) + std::abs(is_inside_dest)) );
        Real3 intrusion_point( from + multiply(displacement, ratio) );
        return boost::make_tuple(true, intrusion_point, new_pos - intrusion_point);
    } 
    else if (is_inside_from < 0 && 0 < is_inside_dest) {
        Real distance_from_surface(std::abs(is_inside_dest));
        Real3 new_pos = temporary_destination + multiply(surface->normal(), (2.0) * distance_from_surface);
        Real ratio( std::abs(is_inside_from) / (std::abs(is_inside_from) + std::abs(is_inside_dest)) );
        Real3 intrusion_point( from + multiply(displacement, ratio) );
        return boost::make_tuple(true, intrusion_point, new_pos - intrusion_point);
    }
}

bool BDPropagator::operator()()
{
    if (queue_.empty())
    {
        return false;
    }

    const ParticleID pid(queue_.back().first);
    queue_.pop_back();
    Particle particle(world_.get_particle(pid).second);

    if (attempt_reaction(pid, particle))
    {
        return true;
    }

    const Real D(particle.D());
    if (D == 0)
    {
        return true;
    }

    const std::vector<boost::shared_ptr<PlanarSurface> > surface_vector = world_.get_surface_container();

    Real3 from = particle.position();
    Real3 displacement(draw_displacement(particle));
    std::size_t bound_surface = -1;
    std::vector<signed int> save_isinside(surface_vector.size());
    for(std::size_t i = 0; i != surface_vector.size(); i++) {
        // the inside or outside status must not be changed by definition.
        save_isinside[i] = extract_sign( surface_vector[i]->is_inside(from) );
    }
    bool refrection_occurance = false;
    do {
        refrection_occurance = false;
        boost::tuple<bool, Real3, Real3> nearest;
        std::size_t nearest_surface = -1;
        for(std::size_t i = 0; i != surface_vector.size(); i++) {
            if (i != bound_surface) {
                boost::tuple<bool, Real3, Real3> t = refrection(surface_vector[i], from, displacement);
                if (t.get<0>() == true) {
                    if (false == refrection_occurance) {
                        refrection_occurance = true;
                        nearest = t;
                        nearest_surface = i;
                    } else {
                        // compare the intrusion point
                        if (length(from - t.get<1>() ) < length(from - nearest.get<1>() )) {
                            nearest = t;
                            nearest_surface = i;
                        }
                    }
                }
            }
        }
        // update or escape from the loop
        if (refrection_occurance == false) {
            break;
        } else {
            from = nearest.get<1>() ;
            displacement = nearest.get<2>();
            bound_surface = nearest_surface;
            // For Debugging
            std::cout << "bound surface ( " << bound_surface << ") at " << from << std::endl;
        }
    } while(true);
    Real3 newpos(world_.apply_boundary(from + displacement));

    // Check for debugging
    for(std::size_t i = 0; i != surface_vector.size(); i++) {
        // the inside or outside status must not be changed before and after moving.
        signed int is_inside = extract_sign( surface_vector[i]->is_inside(from + displacement) );
        if (save_isinside[i] != is_inside) {
            throw IllegalState("Particle moved to the opposite side of the surface");
        }
    }

    /*
    const Real3 newpos(
        world_.apply_boundary(
            particle.position() + draw_displacement(particle)));
    */
    Particle particle_to_update(
        particle.species(), newpos, particle.radius(), particle.D());
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
        overlapped(world_.list_particles_within_radius(
                       newpos, particle.radius(), pid));

    switch (overlapped.size())
    {
    case 0:
        world_.update_particle(pid, particle_to_update);
        return true;
    case 1:
        {
            std::pair<ParticleID, Particle> closest(
                (*(overlapped.begin())).first);
            if (attempt_reaction(
                    pid, particle_to_update, closest.first, closest.second))
            {
                return true;
            }
        }
        return true;
    default:
        return true;
    }
}

bool BDPropagator::attempt_reaction(
    const ParticleID& pid, const Particle& particle)
{
    std::vector<ReactionRule> reaction_rules(
        model_.query_reaction_rules(particle.species()));
    if (reaction_rules.size() == 0)
    {
        return false;
    }

    const Real rnd(rng().uniform(0, 1));
    Real prob(0);
    for (std::vector<ReactionRule>::const_iterator i(reaction_rules.begin());
         i != reaction_rules.end(); ++i)
    {
        const ReactionRule& rr(*i);
        prob += rr.k() * dt();
        if (prob > rnd)
        {
            const ReactionRule::product_container_type& products(rr.products());
            reaction_info_type ri(world_.t() + dt_, reaction_info_type::container_type(1, std::make_pair(pid, particle)), reaction_info_type::container_type());

            switch (products.size())
            {
            case 0:
                remove_particle(pid);
                last_reactions_.push_back(std::make_pair(rr, ri));
                break;
            case 1:
                {
                    const Species species_new(
                        model_.apply_species_attributes(*(products.begin())));
                    const BDWorld::molecule_info_type
                        info(world_.get_molecule_info(species_new));
                    const Real radius_new(info.radius);
                    const Real D_new(info.D);

                    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
                        overlapped(world_.list_particles_within_radius(
                                       particle.position(), radius_new, pid));
                    if (overlapped.size() > 0)
                    {
                        // throw NoSpace("");
                        return false;
                    }

                    Particle particle_to_update(
                        species_new, particle.position(), radius_new, D_new);
                    world_.update_particle(pid, particle_to_update);

                    ri.add_product(std::make_pair(pid, particle_to_update));
                    last_reactions_.push_back(std::make_pair(rr, ri));
                }
                break;
            case 2:
                {
                    ReactionRule::product_container_type::const_iterator
                        it(products.begin());
                    const Species species_new1(
                        model_.apply_species_attributes(*it));
                    const Species species_new2(
                        model_.apply_species_attributes(*(++it)));

                    const BDWorld::molecule_info_type
                        info1(world_.get_molecule_info(species_new1)),
                        info2(world_.get_molecule_info(species_new2));
                    const Real radius1(info1.radius),
                        radius2(info2.radius);
                    const Real D1(info1.D), D2(info2.D);

                    const Real D12(D1 + D2);
                    const Real r12(radius1 + radius2);
                    Real3 newpos1, newpos2;
                    Integer i(max_retry_count_);
                    while (true)
                    {
                        if (--i < 0)
                        {
                            // throw NoSpace("")
                            return false;
                        }

                        const Real3 ipv(draw_ipv(r12, dt(), D12));

                        newpos1 = world_.apply_boundary(
                            particle.position() + ipv * (D1 / D12));
                        newpos2 = world_.apply_boundary(
                            particle.position() - ipv * (D2 / D12));
                        std::vector<
                            std::pair<std::pair<ParticleID, Particle>, Real> >
                            overlapped1(world_.list_particles_within_radius(
                                            newpos1, radius1, pid));
                        std::vector<
                            std::pair<std::pair<ParticleID, Particle>, Real> >
                            overlapped2(world_.list_particles_within_radius(
                                            newpos2, radius2, pid));
                        if (overlapped1.size() == 0 && overlapped2.size() == 0)
                        {
                            break;
                        }
                    }

                    Particle particle_to_update1(
                        species_new1, newpos1, radius1, D1);
                    Particle particle_to_update2(
                        species_new2, newpos2, radius2, D2);
                    world_.update_particle(pid, particle_to_update1);
                    std::pair<std::pair<ParticleID, Particle>, bool> retval = world_.new_particle(particle_to_update2);

                    ri.add_product(std::make_pair(pid, particle_to_update1));
                    ri.add_product(retval.first);
                    last_reactions_.push_back(std::make_pair(rr, ri));
                }
                break;
            default:
                throw NotImplemented(
                    "more than two products are not allowed");
                break;
            }
            return true;
        }
    }

    return false;
}

bool BDPropagator::attempt_reaction(
    const ParticleID& pid1, const Particle& particle1,
    const ParticleID& pid2, const Particle& particle2)
{
    std::vector<ReactionRule> reaction_rules(
        model_.query_reaction_rules(
            particle1.species(), particle2.species()));
    if (reaction_rules.size() == 0)
    {
        return false;
    }

    const Real D1(particle1.D()), D2(particle2.D());
    const Real r12(particle1.radius() + particle2.radius());
    const Real rnd(rng().uniform(0, 1));
    Real prob(0);

    for (std::vector<ReactionRule>::const_iterator i(reaction_rules.begin());
         i != reaction_rules.end(); ++i)
    {
        const ReactionRule& rr(*i);
        prob += rr.k() * dt() / (
            (Igbd_3d(r12, dt(), D1) + Igbd_3d(r12, dt(), D2)) * 4 * M_PI);

        if (prob >= 1)
        {
            // throw std::runtime_error(
            //     "the total reaction probability exceeds 1."
            //     " the step interval is too long");
            std::cerr <<
                "the total reaction probability exceeds 1."
                " the step interval is too long" << std::endl;
        }
        if (prob > rnd)
        {
            const ReactionRule::product_container_type& products(rr.products());
            reaction_info_type ri(world_.t() + dt_, reaction_info_type::container_type(1, std::make_pair(pid1, particle1)), reaction_info_type::container_type());
            ri.add_reactant(std::make_pair(pid2, particle2));

            switch (products.size())
            {
            case 0:
                remove_particle(pid1);
                remove_particle(pid2);

                last_reactions_.push_back(std::make_pair(rr, ri));
                break;
            case 1:
                {
                    const Species sp(*(products.begin()));
                    const BDWorld::molecule_info_type
                        info(world_.get_molecule_info(sp));
                    const Real radius_new(info.radius);
                    const Real D_new(info.D);

                    const Real3 pos1(particle1.position());
                    const Real3 pos2(
                        world_.periodic_transpose(particle2.position(), pos1));
                    const Real D1(particle1.D()), D2(particle2.D());
                    const Real D12(D1 + D2);
                    const Real3 newpos(
                        world_.apply_boundary((pos1 * D2 + pos2 * D1) / D12));

                    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
                        overlapped(world_.list_particles_within_radius(
                                       newpos, radius_new, pid1, pid2));
                    if (overlapped.size() > 0)
                    {
                        // throw NoSpace("");
                        return false;
                    }

                    const Particle particle_to_update(
                        sp, newpos, radius_new, D_new);
                    remove_particle(pid2);
                    // world_.update_particle(pid1, particle_to_update);
                    remove_particle(pid1);
                    std::pair<std::pair<ParticleID, Particle>, bool> retval = world_.new_particle(particle_to_update);

                    ri.add_product(retval.first);
                    last_reactions_.push_back(std::make_pair(rr, ri));
                }
                break;
            default:
                throw NotImplemented(
                    "more than one product is not allowed");
                break;
            }
            return true;
        }
    }

    return false;
}

void BDPropagator::remove_particle(const ParticleID& pid)
{
    world_.remove_particle(pid);
    particle_finder cmp(pid);
    std::vector<std::pair<ParticleID, Particle> >::iterator
        i(std::find_if(queue_.begin(), queue_.end(), cmp));
    if (i != queue_.end())
    {
        queue_.erase(i);
    }
}

} // bd

} // ecell4
