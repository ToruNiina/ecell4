#ifndef ECELL4_PARTICLE_SPACE_HPP
#define ECELL4_PARTICLE_SPACE_HPP

#include <cmath>

#include "get_mapper_mf.hpp"
#include "types.hpp"
#include "functions.hpp"
#include "exceptions.hpp"
#include "Real3.hpp"
#include "Particle.hpp"
#include "Species.hpp"
#include "BoundaryCondition.hpp"
#include "Context.hpp"
#include "comparators.hpp"

#ifdef WITH_HDF5
#include "ParticleSpaceHDF5Writer.hpp"
#endif

namespace ecell4
{

class ParticleSpace
    // : public Space
{
public:

    typedef std::vector<std::pair<ParticleID, Particle> >
    particle_container_type;

public:

    ParticleSpace()
        : t_(0.0)
    {
        ;
    }

    virtual ~ParticleSpace()
    {
        ; // do nothing
    }

    // SpaceTraits

    const Real t() const
    {
        return t_;
    }

    void set_t(const Real& t)
    {
        if (t < 0.0)
        {
            throw std::invalid_argument("the time must be positive.");
        }
        t_ = t;
    }

    // ParticleSpaceTraits

    /**
     * get the axes lengths of a cuboidal region.
     * this function is a part of the trait of ParticleSpace.
     * @return edge lengths Real3
     */
    virtual const Real3& edge_lengths() const
    {
        throw NotImplemented("edge_lengths() not implemented");
    }

    virtual Integer num_molecules(const Species& sp) const
    {
        return num_molecules_exact(sp);
    }

    virtual Integer num_molecules_exact(const Species& sp) const
    {
        throw NotImplemented("num_molecules_exact(const Species&) not implemented");
    }

    /**
     * get the number of particles.
     * this function is a part of the trait of ParticleSpace.
     * @return a number of particles Integer
     */
    virtual Integer num_particles() const
    {
        throw NotImplemented("num_particles() not implemented");
    }

    /**
     * get the number of particles.
     * this function is a part of the trait of ParticleSpace.
     * @param sp a species
     * @return a number of particles Integer
     */
    virtual Integer num_particles(const Species& sp) const
    {
        return num_particles_exact(sp);
    }

    virtual Integer num_particles_exact(const Species& sp) const
    {
        throw NotImplemented("num_particles_exact(const Species&) not implemented");
    }

    /**
     * get all particles.
     * this function is a part of the trait of ParticleSpace.
     * @return a list of particles
     */
    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles() const
    {
        throw NotImplemented("list_particles() not implemented");
    }

    /**
     * get particles.
     * this function is a part of the trait of ParticleSpace.
     * @param sp a species
     * @return a list of particles
     */
    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& sp) const
    {
        return list_particles_exact(sp);
    }

    virtual std::vector<std::pair<ParticleID, Particle> >
    list_particles_exact(const Species& sp) const
    {
        throw NotImplemented("list_particles_exact(const Species&) not implemented");
    }

    /**
     * check if the particle exists.
     * this function is a part of the trait of ParticleSpace.
     * @param pid an ID for the particle
     * @return if the particle exists or not bool
     */
    virtual bool has_particle(const ParticleID& pid) const
    {
        throw NotImplemented("has_particle(const ParticleID&) not implemented.");
    }

    virtual std::vector<Species> list_species() const
    {
        const particle_container_type& pcont(particles());
        std::vector<Species> retval;
        for (particle_container_type::const_iterator i(pcont.begin());
            i != pcont.end(); ++i)
        {
            const Species& sp((*i).second.species());
            if (std::find(retval.begin(), retval.end(), sp)
                == retval.end())
            {
                retval.push_back(sp);
            }
        }
        return retval;
    }

#ifdef WITH_HDF5
    virtual void save_hdf5(H5::Group* root) const = 0;
    virtual void load_hdf5(const H5::Group& root) = 0;
#endif

    // ParticleSpace member functions

    /**
     * update a particle specified by its ID.
     * if the particle does not exist, create a new particle.
     * this function is a member of ParticleSpace
     * @param pid ParticleID
     * @param p Particle
     * @return if the particle does not exist or not bool
     */
    virtual bool update_particle(const ParticleID& pid, const Particle& p) = 0;

    /**
     * get a particle specified with an ID.
     * this function is a member of ParticleSpace
     * @param pid ParticleID
     * @return a pair of ParticleID and Particle
     */
    virtual std::pair<ParticleID, Particle>
    get_particle(const ParticleID& pid) const = 0;

    /**
     * remove a particle
     * this function is a member of ParticleSpace
     * @param pid ParticleID
     */
    virtual void remove_particle(const ParticleID& pid) = 0;

    /**
     * get particles within a spherical region.
     * this function is a part of the trait of ParticleSpace.
     * @param pos a center position of the sphere
     * @param radius a radius of the sphere
     * @return a list of particles
     */
    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Real3& pos, const Real& radius) const = 0;

    /**
     * get particles within a spherical region except for ignore(s).
     * this function is a part of the trait of ParticleSpace.
     * @param pos a center position of the sphere
     * @param radius a radius of the sphere
     * @param ignore an ignored ID
     * @return a list of particles
     */
    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Real3& pos, const Real& radius,
        const ParticleID& ignore) const = 0;

    /**
     * get particles within a spherical region except for ignore(s).
     * this function is a part of the trait of ParticleSpace.
     * @param pos a center position of the sphere
     * @param radius a radius of the sphere
     * @param ignore1 an ignored ID
     * @param ignore2 an ignored ID
     * @return a list of particles
     */
    virtual std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Real3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const = 0;

    /**
     * transpose a position based on the periodic boundary condition.
     * this function is a part of the trait of ParticleSpace.
     * @param pos1 a target position
     * @param pos2 a reference position
     * @return a transposed position Real3
     */
    virtual Real3 periodic_transpose(
        const Real3& pos1, const Real3& pos2) const = 0;

    /**
     * transpose a position based on the periodic boundary condition.
     * if the position is in the region, returns the original position.
     * this function is a part of the trait of ParticleSpace.
     * @param pos a target position
     * @return a transposed position Real3
     */
    virtual Real3 apply_boundary(const Real3& pos) const = 0;

    /**
     * calculate a square of the distance between positions
     * this function is a part of the trait of ParticleSpace.
     * @param pos1
     * @param pos2
     * @return a square of the distance
     */
    Real distance_sq(
        const Real3& pos1, const Real3& pos2) const
    {
        // Real retval(0);
        // const Real3& edges(edge_lengths());
        // for (Real3::size_type dim(0); dim < 3; ++dim)
        // {
        //     const Real edge_length(edges[dim]);
        //     const Real diff(pos2[dim] - pos1[dim]), half(edge_length * 0.5);

        //     if (diff > half)
        //     {
        //         retval += pow_2(diff - edge_length);
        //     }
        //     else if (diff < -half)
        //     {
        //         retval += pow_2(diff + edge_length);
        //     }
        //     else
        //     {
        //         retval += pow_2(diff);
        //     }
        // }
        // return retval;

        return length_sq(subtract(pos1, periodic_transpose(pos2, pos1)));
    }

    /**
     * calculate the distance between positions
     * this function is a part of the trait of ParticleSpace.
     * @param pos1
     * @param pos2
     * @return the distance
     */
    inline Real distance(const Real3& pos1, const Real3& pos2) const
    {
        return std::sqrt(distance_sq(pos1, pos2));
    }

    // Optional members

    virtual const particle_container_type& particles() const = 0;

    virtual Real get_value(const Species& sp) const
    {
        return static_cast<Real>(num_molecules(sp));
    }

    virtual Real get_value_exact(const Species& sp) const
    {
        return static_cast<Real>(num_molecules_exact(sp));
    }

    virtual const Real volume() const
    {
        const Real3& size(edge_lengths());
        return size[0] * size[1] * size[2];
    }

protected:

    Real t_;
};

template<typename BoundaryCondition>
class ParticleSpaceVectorImpl
    : public ParticleSpace
{
public:

    typedef ParticleSpace base_type;
    typedef ParticleSpace::particle_container_type particle_container_type;
    typedef BoundaryCondition boundary_condition_type;

protected:

    typedef utils::get_mapper_mf<
        ParticleID, particle_container_type::size_type>::type particle_map_type;

public:

    ParticleSpaceVectorImpl(const Real3& edge_lengths)
    {
        reset(edge_lengths);
    }

    // ParticleSpaceTraits

    const Real3& edge_lengths() const
    {
        return boundary_.edge_lengths();
    }

    Real3 periodic_transpose(
        const Real3& pos1, const Real3& pos2) const final
    {
        return boundary_.periodic_transpose(pos1, pos2);
    }

    Real3 apply_boundary(const Real3& pos) const final
    {
        return boundary_.apply_boundary(pos);
    }

    Integer num_particles() const;
    Integer num_particles(const Species& sp) const;
    Integer num_particles_exact(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> > list_particles() const;
    std::vector<std::pair<ParticleID, Particle> >
    list_particles(const Species& sp) const;
    std::vector<std::pair<ParticleID, Particle> >
    list_particles_exact(const Species& sp) const;
    bool has_particle(const ParticleID& pid) const;

    // ParticleSpace member functions

    bool update_particle(const ParticleID& pid, const Particle& p);
    std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const;
    void remove_particle(const ParticleID& pid);

    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Real3& pos, const Real& radius) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Real3& pos, const Real& radius,
        const ParticleID& ignore) const;
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
    list_particles_within_radius(
        const Real3& pos, const Real& radius,
        const ParticleID& ignore1, const ParticleID& ignore2) const;

    // CompartmentSpaceTraits

    Integer num_molecules(const Species& sp) const;
    Integer num_molecules_exact(const Species& sp) const;

    // Optional members

    const particle_container_type& particles() const
    {
        return particles_;
    }

    virtual void save(const std::string& filename) const
    {
        throw NotSupported(
            "save(const std::string) is not supported by this space class");
    }

#ifdef WITH_HDF5
    void save_hdf5(H5::Group* root) const
    {
        save_particle_space(*this, root);
    }

    void load_hdf5(const H5::Group& root)
    {
        load_particle_space(root, this);
    }
#endif

    void reset(const Real3& edge_lengths);

protected:

    boundary_condition_type boundary_;
    particle_container_type particles_;
    particle_map_type index_map_;
};

extern template class ParticleSpaceVectorImpl<PeriodicBoundary>;
extern template class ParticleSpaceVectorImpl<UnlimitedBoundary>;

template<typename BoundaryCondition>
Integer ParticleSpaceVectorImpl<BoundaryCondition>::num_particles() const
{
    return static_cast<Integer>(particles_.size());
}

template<typename BoundaryCondition>
Integer ParticleSpaceVectorImpl<BoundaryCondition>::num_particles(const Species& sp) const
{
    return static_cast<Integer>(list_particles(sp).size());
}

template<typename BoundaryCondition>
Integer ParticleSpaceVectorImpl<BoundaryCondition>::num_molecules(const Species& sp) const
{
    Integer retval(0);
    SpeciesExpressionMatcher sexp(sp);
    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        retval += sexp.count((*i).second.species());
    }
    return retval;
}

template<typename BoundaryCondition>
Integer ParticleSpaceVectorImpl<BoundaryCondition>::num_molecules_exact(const Species& sp) const
{
    return num_particles_exact(sp);
}

template<typename BoundaryCondition>
Integer ParticleSpaceVectorImpl<BoundaryCondition>::num_particles_exact(const Species& sp) const
{
    return static_cast<Integer>(list_particles_exact(sp).size());
}

template<typename BoundaryCondition>
bool ParticleSpaceVectorImpl<BoundaryCondition>::has_particle(const ParticleID& pid) const
{
    particle_map_type::const_iterator i(index_map_.find(pid));
    return (i != index_map_.end());
}

/**
 * update or add a particle.
 * @return true if adding a new particle
 */
template<typename BoundaryCondition>
bool ParticleSpaceVectorImpl<BoundaryCondition>::update_particle(
    const ParticleID& pid, const Particle& p)
{
    particle_map_type::const_iterator i(index_map_.find(pid));
    if (i == index_map_.end())
    {
        particle_container_type::size_type idx(particles_.size());
        index_map_[pid] = idx;
        particles_.push_back(std::make_pair(pid, p));
        return true;
    }
    else
    {
        particles_[(*i).second] = std::make_pair(pid, p);
        return false;
    }
}

template<typename BoundaryCondition>
void ParticleSpaceVectorImpl<BoundaryCondition>::remove_particle(const ParticleID& pid)
{
    particle_map_type::const_iterator i(index_map_.find(pid));
    if (i == index_map_.end())
    {
        throw NotFound("particle not found");
    }

    particle_map_type::mapped_type
        idx((*i).second),last_idx(particles_.size() - 1);
    if (idx != last_idx)
    {
        const std::pair<ParticleID, Particle>& last(particles_[last_idx]);
        particles_[idx] = last;
        index_map_[last.first] = idx;
    }

    particles_.pop_back();
    index_map_.erase((*i).first);
}

template<typename BoundaryCondition>
std::pair<ParticleID, Particle> ParticleSpaceVectorImpl<BoundaryCondition>::get_particle(
    const ParticleID& pid) const
{
    particle_map_type::const_iterator i(index_map_.find(pid));
    if (i == index_map_.end())
    {
        throw NotFound("particle not found");
    }

    return particles_[(*i).second];
}

template<typename BoundaryCondition>
std::vector<std::pair<ParticleID, Particle> >
ParticleSpaceVectorImpl<BoundaryCondition>::list_particles() const
{
    return particles_;
}

template<typename BoundaryCondition>
std::vector<std::pair<ParticleID, Particle> >
ParticleSpaceVectorImpl<BoundaryCondition>::list_particles(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    SpeciesExpressionMatcher sexp(sp);

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        if (sexp.match((*i).second.species()))
        {
            retval.push_back(*i);
        }
    }
    return retval;
}

template<typename BoundaryCondition>
std::vector<std::pair<ParticleID, Particle> >
ParticleSpaceVectorImpl<BoundaryCondition>::list_particles_exact(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        if ((*i).second.species() == sp)
        {
            retval.push_back(*i);
        }
    }

    return retval;
}

template<typename BoundaryCondition>
std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl<BoundaryCondition>::list_particles_within_radius(
    const Real3& pos, const Real& radius) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        const Real dist(distance((*i).second.position(), pos) - (*i).second.radius());
        if (dist <= radius)
        {
            retval.push_back(std::make_pair(*i, dist));
        }
    }

    std::sort(retval.begin(), retval.end(),
        utils::pair_second_element_comparator<std::pair<ParticleID, Particle>, Real>());
    return retval;
}

template<typename BoundaryCondition>
std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl<BoundaryCondition>::list_particles_within_radius(
    const Real3& pos, const Real& radius, const ParticleID& ignore) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        const Real dist(distance((*i).second.position(), pos) - (*i).second.radius());
        if (dist <= radius)
        {
            if ((*i).first != ignore)
            {
                retval.push_back(std::make_pair(*i, dist));
            }
        }
    }

    std::sort(retval.begin(), retval.end(),
        utils::pair_second_element_comparator<std::pair<ParticleID, Particle>, Real>());
    return retval;
}

template<typename BoundaryCondition>
std::vector<std::pair<std::pair<ParticleID, Particle>, Real> >
ParticleSpaceVectorImpl<BoundaryCondition>::list_particles_within_radius(
    const Real3& pos, const Real& radius,
    const ParticleID& ignore1, const ParticleID& ignore2) const
{
    std::vector<std::pair<std::pair<ParticleID, Particle>, Real> > retval;

    for (particle_container_type::const_iterator i(particles_.begin());
         i != particles_.end(); ++i)
    {
        const Real dist(distance((*i).second.position(), pos) - (*i).second.radius());
        if (dist <= radius)
        {
            if ((*i).first != ignore1 && (*i).first != ignore2)
            {
                retval.push_back(std::make_pair(*i, dist));
            }
        }
    }

    std::sort(retval.begin(), retval.end(),
        utils::pair_second_element_comparator<std::pair<ParticleID, Particle>, Real>());
    return retval;
}

template<typename BoundaryCondition>
void ParticleSpaceVectorImpl<BoundaryCondition>::reset(const Real3& edge_lengths)
{
    base_type::t_ = 0.0;
    particles_.clear();
    index_map_.clear();

    for (Real3::size_type dim(0); dim < 3; ++dim)
    {
        if (edge_lengths[dim] <= 0)
        {
            throw std::invalid_argument("the edge length must be positive.");
        }
    }

    this->boundary_.update(edge_lengths);
}

} // ecell4
#endif /* ECELL4_PARTICLE_SPACE_HPP */
