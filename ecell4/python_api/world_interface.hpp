#ifndef ECELL4_PYTHON_API_WORLD_INTERFACE_HPP
#define ECELL4_PYTHON_API_WORLD_INTERFACE_HPP

#include <pybind11/pybind11.h>
#include <ecell4/core/WorldInterface.hpp>

namespace py = pybind11;

namespace ecell4
{

namespace python_api
{

    template<class Base = ecell4::WorldInterface>
    class PyWorldInterface: public Base
    {
    public:
        using Base::Base;

        const Real t() const override
        {
            PYBIND11_OVERLOAD(const Real, Base, t,);
        }
        void set_t(const Real& t) override
        {
            PYBIND11_OVERLOAD_PURE(void, Base, set_t, t);
        }

        void save(const std::string& filename) const override
        {
            PYBIND11_OVERLOAD_PURE(void, Base, save, filename);
        }
        void load(const std::string& filename) override
        {
            PYBIND11_OVERLOAD(void, Base, load, filename);
        }

        const Real volume() const override
        {
            PYBIND11_OVERLOAD(const Real, Base, volume,);
        }

        bool has_species(const Species& sp) const override
        {
            PYBIND11_OVERLOAD(bool, Base, has_species, sp);
        }

        std::vector<Species> list_species() const override
        {
            PYBIND11_OVERLOAD(std::vector<Species>, Base, list_species,);
        }

        Integer num_molecules(const Species& sp) const override
        {
            PYBIND11_OVERLOAD(Integer, Base, num_molecules, sp);
        }

        Integer num_molecules_exact(const Species& sp) const override
        {
            PYBIND11_OVERLOAD(Integer, Base, num_molecules_exact, sp);
        }

        Real get_value(const Species& sp) const override
        {
            PYBIND11_OVERLOAD(Real, Base, get_value, sp);
        }

        Real get_value_exact(const Species& sp) const override
        {
            PYBIND11_OVERLOAD(Real, Base, get_value_exact, sp);
        }

        const Real3& edge_lengths() const override
        {
            PYBIND11_OVERLOAD(const Real3&, Base, edge_lengths,);
        }

        Integer num_particles() const override
        {
            PYBIND11_OVERLOAD(Integer, Base, num_particles,);
        }

        Integer num_particles(const Species& sp) const override
        {
            PYBIND11_OVERLOAD(Integer, Base, num_particles, sp);
        }

        Integer num_particles_exact(const Species& sp) const override
        {
            PYBIND11_OVERLOAD(Integer, Base, num_particles_exact, sp);
        }

        bool has_particle(const ParticleID& pid) const override
        {
            PYBIND11_OVERLOAD(bool, Base, has_particle, pid);
        }

        using ParticleWithID = std::pair<ParticleID, Particle>;
        std::pair<ParticleID, Particle> get_particle(const ParticleID& pid) const override
        {
            PYBIND11_OVERLOAD(ParticleWithID, Base, get_particle, pid);
        }

        std::vector<std::pair<ParticleID, Particle>>
        list_particles() const override
        {
            PYBIND11_OVERLOAD(std::vector<ParticleWithID>, Base, list_particles,);
        }

        std::vector<std::pair<ParticleID, Particle>>
        list_particles(const Species& sp) const override
        {
            PYBIND11_OVERLOAD(std::vector<ParticleWithID>, Base, list_particles, sp);
        }

        std::vector<std::pair<ParticleID, Particle>>
        list_particles_exact(const Species& sp) const override
        {
            PYBIND11_OVERLOAD(std::vector<ParticleWithID>, Base, list_particles_exact, sp);
        }
    };

    template<class Base>
    class PyWorldImpl: public PyWorldInterface<Base>
    {
    public:
        using PyWorldInterface<Base>::PyWorldInterface;

        void set_t(const Real& t) override
        {
            PYBIND11_OVERLOAD(void, Base, set_t, t);
        }

        void save(const std::string& filename) const override
        {
            PYBIND11_OVERLOAD(void, Base, save, filename);
        }
    };

}

}
#endif /* ECELL4_PYTHON_API_WORLD_INTERFACE_HPP */
