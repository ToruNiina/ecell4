#include "python_api.hpp"

#include <ecell4/ngfrd/NGFRDFactory.hpp>
#include <ecell4/ngfrd/NGFRDSimulator.hpp>
#include <ecell4/ngfrd/NGFRDWorld.hpp>
#include <ecell4/core/STLFileIO.hpp>

#include "simulator.hpp"
#include "simulator_factory.hpp"
#include "world_interface.hpp"

namespace py = pybind11;

namespace ecell4
{

namespace python_api
{

static inline
void define_ngfrd_factory(py::module& m)
{
    using factory_type = ::ecell4::ngfrd::NGFRDFactory;

    py::class_<::ecell4::ngfrd::NGFRDFactory> factory(m, "NGFRDFactory");
    factory
        .def(py::init<const Integer3&, Real, Real, Real, Real, Integer>(),
             py::arg("matrix_sizes")              = factory_type::default_matrix_sizes(),
             py::arg("bd_dt_factor_3D")           = factory_type::default_bd_dt_factor_3D(),
             py::arg("bd_dt_factor_2D")           = factory_type::default_bd_dt_factor_2D(),
             py::arg("bd_reaction_length_factor") = factory_type::default_bd_reaction_length_factor(),
             py::arg("tree_margin")               = factory_type::default_tree_margin(),
             py::arg("max_retry_moves")           = factory_type::default_max_retry_count()
             )
        .def("rng",     &factory_type::rng)
        .def("polygon", (factory_type& (factory_type::*)(const Real3&, const std::vector<Triangle>&))&factory_type::polygon)
        .def("polygon", (factory_type& (factory_type::*)(const std::string&, const STLFormat))       &factory_type::polygon);

    define_factory_functions(factory);

    m.attr("Factory") = factory;
}

static inline
void define_ngfrd_simulator(py::module& m)
{
    py::class_<::ecell4::ngfrd::NGFRDSimulator, Simulator, PySimulator<::ecell4::ngfrd::NGFRDSimulator>,
        std::shared_ptr<::ecell4::ngfrd::NGFRDSimulator>> simulator(m, "NGFRDSimulator");
    simulator
        .def(py::init<
                const std::shared_ptr<::ecell4::ngfrd::NGFRDWorld>&,
                const std::shared_ptr<::ecell4::Model>&,
                const Real,
                const Real,
                const Real,
                const std::size_t
                >(),
                py::arg("w"),
                py::arg("m"),
                py::arg("bd_dt_factor_3D") = 1e-5,
                py::arg("bd_dt_factor_2D") = 1e-3,
                py::arg("reaction_length") = 1e-1,
                py::arg("max_retry_moves") = 1e-5)
        .def("last_reactions", &::ecell4::ngfrd::NGFRDSimulator::last_reactions)
        .def("set_t", &::ecell4::ngfrd::NGFRDSimulator::set_t);
    define_simulator_functions(simulator);

    m.attr("Simulator") = simulator;
}

static inline
void define_ngfrd_world(py::module& m)
{
    using ::ecell4::ngfrd::NGFRDWorld;
    py::class_<::ecell4::ngfrd::NGFRDWorld, WorldInterface, PyWorldImpl<::ecell4::ngfrd::NGFRDWorld>,
        std::shared_ptr<::ecell4::ngfrd::NGFRDWorld>> world(m, "NGFRDWorld");
    world
        .def(py::init<const Real3&, const Integer3&, const Real>(),
                py::arg("edge_lengths") = Real3(1.0, 1.0, 1.0),
                py::arg("matrix_sizes") = Integer3(1, 1, 1),
                py::arg("margin")       = Real(0.1)
                )
        .def(py::init<const Real3&, const Integer3&, const Real, std::shared_ptr<RandomNumberGenerator>>(),
                py::arg("edge_lengths") = Real3(1.0, 1.0, 1.0),
                py::arg("matrix_sizes") = Integer3(1, 1, 1),
                py::arg("margin")       = Real(0.1),
                py::arg("rng")
                )
        .def(py::init<const Real3&, const Integer3&, const Real, std::shared_ptr<RandomNumberGenerator>, std::shared_ptr<Polygon>>(),
                py::arg("edge_lengths") = Real3(1.0, 1.0, 1.0),
                py::arg("matrix_sizes") = Integer3(1, 1, 1),
                py::arg("margin")       = Real(0.1),
                py::arg("rng"),
                py::arg("polygon")
                )
        .def(py::init<const std::string&>(), py::arg("filename"))

        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool> (NGFRDWorld::*)(const Particle&))
            &NGFRDWorld::new_particle)
        .def("new_particle",
            (std::pair<std::pair<ParticleID, Particle>, bool> (NGFRDWorld::*)(const Species&, const Real3&))
            &NGFRDWorld::new_particle)
        .def("update_particle", &NGFRDWorld::update_particle)
        .def("remove_particle", &NGFRDWorld::remove_particle)
        .def("list_particles_within_radius",
            (std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
             (NGFRDWorld::*)(const Real3&, const Real) const)
            &NGFRDWorld::list_particles_within_radius_XD)
        .def("list_particles_within_radius",
            (std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
             (NGFRDWorld::*)(const Real3&, const Real, const ParticleID&) const)
            &NGFRDWorld::list_particles_within_radius_XD)
        .def("list_particles_within_radius",
            (std::vector<std::pair<std::pair<ParticleID, Particle>, Real>>
             (NGFRDWorld::*)(const Real3&, const Real, const ParticleID&, const ParticleID&) const)
            &NGFRDWorld::list_particles_within_radius_XD)
        .def("periodic_transpose", &NGFRDWorld::periodic_transpose)
        .def("apply_boundary", &NGFRDWorld::apply_boundary)
        .def("distance_sq", &NGFRDWorld::distance_sq)
        .def("distance", &NGFRDWorld::distance)
        .def("add_molecules",
            (void (NGFRDWorld::*)(const Species&, const Integer&)) &NGFRDWorld::add_molecules)
        .def("add_molecules",
            (void (NGFRDWorld::*)(const Species&, const Integer&, const std::shared_ptr<Shape>)) &NGFRDWorld::add_molecules)
        .def("remove_molecules", &NGFRDWorld::remove_molecules)
        .def("bind_to", &NGFRDWorld::bind_to)
        .def("rng", &NGFRDWorld::rng);

    m.attr("World") = world;
}

static inline
void define_reaction_info(py::module& m)
{
    using ::ecell4::ngfrd::ReactionInfo;
    using container_type = std::vector<std::pair<ParticleID, Particle>>;

    py::class_<ReactionInfo>(m, "ReactionInfo")
        .def(py::init<const Real, const container_type&, const container_type&>(),
                py::arg("t"), py::arg("reactants"), py::arg("products"))
        .def("t", &ReactionInfo::t)
        .def("reactants", &ReactionInfo::reactants)
        .def("products", &ReactionInfo::products)
        .def(py::pickle(
            [](const ReactionInfo& self)
            {
                return py::make_tuple(self.t(), self.reactants(), self.products());
            },
            [](py::tuple t)
            {
                if (t.size() != 3)
                    throw std::runtime_error("Invalid state");
                return ReactionInfo(
                    t[0].cast<Real>(),
                    t[1].cast<container_type>(),
                    t[2].cast<container_type>()
                );
            }
        ));
}

void setup_ngfrd_module(py::module& m)
{
    define_ngfrd_factory(m);
    define_ngfrd_simulator(m);
    define_ngfrd_world(m);
    define_reaction_info(m);
}

}

}
