// vim: foldmethod=marker
// copied from Sakamoto EPDP Sample, mymapk
// remove reaction, change output format to xyz, and add polygon structure.

#include <stdexcept>
#include <vector>
#include <string>
#include <numeric>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <cstdlib>
#include <gsl/gsl_roots.h>

#include <boost/format.hpp>

// // epdp headers
// #include <ecell4/epdp/config.h>
// #include <ecell4/epdp/utils/range.hpp>
// #include <ecell4/epdp/World.hpp>
// //#include <ecell4/epdp/ParticleModel.hpp>
// //#include <ecell4/epdp/SpeciesType.hpp>
// //#include <ecell4/epdp/SpeciesTypeID.hpp>
// //#include <ecell4/epdp/CuboidalRegion.hpp>
// //#include <ecell4/epdp/NetworkRules.hpp>
// //#include <ecell4/epdp/ReactionRule.hpp>
// #include <ecell4/epdp/EGFRDSimulator.hpp>
// //#include <ecell4/epdp/NetworkRulesAdapter.hpp>
// //#include <ecell4/epdp/GSLRandomNumberGenerator.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Species.hpp>
#include <ecell4/core/ReactionRule.hpp>

#include <ecell4/egfrd/egfrd.hpp>


// typedef double Real;

int main(int argc, char **argv)
{
    // Traits typedefs
    // {{{
    // typedef ::World< ::CyclicWorldTraits<Real> > world_type;
    // typedef EGFRDSimulator< ::EGFRDSimulatorTraitsBase<world_type> >
    //     simulator_type;
    typedef ecell4::egfrd::EGFRDWorld world_type;
    typedef ecell4::egfrd::EGFRDSimulator simulator_type;
    typedef simulator_type::multi_type multi_type;
    // }}}

    // Constants
    // {{{
    const ecell4::Real L(1e2);
    const ecell4::Real3 edge_lengths(L, L, L);
    const ecell4::Integer3 matrix_sizes(3, 3, 3);
    const ecell4::Real volume(L * L * L);
    const ecell4::Integer N(60);
    const ecell4::Real kd(0.1), U(0.5);
    const ecell4::Real ka(kd * volume * (1 - U) / (U * U * N));
    const ecell4::Real k2(ka), k1(kd);
    const ecell4::Integer dissociation_retry_moves(3);
    // }}}

    boost::shared_ptr<ecell4::NetworkModel>
        model(new ecell4::NetworkModel());

    // add ::SpeciesType to ::ParticleModel
    // {{{
    ecell4::Species sp1(
        std::string("C"), std::string("2.5"), std::string("1e3"));
    model->add_species_attribute(sp1);
    // }}}

    // ReactionRules
    // {{{
    // Nothing!
    // }}}

    // Random Number Generator (Instanciate and Initialize)
    // {{{
    // boost::shared_ptr<ecell4::GSLRandomNumberGenerator>
    boost::shared_ptr<ecell4::RandomNumberGenerator>
        rng(new ecell4::GSLRandomNumberGenerator());
    rng->seed((unsigned long int)0);
    // rng->seed(time(NULL));
    // }}}

    // World Definition
    // {{{
    boost::shared_ptr<world_type>
        world(new world_type(edge_lengths, matrix_sizes, rng));
    world->bind_to(model);
    // }}}

    // Thorow particles into world at random
    // {{{
    world->add_molecules(ecell4::Species("C"), N);

    typedef std::vector<std::pair<ecell4::ParticleID, ecell4::Particle> >
        particle_id_pair_list;
    const particle_id_pair_list particles(world->list_particles());
    for (particle_id_pair_list::const_iterator i(particles.begin());
        i != particles.end(); ++i)
    {
        const ecell4::Real3 pos((*i).second.position());
    }
    // }}}

    // Logger Settings
    // {{{
    boost::shared_ptr< ::LoggerManager> logger_mng(
        new ::LoggerManager("dummy", ::Logger::L_WARNING));
    ::LoggerManager::register_logger_manager(
        "ecell.EGFRDSimulator", logger_mng);
    // }}}

    // EGFRDSimulator instance generated
    // {{{
    boost::shared_ptr<simulator_type> sim(
        new simulator_type(world, model, dissociation_retry_moves));
    // sim->paranoiac() = true;
    sim->initialize();
    // }}}

    // Simulation Executed
    // {{{
    ecell4::Real next_time(0.0), dt(0.2);
    for (int i(0); i < 100; i++)
    {
        next_time += dt;
        while (sim->step(next_time)){}

        ecell4::Integer num_sp1 = world->num_molecules_exact(sp1);
        std::cout << num_sp1 << std::endl;
        std::cout << "now t = " << sim->t() << std::endl;
        const std::vector<
            std::pair<world_type::particle_id_type, world_type::particle_type>
            > particles = world->list_particles();
        for(std::size_t i=0; i< particles.size(); ++i)
        {
            std::cout << particles.at(i).second.sid() << " "
                << particles.at(i).second.position()[0] << " "
                << particles.at(i).second.position()[1] << " "
                << particles.at(i).second.position()[2] << std::endl;
        }
    }
    // }}}

    return 0;
}