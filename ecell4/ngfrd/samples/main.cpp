#include <ecell4/ngfrd/NGFRDWorld.hpp>
#include <ecell4/ngfrd/NGFRDSimulator.hpp>
#include <ecell4/core/NetworkModel.hpp>
#include <ecell4/core/Species.hpp>

void snapshot_output(std::ofstream& ofs,
        const std::shared_ptr<ecell4::ngfrd::NGFRDWorld>& world)
{
    ofs << world->num_particles() << '\n';
    ofs << world->t() << '\n';
    for(const auto& p : world->list_particles())
    {
        ofs << p.second.species().serial() << ' '
            << p.second.position()[0]      << ' '
            << p.second.position()[1]      << ' '
            << p.second.position()[2]      << '\n';
    }
    ofs << std::flush;
    return;
}

int main()
{
    const ecell4::Real     L(10.0);
    const ecell4::Real3    edge_lengths(L, L, L);
    const ecell4::Integer3 matrix_sizes(10, 10, 10);
    const ecell4::Real     volume(L * L * L);

    std::shared_ptr<ecell4::NetworkModel> model(new ecell4::NetworkModel());
    ecell4::Species sp1(std::string("A"), /*r = */ 1e-2, /*D = */ 1.0e-2);
    ecell4::Species sp2(std::string("B"), /*r = */ 1e-2, /*D = */ 1.0e-2);
    ecell4::Species sp3(std::string("C"), /*r = */ 2e-2, /*D = */ 0.5e-2);
    ecell4::Species sp4(std::string("D"), /*r = */ 1e-2, /*D = */ 1e-5);
    model->add_species_attribute(sp1);
    model->add_species_attribute(sp2);
    model->add_species_attribute(sp3);
    model->add_species_attribute(sp4);

    const ecell4::Real k_bind   = 0.1;
    const ecell4::Real k_unbind = 0.001;

    model->add_reaction_rule(ecell4::create_binding_reaction_rule  (sp1, sp2, sp3, k_bind));
    model->add_reaction_rule(ecell4::create_unbinding_reaction_rule(sp3, sp1, sp2, k_unbind));

    std::shared_ptr<ecell4::RandomNumberGenerator> rng =
        std::make_shared<ecell4::GSLRandomNumberGenerator>();
    rng->seed(123456ul);

    std::shared_ptr<ecell4::ngfrd::NGFRDWorld> world =
        std::make_shared<ecell4::ngfrd::NGFRDWorld>(
                edge_lengths, matrix_sizes, /*margin = */ 0.1, rng);
    world->add_molecules_3D(sp1, 20);
    world->add_molecules_3D(sp2, 20);

    world->add_molecules_2D(sp4, 10);

    std::ofstream traj("traj.xyz");
    snapshot_output(traj, world);

    ecell4::ngfrd::NGFRDSimulator sim(world, model);

    const ecell4::Real dt(1.0);
    for(std::size_t i=1; i<1000; ++i)
    {
        if(world->t() <= i * dt)
        {
            sim.initialize();
            while(sim.next_event_time() <= i * dt)
            {
                sim.step();
            }
            assert(sim.next_event_time() > i * dt);
            sim.diagnosis();
            sim.finalize();
        }
        snapshot_output(traj, world);
        std::cout << "t = " << world->t() << " / " << i * dt << std::endl;
    }

    return 0;
}
