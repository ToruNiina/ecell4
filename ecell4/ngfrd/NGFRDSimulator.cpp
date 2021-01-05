#include <ecell4/ngfrd/NGFRDSimulator.hpp>
#include <ecell4/ngfrd/SingleSphericalPropagator.hpp>
#include <ecell4/ngfrd/SingleCircularPropagator.hpp>

namespace ecell4
{
namespace ngfrd
{
constexpr Real NGFRDSimulator::SAFETY_EXPAND;
constexpr Real NGFRDSimulator::SAFETY_SHRINK;
constexpr Real NGFRDSimulator::SINGLE_SHELL_FACTOR;
constexpr Real NGFRDSimulator::DEFAULT_DT_FACTOR;
constexpr Real NGFRDSimulator::CUTOFF_FACTOR;

constexpr std::size_t NGFRDSimulator::SINGLE_SPHERICAL_MAX_RETRY;
constexpr std::size_t NGFRDSimulator::SINGLE_CONICAL_MAX_RETRY;

boost::optional<boost::container::small_vector<DomainID, 4>>
NGFRDSimulator::form_single_domain_2D(
        const ParticleID& pid, const Particle& p, const FaceID& fid)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    // It considers the largest radius of 2D particles when drawing 3D shell,
    // so we don't need to consider the 2D-3D shell overlap here.

    const auto min_shell_radius = p.radius() * SINGLE_SHELL_FACTOR;

    bool min_shell_intruders_contains_multi = false;
    boost::container::small_vector<DomainID, 4> min_shell_intruders;

    for(const auto& sidp : shells_.list_shells_within_radius_2D(
            std::make_pair(p.position(), fid), min_shell_radius))
    {
        const auto did = *sidp.first.second.domain_id();
        unique_push_back(min_shell_intruders, did);
        if(domains_.at(did).second.is_multi())
        {
            min_shell_intruders_contains_multi = true;
        }
    }
    if(min_shell_intruders_contains_multi)
    {
        return min_shell_intruders;
    }

    // burst min shell intruders (if they are non-multi domains)

    Real max_distance = 0.0;
    boost::container::small_vector<DomainID, 4> intruders;
    for(const auto& did : min_shell_intruders)
    {
        for(const auto& result : burst_domain(did))
        {
            const auto& pid2 = result.first;
            const auto& p2   = result.second;
            const auto& fid2 = *(this->world_->on_which_face(pid2));

            const auto dist = ecell4::polygon::distance(this->world_->polygon(),
                std::make_pair(p.position(),  fid),
                std::make_pair(p2.position(), fid2)) - p2.radius();

            max_distance = std::max(max_distance, dist);

            const auto did2 = form_tight_domain_2D(pid2, p2, fid2);

            if(dist < p.radius() * SINGLE_SHELL_FACTOR)
            {
                unique_push_back(intruders, did2);
            }
        }
    }

    if( ! intruders.empty())
    {
        return intruders;
    }

    // if there are no min_shell_intruders, we need to calculate the max distance.
    if(min_shell_intruders.empty())
    {
        max_distance = max_circular_shell_size_at(p.position(), fid);
    }

    // No intruders here. draw single shell.

    const auto shell_size       = max_distance * SAFETY_SHRINK;
    const auto effective_radius = shell_size - p.radius();
    assert(0.0 < effective_radius);

    greens_functions::GreensFunction2DAbsSym gf(p.D(), effective_radius);

    const auto dt_escape   = gf.drawTime(this->world_->rng()->uniform(0.0, 1.0));
    const auto dt_reaction = this->draw_single_reaction_time(p.species());

    const auto dt = std::min(dt_escape,  dt_reaction);
    const auto event_kind = (dt_escape < dt_reaction) ?
            SingleCircularDomain::EventKind::Escape   :
            SingleCircularDomain::EventKind::Reaction ;

    const auto did = didgen_();
    const auto sid = sidgen_();

    ECELL4_NGFRD_LOG("2D domain did = ", did, ", dt_escape = ", dt_escape,
                     ", dt_reaction = ", dt_reaction, ", dt = ",  dt);

    CircularShell sh(p.radius(), Circle(shell_size, p.position(),
                     this->polygon().triangle_at(fid).normal()), fid);
    this->shells_.update_shell(sid, Shell(sh, did));

    SingleCircularDomain dom(event_kind, dt, world_->t(), sid, sh, pid,
            p.D(), effective_radius, std::move(gf));

    // add event with the same domain ID
    const auto evid = this->scheduler_.add(
            std::make_shared<event_type>(this->t() + dt, did));

    // update begin_time and re-insert domain into domains_ container
    this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));

    return boost::none;
}

bool NGFRDSimulator::form_pair_domain_2D(
        const ParticleID&, const Particle&, const FaceID&, const DomainID&)
{
    // TODO
    return false;
}

void NGFRDSimulator::form_domain_2D(
        const ParticleID& pid, const Particle& p, const FaceID& fid)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("form_domain_2D: forming domain for particle ", pid);

    const auto intruders = this->form_single_domain_2D(pid, p, fid);
    if( ! intruders)
    {
        return; // No intruders. Single succeeded.
    }

    // -----------------------------------------------------------------------
    // form_pair

    if(intruders->size() == 1)
    {
        if(form_pair_domain_2D(pid, p, fid, intruders->front()))
        {
            return; // pair is formed.
        }
    }

    // -----------------------------------------------------------------------
    // form_multi

    const auto pbc = this->world_->boundary();
    const Real multi_radius = p.radius() * SINGLE_SHELL_FACTOR;

    const auto host_id = didgen_();
    MultiDomain host(this->t());

    {
        const auto sid = sidgen_();

        // construct shell and assign it to shell container
        CircularShell sh(p.radius(), Circle(multi_radius, p.position(),
                        this->polygon().triangle_at(fid).normal()), fid);
        this->shells_.update_shell(sid, Shell(std::move(sh), host_id));

        host.add_particle(pid);
        host.add_shell(this->shells_.get_shell(sid));
    }

    // merge intruders that are already bursted or multis.
    for(const auto& did : *intruders)
    {
        if(domains_.at(did).second.is_multi())
        {
            // move all the particles and shells to the current host
            for(const auto& pid2 : domains_.at(did).second.as_multi().particle_ids())
            {
                host.add_particle(pid2);
            }
            for(const auto& shidp : domains_.at(did).second.as_multi().shells())
            {
                host.add_shell(shidp);
            }
            auto evid_dom = domains_.at(did);
            this->domains_.erase(did);
            this->scheduler_.remove(evid_dom.first);
        }
        else // 2D domain (non-multi 3D domain never overlap with 2D shell)
        {
            const auto result = this->burst_domain(did);
            assert(result.size() == 1);

            const auto  sid  = sidgen_();
            const auto& pid2 = result.front().first;
            const auto& p2   = result.front().second;
            const auto  fid2 = this->world_->on_which_face(pid2).value();

            const auto  new_shell_radius = p2.radius() * SINGLE_SHELL_FACTOR;

            // construct shell and assign it to shell container
            CircularShell sh(p.radius(), Circle(new_shell_radius, p.position(),
                            this->polygon().triangle_at(fid).normal()), fid);
            this->shells_.update_shell(sid, Shell(std::move(sh), host_id));

            host.add_particle(pid2);
            host.add_shell(this->shells_.get_shell(sid));
        }
    }

    this->recursively_merge_multis(host_id, host);

    host.determine_parameters(*(this->model()), *(this->world()));

    // add event with the same domain ID
    const auto evid = this->scheduler_.add(
            std::make_shared<event_type>(this->t() + host.dt(), host_id));
    ECELL4_NGFRD_LOG("event ", evid," added");

    this->domains_[host_id] = std::make_pair(evid, Domain(std::move(host)));
    return ;
}

boost::optional<std::pair<boost::container::small_vector<DomainID, 4>,
                          boost::container::small_vector<FaceID,   4>>>
NGFRDSimulator::form_single_domain_3D(const ParticleID& pid, const Particle& p)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    const Real largest_2D_particle = world_->largest_particle_radius_2D();
    const Real min_shell_radius = p.radius() * SINGLE_SHELL_FACTOR;

    const auto pbc = this->world_->boundary();

    boost::container::small_vector<DomainID, 4> intrusive_domains;
    boost::container::small_vector<FaceID, 4>   intrusive_faces;

    for(const auto& item : shells_.list_shells_within_radius_3D(
                p.position(), min_shell_radius))
    {
        const auto& shell = item.first.second;
        const auto  did   = shell.domain_id().get();
        if(std::find(intrusive_domains.begin(), intrusive_domains.end(), did) ==
                     intrusive_domains.end())
        {
            ECELL4_NGFRD_LOG("min shell intruder found: Shell ",
                             item.first.first, " in ", did);
            intrusive_domains.push_back(did);
        }
    }
    for(const auto& item : this->world_->list_faces_within_radius(
                p.position(), min_shell_radius + largest_2D_particle))
    {
        const auto fid = item.first.first;
        if(std::find(intrusive_faces.begin(), intrusive_faces.end(), fid) ==
                     intrusive_faces.end())
        {
            ECELL4_NGFRD_LOG("min shell intruder found: Face ", fid);
            intrusive_faces.push_back(fid);
        }
    }

    if( ! intrusive_domains.empty())
    {
        boost::container::small_vector<DomainID, 4> fatal_intruders;
        for(const auto& did : intrusive_domains)
        {
            if(domains_.at(did).second.is_multi())
            {
                fatal_intruders.push_back(did);
            }
            else // non-multi 3D domains
            {
                for(const auto& result : this->burst_domain(did))
                {
                    const auto& pid2 = result.first;
                    const auto& p2   = result.second;
                    const auto  dist = length(
                            pbc.periodic_transpose(p.position(), p2.position()) -
                            p2.position());

                    const auto did2 = this->form_tight_domain_3D(pid2, p2);

                    if(dist <= (p.radius() + p2.radius()) * SINGLE_SHELL_FACTOR)
                    {
                        fatal_intruders.push_back(did2);
                    }
                }
            }
        }
        intrusive_domains = fatal_intruders;
    }
    if( ! intrusive_faces.empty() || ! intrusive_domains.empty())
    {
        return std::make_pair(intrusive_domains, intrusive_faces);
    }

    // Here, it draws a modest sized shell that takes the other particles into
    // account. To maximizes the efficiency, it makes the shell radius
    // proportional to the particles' diffusion coefficient.

    const auto nearest_particle = world_->nearest_particle_3D(p.position(), pid);
    const auto nearest_face     = world_->nearest_face(p.position());

    ECELL4_NGFRD_LOG("nearest particles = ", nearest_particle);
    ECELL4_NGFRD_LOG("nearest face      = ", nearest_face);

    Real max_radius = std::min(pbc.edge_lengths()[0],
            std::min(pbc.edge_lengths()[1], pbc.edge_lengths()[2]));

    if( ! nearest_particle.empty())
    {
        const auto& dist     = nearest_particle.front().second;
        const auto& neighbor = nearest_particle.front().first.second;

        const auto modest_radius = p.radius() + (dist - p.radius()) * p.D() /
                                   (p.D() + neighbor.D());

        max_radius = std::min(max_radius, modest_radius);
    }
    if( ! nearest_face.empty())
    {
        max_radius = std::min(nearest_face.front().second - largest_2D_particle,
                              max_radius);
    }

    // check other 3D shells within max_radius.
    // Note: Here, we already subtract largest_2D_particle that is effective
    //       thickness of the 2D shells. So we just ignore 2D shells here, and
    //       only check 3D shells.
    for(const auto& item : this->shells_.list_shells_within_radius_3D(
                p.position(), max_radius))
    {
        max_radius = std::min(max_radius, item.second);
    }

    const auto shell_size = max_radius * SAFETY_SHRINK;
    assert(min_shell_radius <= shell_size);

    // paranoiac check
    assert(this->shells_.list_shells_within_radius_3D(p.position(), shell_size).empty());

    // -----------------------------------------------------------------------
    // XXX constructing spherical shell

    const auto effective_radius = shell_size - p.radius();
    assert(0.0 < effective_radius);
    greens_functions::GreensFunction3DAbsSym gf(p.D(), effective_radius);

    const auto dt_escape   = gf.drawTime(this->world_->rng()->uniform(0.0, 1.0));
    const auto dt_reaction = draw_single_reaction_time(p.species());

    const auto dt = std::min(dt_escape, dt_reaction);
    const auto event_kind  = (dt_escape < dt_reaction) ?
            SingleSphericalDomain::EventKind::Escape   :
            SingleSphericalDomain::EventKind::Reaction ;

    const auto did = didgen_();
    const auto sid = sidgen_();

    ECELL4_NGFRD_LOG("3D domain did = ", did, ", dt_escape = ", dt_escape,
                     ", dt_reaction = ", dt_reaction, ", dt = ",  dt);

    // construct shell and assign it to shell container
    SphericalShell sh(Sphere(p.position(), shell_size));
    this->shells_.update_shell(sid, Shell(sh, did));

    SingleSphericalDomain dom(event_kind, dt, world_->t(), sid, sh, pid,
            p.D(), effective_radius, std::move(gf));

    // add event with the same domain ID
    const auto evid = this->scheduler_.add(
            std::make_shared<event_type>(this->t() + dom.dt(), did));

    // update begin_time and re-insert domain into domains_ container
    this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));

    return boost::none;
}
bool NGFRDSimulator::form_pair_domain_3D(
        const ParticleID&, const Particle&, const DomainID&)
{
    // TODO
    return false;
}
void NGFRDSimulator::form_domain_3D(const ParticleID& pid, const Particle& p)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("form_domain_3D: forming domain for particle ", pid);

    // -----------------------------------------------------------------------
    // form_single

    const auto intruders = this->form_single_domain_3D(pid, p);
    if( ! intruders)
    {
        return; // No intruders. Single succeeded.
    }

    // -----------------------------------------------------------------------
    // form_pair

    if(intruders->first.size() == 1 && intruders->second.empty())
    {
        if(form_pair_domain_3D(pid, p, intruders->first.front()))
        {
            return; // pair is formed.
        }
    }

    // -----------------------------------------------------------------------
    // form_multi
    // 1. form a new, empty multi domain
    // 2. burst all the non-multis in the intruders
    // 3. merge the resulting particles into the multi domain

    const auto pbc = this->world_->boundary();
    const Real multi_radius = p.radius() * SINGLE_SHELL_FACTOR;

    const auto host_id = didgen_();
    MultiDomain host(this->t());

    // add current shell
    {
        const auto sid = sidgen_();

        // construct shell and assign it to shell container
        SphericalShell sh(Sphere(p.position(), multi_radius));
        this->shells_.update_shell(sid, Shell(std::move(sh), host_id));

        host.add_particle(pid);
        host.add_shell(this->shells_.get_shell(sid));
    }

    // add 3D intrudres collected by form_single_domain_3D.
    // It contains only multi or tight domains.
    for(const auto& did : intruders->first)
    {
        if(domains_.at(did).second.is_multi())
        {
            // move all the particles and shells to the current host
            for(const auto& pid2 : domains_.at(did).second.as_multi().particle_ids())
            {
                host.add_particle(pid2);
            }
            for(const auto& shidp : domains_.at(did).second.as_multi().shells())
            {
                host.add_shell(shidp);
            }
            auto evid_dom = domains_.at(did);
            this->domains_.erase(did);
            this->scheduler_.remove(evid_dom.first);
        }
        else // single tight domain. burst and add it to multi.
        {
            const auto result = this->burst_domain(did);
            assert(result.size() == 1);

            const auto  sid  = sidgen_();
            const auto& pid2 = result.front().first;
            const auto& p2   = result.front().second;

            const auto  new_shell_radius = p2.radius() * SINGLE_SHELL_FACTOR;

            // construct shell and assign it to shell container
            SphericalShell sh(Sphere(p2.position(), new_shell_radius));
            this->shells_.update_shell(sid, Shell(std::move(sh), host_id));

            host.add_particle(pid2);
            host.add_shell(this->shells_.get_shell(sid));
        }
    }

    // check overlap between the formed multi and other shells.
    // 2D shells are also collected inside this.
    this->recursively_merge_multis(host_id, host);

    host.determine_parameters(*(this->model()), *(this->world()));

    // add event with the same domain ID
    const auto evid = this->scheduler_.add(
            std::make_shared<event_type>(this->t() + host.dt(), host_id));
    ECELL4_NGFRD_LOG("event ", evid," added");

    this->domains_[host_id] = std::make_pair(evid, Domain(std::move(host)));
    return;
}

void NGFRDSimulator::recursively_merge_multis(const DomainID& host_id, MultiDomain& host)
{
    const Real largest_2D_particle = world_->largest_particle_radius_2D();
    std::vector<DomainID> intrusive_domains;

    // collect intrusive domains

    for(const auto& sidp : host.shells())
    {
        const auto& shid = sidp.first;
        const auto& sh   = sidp.second;

        if(sh.is_spherical())
        {
            const auto& shape = sh.as_spherical().shape();

            // collect 3D shell
            for(const auto& item : shells_.list_shells_within_radius_3D(
                        shape.position(), shape.radius()))
            {
                const auto& shell = item.first.second;
                const auto  did   = shell.domain_id().get();
                if(did == host_id) {continue;}
                unique_push_back(intrusive_domains, did);
            }

            // collect 2D shell
            for(const auto& item : this->world_->list_faces_within_radius(
                        shape.position(), shape.radius() + largest_2D_particle))
            {
                const auto fid = item.first.first;
                this->collect_possibly_overlapping_2D_domains(
                    shape.position(), shape.radius(), fid, host_id, intrusive_domains);
            }
        }
        else if(sh.is_circular())
        {
            const auto& fid   = sh.as_circular().fid();
            const auto& shape = sh.as_circular().shape();

            // collect 3D shell
            for(const auto& item : shells_.list_shells_within_radius_3D(
                        shape.position(), shape.radius()))
            {
                const auto& shell = item.first.second;
                const auto  did   = shell.domain_id().get();
                if(did == host_id) {continue;}
                // non-multi 3D domains never overlap with a 2D domain because
                // they do not intersect with a face (with a tickness of
                // largest_particle_radius_2D).
                if(domains_.at(did).second.is_multi())
                {
                    unique_push_back(intrusive_domains, did);
                }
            }

            // collect 2D shell
            for(const auto& item : shells_.list_shells_within_radius_2D(
                    std::make_pair(shape.position(), fid), shape.radius()))
            {
                const auto& shell = item.first.second;
                const auto  did   = shell.domain_id().get();
                if(did == host_id) {continue;}
                unique_push_back(intrusive_domains, did);
            }
        }
        else
        {
            throw_exception<IllegalState>(
                    "Multi has a shell that is not a spherical nor circular!");
        }
    }

    // burst all the domains and merge them into multi

    const auto n_shells = host.shells().size();
    for(const auto& did : intrusive_domains)
    {
        for(const auto& result : burst_domain(did))
        {
            const auto& pid = result.first;
            const auto& p   = result.second;
            if(const auto fid = world_->on_which_face(pid))
            {
                merge_into_multi_or_form_tight_domain_2D(host_id, host, pid, p, *fid);
            }
            else // 3D particle
            {
                merge_into_multi_or_form_tight_domain_3D(host_id, host, pid, p);
            }
        }
    }
    if(host.shells().size() != n_shells) // multi expanded
    {
        return recursively_merge_multis(host_id, host);
    }
    return ;
}

void NGFRDSimulator::merge_into_multi_or_form_tight_domain_3D(
        const DomainID& host_id, MultiDomain& host,
        const ParticleID& pid, const Particle& p)
{
    const auto min_radius = p.radius() * SINGLE_SHELL_FACTOR;

    const auto pbc = this->world_->boundary();
    for(const auto& sidp : host.shells())
    {
        const auto sph  = sidp.second.bounding_sphere();
        const auto dist = length(
                pbc.periodic_transpose(p.position(), sph.position()) -
                sph.position());

        if(dist <= min_radius + sph.radius())
        {
            const auto sid = sidgen_();
            SphericalShell sh(Sphere(p.position(), min_radius));
            this->shells_.update_shell(sid, Shell(sh, host_id));

            host.add_particle(pid);
            host.add_shell(this->shells_.get_shell(sid));
        }
        else
        {
            form_tight_domain_3D(pid, p);
        }
    }
    return ;
}
void NGFRDSimulator::merge_into_multi_or_form_tight_domain_2D(
        const DomainID& host_id, MultiDomain& host,
        const ParticleID& pid, const Particle& p, const FaceID& fid)
{
    const auto min_radius = p.radius() * SINGLE_SHELL_FACTOR;

    const auto pbc = this->world_->boundary();
    for(const auto& sidp : host.shells())
    {
        Real dist = 0;

        const auto& sh = sidp.second;
        if(sh.is_spherical())
        {
            const auto& sph = sh.as_spherical().shape();
            dist = length(pbc.periodic_transpose(p.position(), sph.position()) -
                sph.position()) - sph.radius();
        }
        else if(sh.is_circular())
        {
            dist = ecell4::polygon::distance(this->world_->polygon(),
                std::make_pair(p.position(), fid),
                std::make_pair(sh.as_circular().position(), sh.as_circular().fid())
                ) - sh.as_circular().shape().radius();
        }
        else
        {
            throw_exception<IllegalState>(
                    "Multi has a shell that is not a spherical nor circular!");
        }

        if(dist <= min_radius)
        {
            const auto sid = sidgen_();
            CircularShell sh(p.radius(), Circle(min_radius, p.position(),
                            this->polygon().triangle_at(fid).normal()), fid);
            this->shells_.update_shell(sid, Shell(sh, host_id));

            host.add_particle(pid);
            host.add_shell(this->shells_.get_shell(sid));
        }
        else
        {
            form_tight_domain_2D(pid, p, fid);
        }
    }
    return ;
}

DomainID NGFRDSimulator::form_tight_domain_2D(
        const ParticleID& pid, const Particle& p, const FaceID& fid)
{
    const auto did = didgen_();
    const auto sid = sidgen_();

    CircularShell sh(p.radius(), Circle(p.radius(), p.position(),
                    this->polygon().triangle_at(fid).normal()), fid);

    this->shells_.update_shell(sid, Shell(sh, did));

    SingleCircularDomain dom(SingleCircularDomain::EventKind::Escape,
        /*dt = */ 0.0, world_->t(), sid, sh, pid, p.D(), /*radius = */ 0.0,
        greens_functions::GreensFunction2DAbsSym(p.D(), 0.0));

    // add event with the same domain ID
    const auto evid = this->scheduler_.add(
            std::make_shared<event_type>(this->t(), did));

    // update begin_time and re-insert domain into domains_ container
    this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));
    return did;
}

DomainID NGFRDSimulator::form_tight_domain_3D(const ParticleID& pid, const Particle& p)
{
    const auto did = didgen_();
    const auto sid = sidgen_();

    // construct shell and assign it to shell container
    SphericalShell sh(Sphere(p.position(), p.radius()));

    this->shells_.update_shell(sid, Shell(sh, did));

    SingleSphericalDomain dom(SingleSphericalDomain::EventKind::Escape,
        /*dt = */ 0.0, world_->t(), sid, sh, pid, p.D(), /*radius = */ 0.0,
        greens_functions::GreensFunction3DAbsSym(p.D(), 0.0));

    // add event with the same domain ID
    const auto evid = this->scheduler_.add(
            std::make_shared<event_type>(this->t(), did));

    // update begin_time and re-insert domain into domains_ container
    this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));
    return did;
}

boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
NGFRDSimulator::fire_single_circular(const DomainID& did, SingleCircularDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("firing single circular: ", did);
    ECELL4_NGFRD_LOG("included shell: ", dom.shell_id());

    if(dom.dt() == 0.0) // means it is a tight domain. do nothing.
    {
        this->shells_.remove_shell(dom.shell_id());
        return boost::container::small_vector<std::pair<ParticleID, Particle>, 4>{
                world_->get_particle(dom.particle_id())
            };
    }

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions;
    SingleCircularPropagator prop(did,
            *(this->model_), *(this->world_), *this, *(this->world_->rng()),
            SINGLE_CONICAL_MAX_RETRY, last_reactions);

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4> results;
    for(const auto& pid : prop(dom))
    {
        results.push_back(world_->get_particle(pid));
    }

    // remove shell from shell container
    this->shells_.remove_shell(dom.shell_id());

    if( ! last_reactions.empty())
    {
        std::copy(last_reactions.begin(), last_reactions.end(),
                  std::back_inserter(last_reactions_));
    }
    return results;
}

boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
NGFRDSimulator::fire_single_spherical(const DomainID& did, SingleSphericalDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("firing single spherical: ", did);
    ECELL4_NGFRD_LOG("included shell: ", dom.shell_id());

    if(dom.dt() == 0.0) // means it is a tight domain
    {
        this->shells_.remove_shell(dom.shell_id());
        return boost::container::small_vector<std::pair<ParticleID, Particle>, 4>{
                world_->get_particle(dom.particle_id())
            };
    }

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions;
    SingleSphericalPropagator prop(
            did, *(this->model_), *(this->world_), *this,
            *(this->world_->rng()), SINGLE_SPHERICAL_MAX_RETRY, last_reactions);

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4> results;
    for(const auto& pid : prop(dom))
    {
        results.push_back(world_->get_particle(pid));
    }

    // remove shell from shell container
    this->shells_.remove_shell(dom.shell_id());

    if( ! last_reactions.empty())
    {
        std::copy(last_reactions.begin(), last_reactions.end(),
                  std::back_inserter(last_reactions_));
    }

    return results;
}

boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
NGFRDSimulator::fire_multi(const DomainID& did, MultiDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("firing multi: ", did);
    ECELL4_NGFRD_LOG("included shells: ", dom.shells());

    dom.step(did, *(this->model_), *this, *(this->world_));

    // XXX: If no (reaction, escapement) happens, we don't need to break it
    //      down into independent domains. For the efficiency, it re-inserts
    //      this domain into scheduler and returns nothing.
    //          Since nothing is returned, no domain will be formed at the end
    //      of this step.
    if(dom.eventkind() == MultiDomain::EventKind::None)
    {
        // add event with the same domain ID
        const auto evid = scheduler_.add(
                // it performed stepping. the next event is at t+dt
                std::make_shared<event_type>(this->t() + dom.dt(), did));

        // update begin_time and re-insert domain into domains_ container
        dom.begin_time() = this->t() + dom.dt();
        domains_[did] = std::make_pair(evid, Domain(std::move(dom)));

        return {/* All particles are re-inserted as Multi! */};
    }
    // something happens. remove multi and return resulting particles
    // to re-wrap the particles by new domain

    const auto& last_reactions = dom.last_reactions();
    if( ! last_reactions.empty())
    {
        std::copy(last_reactions.begin(), last_reactions.end(),
                  std::back_inserter(last_reactions_));
    }

    // remove shells
    for(const auto& sidp : dom.shells())
    {
        ECELL4_NGFRD_LOG("removing shell: ", sidp);
        this->shells_.remove_shell(sidp.first);
    }

    // collect resulting particles
    boost::container::small_vector<std::pair<ParticleID, Particle>, 4> retval;
    for(const auto& pid : dom.particle_ids())
    {
        retval.push_back(world_->get_particle(pid));
    }
    return retval;
}

boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
NGFRDSimulator::burst_single_circular(const DomainID& did, SingleCircularDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("bursting single circular: ", did);
    ECELL4_NGFRD_LOG("included shell: ", dom.shell_id());

    if(dom.dt() == 0.0) // means it is a tight domain. do nothing.
    {
        this->shells_.remove_shell(dom.shell_id());
        return boost::container::small_vector<std::pair<ParticleID, Particle>, 4>{
                world_->get_particle(dom.particle_id())
            };
    }

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions;
    SingleCircularPropagator prop(did,
            *(this->model_), *(this->world_), *this, *(this->world_->rng()),
            SINGLE_CONICAL_MAX_RETRY, last_reactions);

    // While bursting single domain, it never dissociates (becuase the time when
    // it dissociates is already calculated and considered in the domain dt).
    boost::container::small_vector<std::pair<ParticleID, Particle>, 4> results;
    results.push_back(world_->get_particle(prop.burst(dom, this->t())));

    // remove shell from shell container
    this->shells_.remove_shell(dom.shell_id());

    assert(last_reactions.empty());
    return results;
}

boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
NGFRDSimulator::burst_single_spherical(const DomainID& did, SingleSphericalDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("bursting single spherical: ", did);
    ECELL4_NGFRD_LOG("included shell: ", dom.shell_id());

    if(dom.dt() == 0.0) // it means this is a tight domain
    {
        this->shells_.remove_shell(dom.shell_id());
        return boost::container::small_vector<std::pair<ParticleID, Particle>, 4>{
                world_->get_particle(dom.particle_id())
            };
    }

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions;
    SingleSphericalPropagator prop(
            did, *(this->model_), *(this->world_), *this,
            *(this->world_->rng()), SINGLE_SPHERICAL_MAX_RETRY, last_reactions);

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4> results;
    results.push_back(world_->get_particle(prop.burst(dom, this->t())));

    // remove shell from shell container
    this->shells_.remove_shell(dom.shell_id());

    assert(last_reactions.empty());
    return results;
}

boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
NGFRDSimulator::burst_multi(const DomainID& did, MultiDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("bursting multi: ", did);
    ECELL4_NGFRD_LOG("included shells: ", dom.shells());

    // step until this->t()
    const auto dt = this->t() - dom.begin_time();
    dom.step(did, *(this->model_), *this, *(this->world_), dt);

    assert(shells_.diagnosis()); // XXX

    // The multi is going to be bursted, so we don't need to check if
    // any reaction or escape happened. We just remove it.

    // remove shells
    for(const auto& sidp : dom.shells())
    {
        ECELL4_NGFRD_LOG("removing shell: ", sidp);
        this->shells_.remove_shell(sidp.first);
    }

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4> retval;
    for(const auto& pid : dom.particle_ids())
    {
        retval.push_back(world_->get_particle(pid));
    }
    return retval;
}

Real NGFRDSimulator::draw_single_reaction_time(const Species& sp) const
{
    const auto& rules = this->model()->query_reaction_rules(sp);
    if(rules.empty())
    {
        return std::numeric_limits<Real>::infinity();
    }

    const auto k_tot = std::accumulate(rules.begin(), rules.end(), Real(0.0),
            [](const Real acc, const ReactionRule& rule) noexcept -> Real {
                return acc + rule.k();
            });

    if(k_tot == std::numeric_limits<Real>::infinity())
    {
        return 0.0;
    }
    else
    {
        const Real rnd(this->rng()->uniform(0., 1.));
        if(rnd <= 0.0)
        {
            return std::numeric_limits<Real>::infinity();
        }
        else
        {
            return -std::log(rnd) / k_tot;
        }
    }
}

} // ngfrd
} // ecell4
