#include <ecell4/ngfrd/NGFRDSimulator.hpp>
#include <ecell4/ngfrd/SingleSphericalPropagator.hpp>

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

void NGFRDSimulator::form_domain_2D(
        const ParticleID& pid, const Particle& p, const FaceID& fid)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    // TODO: Currently we always form a multi domain.
    ECELL4_NGFRD_LOG("form_domain_2D: forming domain for particle ", pid);

    // -----------------------------------------------------------------------
    // form_multi
    const Real multi_radius = p.radius() * SINGLE_SHELL_FACTOR;

    // list 2D domains within the multi shell
    std::vector<DomainID> intruders;
    for(const auto& item : shells_.list_shells_within_radius_2D(
                std::make_pair(p.position(), fid), multi_radius))
    {
        const auto& shell = item.first.second;
        const auto  did   = shell.domain_id().get();
        if(std::find(intruders.begin(), intruders.end(), did) == intruders.end())
        {
            ECELL4_NGFRD_LOG("intruder found: 2DShell ", item.first.first, " in ", did);
            intruders.push_back(did);
        }
    }
    ECELL4_NGFRD_LOG("form_domain_2D: ", intruders.size(), " intrusive shells found");

    // list 3D multi domains that overlap with the shell.
    // Here, we list 3D shells that are within the bounding sphere and overlap
    // with the polygon.
    for(const auto& item : shells_.list_shells_within_radius_3D(
                p.position(), multi_radius))
    {
        const auto& shell = item.first.second;
        const auto  did   = shell.domain_id().get();

        if(!this->domains_.at(did).second.is_multi())
        {
            // non-Multi domain never overlaps with the polygon (?)
            continue;
        }
        // Multi 3D Shell is always spherical.
        const auto& sh = shell.as_spherical();
        if(this->polygon().has_overlapping_faces(sh.shape().position(), sh.shape().radius()))
        {
            // This shell is within the bounding sphere of 2D shell and overlaps
            // with Polygon. Insert it to multi (that does not always mean that
            // the domain overlaps with 2D shell, but a nice approximation (sort of))
            if(std::find(intruders.begin(), intruders.end(), did) == intruders.end())
            {
                ECELL4_NGFRD_LOG("intruder found: 3DShell ", item.first.first, " in ", did);
                intruders.push_back(did);
            }
        }
    }

    if(intruders.empty())
    {
        // form a new multi!
        const auto did = didgen_();
        const auto sid = sidgen_();

        Shell sh(CircularShell(p.radius(), Circle(multi_radius, p.position(),
                        this->polygon().triangle_at(fid).normal()), fid), did);

        ECELL4_NGFRD_LOG("shell created: ", sid);
        this->shells_.update_shell(sid, sh);
        ECELL4_NGFRD_LOG("shell container updated");

        MultiDomain dom(this->t());
        dom.add_particle(pid);
        dom.add_shell(shells_.get_shell(sid));
        dom.determine_parameters(*(this->model()), *(this->world()));

        ECELL4_NGFRD_LOG("multi domain created ", did);

        // add event with the same domain ID
        const auto evid = this->scheduler_.add(
                std::make_shared<event_type>(this->t() + dom.dt(), did));

        ECELL4_NGFRD_LOG("event added");

        // update begin_time and re-insert domain into domains_ container
        this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));
        return;
    }

    ECELL4_NGFRD_LOG("intruder found. merge all domains");
    // XXX Currently all the domains are multi domains. merge all those multi domains.
    // Later we need to burst domains and check if the resulting particle should
    // be in Multi or form Single

    const auto host_id = intruders.back();
    intruders.pop_back();

    const auto sid = sidgen_();
    Shell sh(CircularShell(p.radius(), Circle(multi_radius, p.position(),
                    this->polygon().triangle_at(fid).normal()), fid), host_id);
    ECELL4_NGFRD_LOG("shell created ", sid);
    this->shells_.update_shell(sid, sh);
    ECELL4_NGFRD_LOG("shell inserted");

    assert(domains_.at(host_id).second.is_multi());
    auto& host = domains_.at(host_id).second.as_multi();
    host.add_particle(pid);
    host.add_shell(shells_.get_shell(sid));
    for(const auto& did : intruders)
    {
        const auto dom_iter = domains_.find(did);

        const auto evid = dom_iter->second.first;
        scheduler_.remove(evid);

        assert(dom_iter->second.second.is_multi());
        const auto& dom = dom_iter->second.second.as_multi();
        for(const auto& pid : dom.particle_ids())
        {
            host.add_particle(pid);
        }
        for(const auto& sidp : dom.shells())
        {
            this->shells_.at(sidp.first).second.domain_id() = host_id;
            host.add_shell(sidp);
        }
        domains_.erase(dom_iter);
    }
    host.determine_parameters(*(this->model()), *(this->world()));

    return ;
}

boost::optional<std::pair<boost::container::small_vector<DomainID, 4>,
                          boost::container::small_vector<FaceID,   4>>>
NGFRDSimulator::form_single_domain_3D(const ParticleID& pid, const Particle& p)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    const Real largest_2D_particle = world_->largest_particle_radius_2D();
    const Real min_shell_radius = p.radius() * SINGLE_SHELL_FACTOR;

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
    if( ! intrusive_faces.empty())
    {
        return std::make_pair(intrusive_domains, intrusive_faces);
    }
    if( ! intrusive_domains.empty())
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

    const auto pbc = this->world_->boundary();
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

    // construct shell and assign it to shell container
    SphericalShell sh(Sphere(p.position(), shell_size));
    this->shells_.update_shell(sid, Shell(sh, did));

    SingleSphericalDomain dom(event_kind, dt, world_->t(), sid, sh, pid,
            p.D(), shell_size - p.radius(), std::move(gf));

    // add event with the same domain ID
    const auto evid = this->scheduler_.add(
            std::make_shared<event_type>(this->t() + dom.dt(), did));

    // update begin_time and re-insert domain into domains_ container
    this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));

    return boost::none;
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
    // form_pair (if size of intrusive particles == 2 && no intrusive faces)
    //
    // TODO

    // -----------------------------------------------------------------------
    // form_multi
    const Real multi_radius = p.radius() * SINGLE_SHELL_FACTOR;

    // list 3D domains within the multi shell
    std::vector<DomainID> intrusive_domains;
    std::copy(intruders->first.begin(), intruders->first.end(),
              std::back_inserter(intrusive_domains));

    // list faces
    const auto pbc = this->world_->boundary();
    for(const auto& fid : intruders->second)
    {
        // shells on this triangle
        if(const auto shids = shells_.shells_on(fid))
        {
            for(const auto shid : *shids)
            {
                // check overlap with the bounding sphere of 2D shell
                // as an approximation
                const auto sh = this->shells_.get_shell(shid).second;
                const auto bs = sh.bounding_sphere();
                const auto dist = length(
                        pbc.periodic_transpose(p.position(), bs.position()) -
                        bs.position());
                if(dist <= multi_radius + bs.radius())
                {
                    intrusive_domains.push_back(*sh.domain_id());
                }
            }
        }
        // shells on the vertices of the triangle
        for(const auto vid : this->polygon().vertices_of(fid))
        {
            if(const auto shids = shells_.shells_on(vid))
            {
                for(const auto shid : *shids)
                {
                    // check overlap with the bounding sphere of 2D shell
                    // as an approximation
                    const auto sh = this->shells_.get_shell(shid).second;
                    const auto bs = sh.bounding_sphere();
                    const auto dist = length(
                            pbc.periodic_transpose(p.position(), bs.position()) -
                            bs.position());
                    if(dist <= multi_radius + bs.radius())
                    {
                        intrusive_domains.push_back(*sh.domain_id());
                    }
                }
            }
        }
        // shells on the neighboring triangles
        for(const auto nfid : this->polygon().neighbor_faces_of(fid))
        {
            if(const auto shids = shells_.shells_on(nfid))
            {
                for(const auto shid : *shids)
                {
                    // check overlap with the bounding sphere of 2D shell
                    // as an approximation
                    const auto sh = this->shells_.get_shell(shid).second;
                    const auto bs = sh.bounding_sphere();
                    const auto dist = length(
                            pbc.periodic_transpose(p.position(), bs.position()) -
                            bs.position());
                    if(dist <= multi_radius + bs.radius())
                    {
                        intrusive_domains.push_back(*sh.domain_id());
                    }
                }
            }
        }
    }

    assert( ! intrusive_domains.empty());
    ECELL4_NGFRD_LOG("intruder found. merge all domains");

    // 1. form a new, empty multi domain
    // 2. burst all the non-multis in the intruders
    // 3. merge all the resulting particles into the multi domain

    const auto host_id = didgen_();
    MultiDomain host(this->t());
    for(const auto did : intrusive_domains)
    {
        const auto& dom = this->domains_.at(did).second;

        // non-Multi domain will never leak particle from its shell.
        //
        // In case of (single-)reaction, the resulting particles might stick
        // out from the shell. But in the case of bursting, the particle
        // cannot react because the bursting always occurs before the
        // predefined event time. If reaction happens before the bursting,
        // the event kind should be SingleReaction.
        //
        // Thus we can safely burst the domains.

        if(dom.is_multi())
        {
            for(const auto& pid : dom.as_multi().particle_ids())
            {
                host.add_particle(pid);
            }
            for(const auto& shidp : dom.as_multi().shells())
            {
                host.add_shell(shidp);
            }

            using std::swap;
            // remove the domain (Note that after this `dom` will point an empty domain)
            std::pair<event_id_type, Domain> evid_dom;
            swap(evid_dom, domains_.at(did));

            // remove domain and event from container
            this->domains_.erase(did);
            this->scheduler_.remove(evid_dom.first);

            // `dom` is merged into `host`.
        }
        else if(dom.is_2D())
        {
            // burst the 2D domain. Here it will never react (because it is
            // bursted before any event).
            // It bursts the domain and if the resulting particle is too close
            // to any of the shells in the host domain, it merge the particle
            // into the domain. Otherwise, it forms tight shell.

            for(const auto& pidp : this->burst_domain(did))
            {
                const auto& pid = pidp.first;
                const auto& p   = pidp.second;
                const auto& fid = this->world_->on_which_face(pid).value();

                const auto multi_radius = p.radius() * SINGLE_SHELL_FACTOR;

                // Here, it calculates 3D distance. 2D distance is always longer,
                // so there could be a "needless" multi. But I just ignore this
                // inefficiency for the sake of simplicity, currently.
                //
                // TODO(?): add ShellDistanceExactCalculator that takes optional
                //          FaceID for 2D position
                const ShellDistanceCalculator shell_distance(p.position(), pbc);

                // This min_dist is not the distance between the surfaces, but
                // distance between "shell surface" and "a point (center of the
                // particle)".
                Real min_dist = std::numeric_limits<Real>::infinity();
                for(const auto& shidp : host.shells())
                {
                    const auto dist = visit(shell_distance, shidp.second);
                    min_dist = std::min(min_dist, dist);
                }
                if(min_dist <= multi_radius)
                {
                    // form a (2D) multi shell
                    const auto sid = sidgen_();

                    Shell sh(CircularShell(p.radius(), Circle(multi_radius, p.position(),
                        this->polygon().triangle_at(fid).normal()), fid), host_id);

                    this->shells_.update_shell(sid, sh);

                    host.add_particle(pid);
                    host.add_shell(std::make_pair(sid, sh));
                }
                else
                {
                    this->form_tight_domain_2D(pid, p, fid);
                }
            }
        }
        else // dom.is_3D
        {
            // burst the 3D domain. Here it will never react (because it is
            // bursted before any event).
            // It bursts the domain and if the resulting particle is too close
            // to any of the shells in the host domain, it merge the particle
            // into the domain. Otherwise, it forms tight shell.
            for(const auto& pidp : this->burst_domain(did))
            {
                const auto& pid = pidp.first;
                const auto& p   = pidp.second;

                const auto multi_radius = p.radius() * SINGLE_SHELL_FACTOR;

                const ShellDistanceCalculator shell_distance(p.position(), pbc);

                // This min_dist is not the distance between the surfaces, but
                // distance between "shell surface" and "a point (center of the
                // particle)".
                Real min_dist = std::numeric_limits<Real>::infinity();
                for(const auto& shidp : host.shells())
                {
                    const auto dist = visit(shell_distance, shidp.second);
                    min_dist = std::min(min_dist, dist);
                }
                if(min_dist <= multi_radius)
                {
                    // form a (3D) multi shell
                    const auto sid = sidgen_();

                    Shell sh(SphericalShell(Sphere(p.position(), multi_radius)),
                             host_id);

                    this->shells_.update_shell(sid, sh);

                    host.add_particle(pid);
                    host.add_shell(std::make_pair(sid, sh));
                }
                else
                {
                    this->form_tight_domain_3D(pid, p);
                }
            }
        }
    }

    host.determine_parameters(*(this->model()), *(this->world()));

    // add event with the same domain ID
    const auto evid = this->scheduler_.add(
            std::make_shared<event_type>(this->t() + host.dt(), host_id));
    ECELL4_NGFRD_LOG("event ", evid," added");

    this->domains_[host_id] = std::make_pair(evid, Domain(std::move(host)));
    return;
}

void NGFRDSimulator::form_tight_domain_2D(
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
    return;
}

void NGFRDSimulator::form_tight_domain_3D(const ParticleID& pid, const Particle& p)
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
    return;
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

    throw_exception<NotImplemented>("TODO: ", ECELL4_NGFRD_LOG_FUNCTION_NAME);
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

    throw_exception<NotImplemented>("TODO: ", ECELL4_NGFRD_LOG_FUNCTION_NAME);
}

boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
NGFRDSimulator::burst_single_spherical(const DomainID& did, SingleSphericalDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("bursting single spherical: ", did);
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
    // any reaction or escape happened.

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
