#include <ecell4/ngfrd/NGFRDSimulator.hpp>
#include <ecell4/ngfrd/SingleSphericalPropagator.hpp>
#include <ecell4/ngfrd/SingleCircularPropagator.hpp>
#include <ecell4/ngfrd/SingleConicalPropagator.hpp>
#include <ecell4/ngfrd/PairSphericalPropagator.hpp>
#include <ecell4/ngfrd/PairCircularPropagator.hpp>

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
constexpr std::size_t NGFRDSimulator::SINGLE_CIRCULAR_MAX_RETRY;
constexpr std::size_t NGFRDSimulator::PAIR_SPHERICAL_MAX_RETRY;
constexpr std::size_t NGFRDSimulator::PAIR_CIRCULAR_MAX_RETRY;

boost::optional<boost::container::small_vector<DomainID, 4>>
NGFRDSimulator::form_single_conical_domain_2D(
        const ParticleID& pid, const Particle& p, const FaceID& fid)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("pid = ", pid, ", p = ", p, ", fid = ", fid);

    ECELL4_NGFRD_LOG("there are no intrusive domains, but an intrusive vertex exists. "
                     "trying conical shell and then try pair or multi.");

    // there are no intruders, but max_distance (came from the geometric
    // constraint) is too small. It means that there is a vertex too close
    // to the particle. Try forming a conical shell.

    const auto& poly = world_->polygon();

    // 1. detect the nearest vertex.
    Real dist_to_vtx = std::numeric_limits<Real>::infinity();
    VertexID vtxid;
    for(const auto& vid : poly.vertices_of(fid))
    {
        const auto dist2 = ecell4::polygon::distance_sq(poly,
                std::make_pair(p.position(), fid),
                std::make_pair(poly.position_at(vid), vid));
        if(dist2 < dist_to_vtx)
        {
            dist_to_vtx = dist2;
            vtxid = vid;
        }
    }
    dist_to_vtx = std::sqrt(dist_to_vtx);

    // TODO: consider changing this factor by the shell kind
    const Real min_conical_shell_size = dist_to_vtx * SINGLE_SHELL_FACTOR + p.radius();

    // 2. determine the shell size
    //   - find the shortest edge extending from the vertex.
    //   - check the intruders

    Real max_shell_size = std::numeric_limits<Real>::infinity();
    for(const auto& eid : poly.outgoing_edges(vtxid))
    {
        // To avoid conical-conical collision, make max_shell_size smaller
        // than the half length of edges.
        max_shell_size = std::min(max_shell_size, poly.length_of(eid) * 0.5);

        // check shells on faces which connect to the vertex
        const auto neighbor_fid = poly.face_of(eid);
        if(const auto possible_intruders = world_->particles_on(neighbor_fid))
        {
            for(const auto& possible_intruder_id : *possible_intruders)
            {
                const auto possible_intruder =
                    world_->get_particle(possible_intruder_id).second;

                max_shell_size = std::min(max_shell_size,
                    ecell4::polygon::distance_sq(poly,
                        std::make_pair(possible_intruder.position(), neighbor_fid),
                        std::make_pair(poly.position_at(vtxid), vtxid)));
            }
        }
    }
    if(max_shell_size < min_conical_shell_size ||
       max_shell_size <= dist_to_vtx + p.radius())
    {
        // maximum possible shell size is smaller than the minimum conical.
        // Conical shell cannot be formed.
        return {/*empty*/};
    }

    // 3. form conical shell.
    const auto shell_size       = max_shell_size * SAFETY_SHRINK;
    const auto effective_radius = shell_size - p.radius();
    const auto initial_position = dist_to_vtx;
    const auto phi              = poly.apex_angle_at(vtxid);
    assert(0.0 < effective_radius);

    greens_functions::GreensFunction2DRefWedgeAbs
        gf(p.D(), initial_position, effective_radius, phi);

    const auto dt_escape   = gf.drawTime(this->world_->rng()->uniform(0.0, 1.0));
    const auto dt_reaction = this->draw_single_reaction_time(p.species());

    const auto dt = std::min(dt_escape,  dt_reaction);
    const auto event_kind = (dt_escape < dt_reaction) ?
            SingleConicalDomain::EventKind::Escape   :
            SingleConicalDomain::EventKind::Reaction ;

    const auto did = didgen_();
    const auto sid = sidgen_();

    ECELL4_NGFRD_LOG("2D Conical domain: did = ", did,
            ", dt_escape = ", dt_escape, ", dt_reaction = ", dt_reaction,
            ", dt = ",  dt);

    ConicalShell sh(p.radius(), ConicalSurface(p.position(), phi, shell_size), vtxid);
    this->shells_.update_shell(sid, Shell(sh, did));

    SingleConicalDomain dom(event_kind, dt, world_->t(), sid, sh,
            pid, p.D(), effective_radius, std::move(gf));

    // add event with the same domain ID
    const auto evid = this->scheduler_.add(
            std::make_shared<event_type>(this->t() + dt, did));

    // update begin_time and re-insert domain into domains_ container
    this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));

    ECELL4_NGFRD_LOG("2D single circular domain formed!");

    return boost::none;
}


boost::optional<boost::container::small_vector<DomainID, 4>>
NGFRDSimulator::form_single_domain_2D(
        const ParticleID& pid, const Particle& p, const FaceID& fid)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("pid = ", pid, ", p = ", p, ", fid = ", fid);

    assert(world_->on_which_face(pid).value() == fid);

    // It considers the largest radius of 2D particles when drawing 3D shell,
    // so we don't need to consider the 2D-3D shell overlap here.

    const auto min_shell_radius = p.radius() * SINGLE_SHELL_FACTOR;
    ECELL4_NGFRD_LOG("min_shell_radius = ", min_shell_radius);

    // first, list all the min-shell intruders.

    bool min_shell_intruders_contains_multi = false;
    boost::container::small_vector<DomainID, 4> min_shell_intruders;

    for(const auto& sidp : shells_.list_shells_within_radius_2D(
            std::make_pair(p.position(), fid), min_shell_radius))
    {
        const auto did = *sidp.first.second.domain_id();

        ECELL4_NGFRD_LOG("domain ", did, " is at the ", sidp.second, " distant");

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

    // Then, burst min shell intruders (if they are non-multi domains).
    // The bursted particles will have a tight (zero-radius) domain.

    this->world_->poly_con_.diagnosis();

    boost::container::small_vector<DomainID, 4> intruders;
    for(const auto& did : min_shell_intruders)
    {
        ECELL4_NGFRD_LOG("bursting min shell intruder ", did);

        for(const auto& result : burst_domain(did))
        {
            const auto& pid2 = result.first;
            const auto& p2   = result.second;
            const auto  fid2 = this->world_->on_which_face(pid2).value();
            const auto  did2 = form_tight_domain_2D(pid2, p2, fid2);

            ECELL4_NGFRD_LOG("bursted particle = ", pid2, " at ", p2.position(),
                             " on ", fid2, " in ", did2);

            const auto dist = ecell4::polygon::distance(this->world_->polygon(),
                std::make_pair(p.position(),  fid),
                std::make_pair(p2.position(), fid2));

            if(dist < (p.radius() + p2.radius()) * SINGLE_SHELL_FACTOR)
            {
                unique_push_back(intruders, did2);
            }
            ECELL4_NGFRD_LOG("resulting particle is at ", dist, " distant");
        }
    }
    ECELL4_NGFRD_LOG("all the min-shell-intruders are bursted.");

    if( ! intruders.empty())
    {
        ECELL4_NGFRD_LOG("there is an intruder: ", intruders,
                         ". trying to form pair or multi.");
        return intruders;
    }

    ECELL4_NGFRD_LOG("No intruders here. draw single shell");

    // dispatch shell shape

    const Real max_circle_size = max_circular_shell_size_at(p.position(), fid);
    if(max_circle_size < min_shell_radius)
    {
        // no intruders, but there is an intrusive vertex. try to form conical one.
        return form_single_conical_domain_2D(pid, p, fid);
    }

    // find maximum possible shell size.

    Real max_shell_size = max_circle_size;
    for(const auto& sidp : shells_.list_shells_within_radius_2D(
            std::make_pair(p.position(), fid), max_circle_size))
    {
        const auto did  = sidp.first.second.domain_id().value();
        const auto& dom = domains_.at(did).second;

        Real dist = sidp.second;
        if(dom.is_single_circular())
        {
            const auto& sh   = dom.as_single_circular().shell();
            const auto& pid2 = dom.as_single_circular().particle_id();
            const auto  p2   = world_->get_particle(pid2).second;

            const Real modest_dist = p.radius() +
                (dist + sh.shape().radius() - p.radius() - p2.radius()) *
                p.D() / (p.D() + p2.D());

            if(min_shell_radius < modest_dist)
            {
                dist = modest_dist;
            }
            else if(dom.as_single_circular().dt() == 0)
            {
                // take care about tight domains to avoid infinite loop...
                dist = dist + sh.shape().radius() - p2.radius() * SINGLE_SHELL_FACTOR;
            }
        }
        max_shell_size = std::min(max_shell_size, dist);
    }
    assert(min_shell_radius <= max_shell_size);

    const auto shell_size       = max_shell_size * SAFETY_SHRINK;
    const auto effective_radius = shell_size - p.radius();
    assert(0.0 < effective_radius);

    ECELL4_NGFRD_LOG("shell size = ", shell_size, ", effective_radius = ", effective_radius);

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

    ECELL4_NGFRD_LOG("2D single circular domain formed!");

    return boost::none;
}

bool NGFRDSimulator::form_pair_domain_2D(
        const ParticleID& pid1, const Particle& p1, const FaceID& fid1,
        const DomainID& intruder_id)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("form_pair_domain_2D: forming domain for particle ", pid1);

    if( ! domains_.at(intruder_id).second.is_single_circular())
    {
        return false;
    }
    const auto& partner = domains_.at(intruder_id).second.as_single_circular();

    // intruders should be a tight domain that is bursted while trying to form a
    // single domain.
    if(partner.dt() != 0.0)
    {
        return false;
    }

    const auto& pid2 = partner.particle_id();
    const auto& p2   = world_->get_particle(pid2).second;
    const auto& fid2 = world_->on_which_face(pid2).value();

    const Real D1  = p1.D();
    const Real D2  = p2.D();
    const Real D12 = D1 + D2;

    const Real r1  = p1.radius();
    const Real r2  = p2.radius();
    const Real r12 = r1 + r2;

    const auto& poly = this->polygon();

    // from p1 t op2
    const Real3 ipv = ecell4::polygon::direction(poly,
            std::make_pair(p1.position(), fid1),
            std::make_pair(p2.position(), fid2));
    const auto  com = ecell4::polygon::travel(poly,
            std::make_pair(p1.position(), fid1), ipv * (D1 / D12));

    const Real ipv_len = length(ipv);
    assert(r12 <= ipv_len);

    const Real min_shell_radius =
        std::max(ipv_len * D1 / D12 + r1, ipv_len * D2 / D12 + r2) * 3;

    Real largest_possible_shell_size = max_circular_shell_size_at(p1.position(), fid1);
    for(const auto& item : this->world_->list_particles_within_radius_2D(
                com, largest_possible_shell_size))
    {
        const auto dist = item.second;
        if(dist < min_shell_radius)
        {
            // it is impossible to form a domain.
            return false;
        }
        else
        {
            largest_possible_shell_size = dist;
        }
    }

    const Real shell_size = largest_possible_shell_size * SAFETY_SHRINK;
    assert(min_shell_radius < shell_size);

    const Real Dipv = D1 + D2;
    const Real Dcom = (D1 * D2) / (D1 + D2);

    // XXX: those shell sizes are too small. To maximize the efficiency, we need
    //      to be a bit more smart
    const Real shell_radius_ipv = (shell_size - std::max(r1, r2)) * Dipv / (Dcom + Dipv);
    const Real shell_radius_com = (shell_size - std::max(r1, r2)) * Dcom / (Dcom + Dipv);

    ECELL4_NGFRD_LOG("shell size   = ", shell_size);
    ECELL4_NGFRD_LOG("CoM boundary = ", shell_radius_com);
    ECELL4_NGFRD_LOG("ipv boundary = ", shell_radius_ipv);
    ECELL4_NGFRD_LOG("ipv half + radius    = ", shell_radius_ipv * D1 / D12 + r1);
    ECELL4_NGFRD_LOG("ipv another + radius = ", shell_radius_ipv * D2 / D12 + r2);

    assert(shell_radius_com + std::max(shell_radius_ipv * D1 / D12 + r1, shell_radius_ipv * D2 / D12 + r2) < shell_size);

    {
        // remove shell
        this->shells_.remove_shell(partner.shell_id());

        // remove domain and event from container

        using std::swap; // zero clear the element in the map to avoid bugs
        std::pair<event_id_type, Domain> evid_dom;
        swap(evid_dom, domains_.at(intruder_id));

        this->domains_.erase(intruder_id);
        this->scheduler_.remove(evid_dom.first);
    }

    const auto dt_reaction1 = draw_single_reaction_time(p1.species());
    const auto dt_reaction2 = draw_single_reaction_time(p2.species());

    const auto rules = this->model_->query_reaction_rules(p1.species(), p2.species());
    const auto k_tot = std::accumulate(rules.begin(), rules.end(), Real(0),
            [](const Real k, const ReactionRule& rule) -> Real {
                return k + rule.k();
            });

    greens_functions::GreensFunction2DAbsSym gf_com(Dcom, shell_radius_com);
    greens_functions::GreensFunction2DRadAbs gf_ipv(Dipv, k_tot, ipv_len, r12, shell_radius_ipv);

    const auto dt_escape_com = gf_com.drawTime(this->world_->rng()->uniform(0.0, 1.0));
    const auto dt_event_ipv  = gf_ipv.drawTime(this->world_->rng()->uniform(0.0, 1.0));

    Real dt;
    PairCircularDomain::EventKind event_kind;
    if(std::min(dt_reaction1, dt_reaction2) < std::min(dt_escape_com, dt_event_ipv))
    {
        if(dt_reaction1 < dt_reaction2)
        {
            dt = dt_reaction1;
            event_kind = PairCircularDomain::EventKind::SingleReaction1;
        }
        else
        {
            dt = dt_reaction2;
            event_kind = PairCircularDomain::EventKind::SingleReaction2;
        }
    }
    else
    {
        if(dt_escape_com < dt_event_ipv)
        {
            dt = dt_escape_com;
            event_kind = PairCircularDomain::EventKind::ComEscape;
        }
        else
        {
            dt = dt_event_ipv;
            event_kind = PairCircularDomain::EventKind::PairEvent;
        }
    }

    const auto did = didgen_();
    const auto sid = sidgen_();

    // construct shell and assign it to shell container
    CircularShell sh(/* thick = */std::max(r1, r2),
            Circle(shell_size, com.first, poly.triangle_at(com.second).normal()), com.second);

    this->shells_.update_shell(sid, Shell(sh, did));

    PairCircularDomain dom(event_kind, dt, world_->t(), sid, sh,
            shell_radius_com, shell_radius_ipv, pid1, pid2, ipv,
            std::move(gf_com), std::move(gf_ipv));

    // add event with the same domain ID
    const auto evid = this->scheduler_.add(
            std::make_shared<event_type>(this->t() + dom.dt(), did));

    // update begin_time and re-insert domain into domains_ container
    this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));

    // pair is formed.
    return true;
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

    ECELL4_NGFRD_LOG("forming Multi domain for 2D particle ", pid);

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

    // If there are several domains within the min_shell_radius, burst them and
    // list fatal intruders (particles within min_shell)
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

        ECELL4_NGFRD_LOG("nearest one is at ", dist, " distant.");

        // To maximize the efficiency, it is better to consider the neighbor's
        // diffusion coefficient because that makes the shell size best balance.
        // But in some cases, the `modest_radius` becomes less than the
        // min_shell_radius. This is unacceptable. So we need to check it and
        // use min_shell_radius.
        //     Of course, if two particles are closer than the sum of their
        // min_shell_radius, then it is the time for Multi. But we already
        // checked shells that are within the min shell radius.
        const auto modest_radius = p.radius() + (dist - p.radius()) * p.D() /
                                   (p.D() + neighbor.D());
        if(min_shell_radius < modest_radius)
        {
            max_radius = std::min(max_radius, modest_radius);
            ECELL4_NGFRD_LOG("the radius is bound by a particle of which ID is ",
                    nearest_particle.front().first.first,
                    ". The possible max radius is ", max_radius);
        }
        else
        {
            max_radius = min_shell_radius;
        }
    }
    if( ! nearest_face.empty())
    {
        max_radius = std::min(nearest_face.front().second - largest_2D_particle,
                              max_radius);
        ECELL4_NGFRD_LOG("the radius is bound by a face of which ID is ",
                nearest_face.front().first.first, ". The possible max radius is ", max_radius);
    }
    ECELL4_NGFRD_LOG("the candidate of shell radius is ", max_radius);

    // check other 3D shells within max_radius.
    // Note: Here, we already subtract largest_2D_particle that is effective
    //       thickness of the 2D shells. So we just ignore 2D shells here, and
    //       only check 3D shells.
    for(const auto& item : this->shells_.list_shells_within_radius_3D(
                p.position(), max_radius))
    {
        max_radius = std::min(max_radius, item.second);
    }
    ECELL4_NGFRD_LOG("After checking intruders, the max possible radius = ", max_radius);

    // after SAFETY_SHRINK, the shell size could be (slightly) smaller than the
    // min_shell_radius, but I don't think it is a problem.
    ensure(min_shell_radius <= max_radius,
           "min_shell_radius (", min_shell_radius, ") should be smaller than ",
           "the max radius (", max_radius, ").");

    const auto shell_size = max_radius * SAFETY_SHRINK;

    // paranoiac check: no shells are within the current shell size.
    ensure(this->shells_.list_shells_within_radius_3D(p.position(), shell_size).empty());

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
        const ParticleID& pid1, const Particle& p1, const DomainID& intruder_id)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("form_pair_domain_3D: forming domain for particle ", pid1);

    // intruder is a single spherical domain
    if( ! domains_.at(intruder_id).second.is_single_spherical())
    {
        return false;
    }
    const auto& partner = domains_.at(intruder_id).second.as_single_spherical();

    // intruders should be a tight domain that is bursted while forming single,
    // to form a pair domain.
    if(partner.dt() != 0.0)
    {
        return false;
    }

    const auto& pid2 = partner.particle_id();
    const auto& p2   = world_->get_particle(pid2).second;

    const Real D1  = p1.D();
    const Real D2  = p2.D();
    const Real D12 = D1 + D2;

    const Real r1  = p1.radius();
    const Real r2  = p2.radius();
    const Real r12 = r1 + r2;

    const auto pbc = this->world_->boundary();

    const Real3 com = pbc.apply_boundary(p1.position() * (D2 / D12) +
         pbc.periodic_transpose(p2.position(), p1.position()) * (D1 / D12));

    // from p1 -> p2
    const Real3    ipv = pbc.periodic_transpose(p2.position(), p1.position()) - p1.position();
    const Real ipv_len = length(ipv);

    const Real min_shell_size = std::max(ipv_len * D1 / D12 + r1 * 3,
                                         ipv_len * D2 / D12 + r2 * 3);

    // we first check if there is any particles around the pair because
    // there can be particles that is bursted while trying to form a single
    // with high probability.
    const boost::container::static_vector<std::pair<std::pair<ParticleID, Particle>, Real>, 1>
         nearest_particles = world_->nearest_particle_3D(com, pid1, pid2);
    const auto& nearest_particle = nearest_particles.at(0).first.second;

    const auto dnearest = length(com -
            pbc.periodic_transpose(nearest_particle.position(), com)) -
            nearest_particle.radius() * SINGLE_SHELL_FACTOR * SAFETY_EXPAND;
    if(dnearest < min_shell_size)
    {
        // there is a particle locating too close to the center of the pair.
        return false;
    }
    Real largest_possible_shell_size = dnearest;

    // then check if there is a intrusive domain or a polygon face.

    const Real largest_2D_particle = world_->largest_particle_radius_2D();
    for(const auto& item : this->world_->list_faces_within_radius(
                com, largest_possible_shell_size + largest_2D_particle))
    {
        const auto dist = item.second;
        if(dist < min_shell_size)
        {
            // it is impossible to form a domain.
            return false;
        }
        else
        {
            largest_possible_shell_size = dist;
        }
    }

    for(const auto& item : shells_.list_shells_within_radius_3D(
                com, largest_possible_shell_size, partner.shell_id()))
    {
        const auto dist = item.second;
        if(dist < min_shell_size)
        {
            // it is impossible to form a domain.
            return false;
        }
        else
        {
            largest_possible_shell_size = dist;
        }
    }

    const Real shell_size = largest_possible_shell_size * SAFETY_SHRINK;

    // we have already tried to form a single, so we are sure that pair is worth
    // forming compared to two singles. And we assume that pair is always faster
    // than a multi.

    // -----------------------------------------------------------------------
    // construct a pair domain.

    // remove partner domain because it is no longer needed after forming pair.
    {
        // remove shell
        this->shells_.remove_shell(partner.shell_id());

        // remove domain and event from container

        using std::swap; // zero clear the element in the map to avoid bugs
        std::pair<event_id_type, Domain> evid_dom;
        swap(evid_dom, domains_.at(intruder_id));

        this->domains_.erase(intruder_id);
        this->scheduler_.remove(evid_dom.first);
    }

    //              com,  ipv
    const std::pair<Real, Real> boundaries =
        [](Real r1, Real r2, Real D1, Real D2, const Real D12,
           const Real ipv_len, const Real shell_size, const Real SAFETY_SHRINK)
        {
            const auto D_geom = std::sqrt(D1 * D2);
            if((D_geom - D1) * ipv_len / D12 + shell_size + std::sqrt(D1 / D2) * (r2 - shell_size) - r1 < 0)
            {
                std::swap(r1, r2);
                std::swap(D1, D2);
            }

            const Real a_com =
                D_geom * (D1 * (shell_size - r2) + D2 * (shell_size - ipv_len - r2)) /
                (D2 * D2 + D1 * D2 + D_geom * D12);
            const Real a_ipv =
                (D_geom * ipv_len + D12 * (shell_size - r2)) / (D2 + D_geom);

            return std::make_pair(a_com * SAFETY_SHRINK, a_ipv * SAFETY_SHRINK);

        }(r1, r2, D1, D2, D12, ipv_len, shell_size, SAFETY_SHRINK);

    ECELL4_NGFRD_LOG("shell size   = ", shell_size);
    ECELL4_NGFRD_LOG("CoM boundary = ", boundaries.first);
    ECELL4_NGFRD_LOG("ipv boundary = ", boundaries.second);
    ECELL4_NGFRD_LOG("ipv half + radius    = ", boundaries.second * D1 / D12 + r1);
    ECELL4_NGFRD_LOG("ipv another + radius = ", boundaries.second * D2 / D12 + r2);

    assert(boundaries.first + boundaries.second * D1 / D12 + r1 <= shell_size);
    assert(boundaries.first + boundaries.second * D2 / D12 + r2 <= shell_size);

    const auto dt_reaction1 = draw_single_reaction_time(p1.species());
    const auto dt_reaction2 = draw_single_reaction_time(p2.species());

    const auto rules = this->model_->query_reaction_rules(p1.species(), p2.species());
    const auto k_tot = std::accumulate(rules.begin(), rules.end(), Real(0),
            [](const Real k, const ReactionRule& rule) -> Real {
                return k + rule.k();
            });

    greens_functions::GreensFunction3DAbsSym gf_com(D1 * D2 / D12, boundaries.first);
    greens_functions::GreensFunction3DRadAbs gf_ipv(D12, k_tot, ipv_len, r12, boundaries.second);

    const auto dt_escape_com = gf_com.drawTime(this->world_->rng()->uniform(0.0, 1.0));
    const auto dt_event_ipv  = gf_ipv.drawTime(this->world_->rng()->uniform(0.0, 1.0));

    Real dt;
    PairSphericalDomain::EventKind event_kind;

    if(std::min(dt_reaction1, dt_reaction2) < std::min(dt_escape_com, dt_event_ipv))
    {
        if(dt_reaction1 < dt_reaction2)
        {
            dt = dt_reaction1;
            event_kind = PairSphericalDomain::EventKind::SingleReaction1;
        }
        else
        {
            dt = dt_reaction2;
            event_kind = PairSphericalDomain::EventKind::SingleReaction2;
        }
    }
    else
    {
        if(dt_escape_com < dt_event_ipv)
        {
            dt = dt_escape_com;
            event_kind = PairSphericalDomain::EventKind::ComEscape;
        }
        else
        {
            dt = dt_event_ipv;
            event_kind = PairSphericalDomain::EventKind::PairEvent;
        }
    }

    const auto did = didgen_();
    const auto sid = sidgen_();

    // construct shell and assign it to shell container
    SphericalShell sh(Sphere(com, shell_size));
    this->shells_.update_shell(sid, Shell(sh, did));

    PairSphericalDomain dom(event_kind, dt, world_->t(), sid, sh,
            boundaries.first, boundaries.second, pid1, pid2, ipv,
            std::move(gf_com), std::move(gf_ipv));

    // add event with the same domain ID
    const auto evid = this->scheduler_.add(
            std::make_shared<event_type>(this->t() + dom.dt(), did));

    // update begin_time and re-insert domain into domains_ container
    this->domains_[did] = std::make_pair(evid, Domain(std::move(dom)));

    // pair is formed.
    return true;
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
    ECELL4_NGFRD_LOG("firing single circular: ", did, " at t = ", this->t());
    ECELL4_NGFRD_LOG("included shell: ", dom.shell_id(), ", particle: ", dom.particle_id());

    if(dom.dt() == 0.0)
    {
        ECELL4_NGFRD_LOG("This domain is a tight domain. Fire without move.");
        this->shells_.remove_shell(dom.shell_id());
        return boost::container::small_vector<std::pair<ParticleID, Particle>, 4>{
                world_->get_particle(dom.particle_id())
            };
    }

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions;
    SingleCircularPropagator prop(did,
            *(this->model_), *(this->world_), *this, *(this->world_->rng()),
            SINGLE_CIRCULAR_MAX_RETRY, last_reactions);

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
NGFRDSimulator::fire_single_conical(const DomainID& did, SingleConicalDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("firing single conical: ", did, " at t = ", this->t());
    ECELL4_NGFRD_LOG("included shell: ", dom.shell_id(), ", particle: ", dom.particle_id());

    // conical never be a tight domain.

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions;
    SingleConicalPropagator prop(did,
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
    ECELL4_NGFRD_LOG("firing single spherical: ", did, " at t = ", this->t());
    ECELL4_NGFRD_LOG("included shell: ", dom.shell_id(), ", particle: ", dom.particle_id());

    if(dom.dt() == 0.0) // means it is a tight domain
    {
        ECELL4_NGFRD_LOG("This domain is a tight domain. Fire without move.");
        this->shells_.remove_shell(dom.shell_id());
        return boost::container::small_vector<std::pair<ParticleID, Particle>, 4>{
                world_->get_particle(dom.particle_id())
            };
    }

    ECELL4_NGFRD_LOG("Propagating for ", dom.dt());
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
NGFRDSimulator::fire_pair_spherical(const DomainID& did, PairSphericalDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("firing pair spherical: ", did, " at t = ", this->t());
    ECELL4_NGFRD_LOG("included shell: ", dom.shell_id(),
                     ", particle1: ", dom.particle1_id(),
                     ", particle2: ", dom.particle2_id());

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions;
    PairSphericalPropagator prop(
            did, *(this->model_), *(this->world_), *this,
            *(this->world_->rng()), PAIR_SPHERICAL_MAX_RETRY, last_reactions);

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
NGFRDSimulator::fire_pair_circular(const DomainID& did, PairCircularDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("firing pair circular: ", did, " at t = ", this->t());
    ECELL4_NGFRD_LOG("included shell: ", dom.shell_id(),
                     ", particle1: ", dom.particle1_id(),
                     ", particle2: ", dom.particle2_id());

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions;
    PairCircularPropagator prop(
            did, *(this->model_), *(this->world_), *this,
            *(this->world_->rng()), PAIR_CIRCULAR_MAX_RETRY, last_reactions);

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
    ECELL4_NGFRD_LOG("firing multi: ", did, " at t = ", this->t());
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
        ECELL4_NGFRD_LOG("It is a tight shell. removing without movement...");
        this->shells_.remove_shell(dom.shell_id());
        return boost::container::small_vector<std::pair<ParticleID, Particle>, 4>{
                world_->get_particle(dom.particle_id())
            };
    }
    if(this->t() == dom.begin_time())
    {
        ECELL4_NGFRD_LOG("Now is the time when this domain is formed. "
                         "removing without movement...");
        this->shells_.remove_shell(dom.shell_id());
        return boost::container::small_vector<std::pair<ParticleID, Particle>, 4>{
                world_->get_particle(dom.particle_id())
            };
    }
    ECELL4_NGFRD_LOG("Propagating...");

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions;
    SingleCircularPropagator prop(did,
            *(this->model_), *(this->world_), *this, *(this->world_->rng()),
            SINGLE_CIRCULAR_MAX_RETRY, last_reactions);

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
NGFRDSimulator::burst_single_conical(const DomainID& did, SingleConicalDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("bursting single conical: ", did);
    ECELL4_NGFRD_LOG("included shell: ", dom.shell_id());

    // conical shell never be a tight domain.

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions;
    SingleConicalPropagator prop(did,
            *(this->model_), *(this->world_), *this, *(this->world_->rng()),
            SINGLE_CONICAL_MAX_RETRY, last_reactions);

    // While bursting single domain, it never reacts (becuase the time when
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
        ECELL4_NGFRD_LOG("It is a tight shell. removing without movement...");
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
NGFRDSimulator::burst_pair_spherical(const DomainID& did, PairSphericalDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("bursting pair spherical: ", did);
    ECELL4_NGFRD_LOG("included shell: ", dom.shell_id());

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions;
    PairSphericalPropagator prop(
            did, *(this->model_), *(this->world_), *this,
            *(this->world_->rng()), PAIR_SPHERICAL_MAX_RETRY, last_reactions);

    const std::array<ParticleID, 2> resulting_particles = prop.burst(dom, this->t());

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4> results;
    results.push_back(world_->get_particle(resulting_particles[0]));
    results.push_back(world_->get_particle(resulting_particles[1]));

    // remove shell from shell container
    this->shells_.remove_shell(dom.shell_id());

    assert(last_reactions.empty());
    return results;
}

boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
NGFRDSimulator::burst_pair_circular(const DomainID& did, PairCircularDomain dom)
{
    ECELL4_NGFRD_LOG_FUNCTION();
    ECELL4_NGFRD_LOG("bursting pair circular: ", did);
    ECELL4_NGFRD_LOG("included shell: ", dom.shell_id());

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions;
    PairCircularPropagator prop(
            did, *(this->model_), *(this->world_), *this,
            *(this->world_->rng()), PAIR_SPHERICAL_MAX_RETRY, last_reactions);

    const std::array<ParticleID, 2> resulting_particles = prop.burst(dom, this->t());

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4> results;
    results.push_back(world_->get_particle(resulting_particles[0]));
    results.push_back(world_->get_particle(resulting_particles[1]));

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
