#ifndef ECELL4_NGFRD_NGFRDSIMULATOR_HPP
#define ECELL4_NGFRD_NGFRDSIMULATOR_HPP

#include <boost/format.hpp>
#include <boost/none_t.hpp>
#include <boost/optional.hpp>
#include <boost/variant.hpp>
#include <boost/container/small_vector.hpp>

#include <ecell4/core/Model.hpp>
#include <ecell4/core/EventScheduler.hpp>
#include <ecell4/core/SerialIDGenerator.hpp>
#include <ecell4/core/SimulatorBase.hpp>

#include <ecell4/ngfrd/ShellID.hpp>
#include <ecell4/ngfrd/Shell.hpp>
#include <ecell4/ngfrd/DomainID.hpp>
#include <ecell4/ngfrd/Domain.hpp>
#include <ecell4/ngfrd/ShellContainer.hpp>
#include <ecell4/ngfrd/NGFRDEvent.hpp>
#include <ecell4/ngfrd/NGFRDWorld.hpp>

#include <ecell4/ngfrd/Logger.hpp>

#include <greens_functions/PairGreensFunction.hpp>
#include <greens_functions/GreensFunction3DRadAbs.hpp>
#include <greens_functions/GreensFunction3DRadInf.hpp>
#include <greens_functions/GreensFunction3DAbsSym.hpp>
#include <greens_functions/GreensFunction3DAbs.hpp>
#include <greens_functions/GreensFunction3D.hpp>

namespace ecell4
{
namespace ngfrd
{

class NGFRDSimulator final: public SimulatorBase<NGFRDWorld>
{
public:
    static constexpr Real SAFETY              = 1.0 + 1e-5;
    static constexpr Real SINGLE_SHELL_FACTOR = 1.0 + 0.1;
    static constexpr Real MULTI_SHELL_FACTOR  = 1.0 + 0.05;
    static constexpr Real DEFAULT_DT_FACTOR   = 1e-5;
    static constexpr Real CUTOFF_FACTOR       = 5.6;

    static constexpr std::size_t SINGLE_3D_MAX_RETRY = 4;


public:

    using base_type       = SimulatorBase<NGFRDWorld>;
    using model_type      = typename base_type::model_type;
    using world_type      = typename base_type::world_type;
    using event_type      = NGFRDEvent;
    using scheduler_type  = EventSchedulerBase<event_type>;
    using event_id_type   = typename scheduler_type::identifier_type;

    NGFRDSimulator(const std::shared_ptr<world_type>& world,
                   const std::shared_ptr<model_type>& model,
                   const Real bd_dt_factor_3D = 1e-5,
                   const Real bd_dt_factor_2D = 1e-3,
                   const Real reaction_length = 1e-1,
                   const std::size_t max_retry_moves = 1)
        : base_type(world, model), dt_factor_3D_(bd_dt_factor_3D),
          dt_factor_2D_(bd_dt_factor_2D), reaction_length_(reaction_length),
          max_retry_(max_retry_moves),
          shells_(world->edge_lengths(), world->polygon_ptr()),
          is_uninitialized_(true)
    {}
    ~NGFRDSimulator() override = default;

    void initialize() override
    {
        ECELL4_NGFRD_LOG_FUNCTION();
        // --------------------------------------------------------------------
        // clear everything

        scheduler_.clear();
        domains_  .clear();
        shells_   .clear();
        shells_.reset_boundary(this->world_->edge_lengths());
        last_reactions_.clear();

        // --------------------------------------------------------------------
        // form domain for all particles

        for(const auto& pp : world_->particles())
        {
            ECELL4_NGFRD_LOG("initializing domain for ", pp.first);
            form_domain(pp.first, pp.second);
        }

        // --------------------------------------------------------------------
        // form birth domain if needed

        // TODO!

        this->is_uninitialized_ = false;
        return;
    }
    void step() override
    {
        if(this->is_uninitialized_)
        {
            this->initialize();
        }
        this->step_unchecked();
        return;
    }

    bool step(const Real& upto) override
    {
        if(this->is_uninitialized_)
        {
            this->initialize();
        }

        if(upto <= this->t())
        {
            return false;
        }

        if(scheduler_.next_time() <= upto)
        {
            this->step_unchecked();
            return true;
        }
        else
        {
            this->set_t(upto);
            this->finalize();
            return false;
        }
    }
    void finalize()
    {
        ECELL4_NGFRD_LOG_FUNCTION();
        std::vector<DomainID> non_singles;

        // Single does not stick out from shell when it is bursted.
        // So burst single first.
        // Events are cleared later, so we don't need to pop it in the loop
        for(const auto& eidp: scheduler_.events())
        {
            const auto& eid = eidp.first;
            const auto& ev  = eidp.second;
            const auto& did = ev->domain_id();

            if(!this->domains_.at(did).second.is_single_sphecal())
            {
                non_singles.push_back(ev->domain_id());
                continue;
            }

            // move out the domain from domtain conteiner
            auto dom = std::move(this->domains_.at(did).second);
            this->domains_.erase(did); // remove the element from container

            this->burst_domain(did, std::move(dom));

            assert(this->diagnosis());
        }

        // then burst non-single domains.
        for(const auto& did : non_singles)
        {
            auto dom = std::move(this->domains_.at(did).second);
            this->domains_.erase(did);

            this->burst_domain(did, std::move(dom));
            assert(this->diagnosis());
        }
        assert(this->domains_.empty());

        this->dt_ = 0.0;
        this->is_uninitialized_ = true;

        scheduler_.clear();
        return;
    }

    Real next_event_time() const
    {
        return scheduler_.next_time();
    }
    Real dt() const override
    {
        return scheduler_.next_time() - world_->t();
    }

    std::shared_ptr<RandomNumberGenerator> const& rng() const noexcept
    {
        return this->world_->rng();
    }

    // shrink domains that overlap with the sphere centered at pos with radius.
    // Multi domains are kept intact.
    void determine_positions_3D(const Real3& pos, const Real radius)
    {
        this->determine_positions_3D_impl(pos, radius,
                [](const DomainID&){return false;});
        return;
    }
    void determine_positions_3D(const Real3& pos, const Real radius,
                                const DomainID& ignored)
    {
        this->determine_positions_3D_impl(pos, radius,
                [&](const DomainID& did){return did == ignored;});
        return;
    }

    void determine_positions_2D(const std::pair<Real3, FaceID>& pos,
                                const Real radius)
    {
        // XXX: currently, all the domains are multi.
        return;
    }
    void determine_positions_2D(const std::pair<Real3, FaceID>& pos,
                                const Real radius, const DomainID& ignored)
    {
        // XXX: currently, all the domains are multi.
        return;
    }


    // health check
    bool diagnosis() const
    {
        return shells_.diagnosis();
    }

private:

    void step_unchecked()
    {
        ECELL4_NGFRD_LOG_FUNCTION();
        this->num_steps_ += 1;

        if(scheduler_.size() == 0)
        {
            ECELL4_NGFRD_LOG("no event found");
            this->set_t(scheduler_.next_time());
            return;
        }

        const auto eidp = scheduler_.pop();
        this->set_t(eidp.second->time());

        const auto fired = fire_event(*eidp.second);

        for(const auto& pidp : fired)
        {
            ECELL4_NGFRD_LOG("forming domain for resulting particle ", pidp.first);
            this->form_domain(pidp.first, pidp.second);
        }

        const auto next_time = scheduler_.top().second->time();
        this->dt_ = next_time - this->t();

        if(this->dt_ == 0.0)
        {
            this->zero_step_count_ += 1;
        }
        return;
    }

    // -----------------------------------------------------------------------
    // form_domain

    void form_domain(const ParticleID& pid, const Particle& p)
    {
        if(const auto fid = this->world_->on_which_face(pid))
        {
            this->form_domain_2D(pid, p, *fid);
        }
        else
        {
            this->form_domain_3D(pid, p);
        }
    }

    void form_domain_2D(const ParticleID& pid, const Particle& p, const FaceID& fid);
    void form_domain_3D(const ParticleID& pid, const Particle& p);

    void form_tight_domain_2D(const ParticleID& pid, const Particle& p, const FaceID& fid);
    void form_tight_domain_3D(const ParticleID& pid, const Particle& p);

    // -----------------------------------------------------------------------
    // fire_event
    // - it checks if domain knows its event
    // - it moves out the domain from domain container

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
    fire_event(const NGFRDEvent& ev)
    {
        // pop domain from domains_ container
        auto didp_iter = domains_.find(ev.domain_id());
        assert(ev.domain_id() == didp_iter->first);

        auto dom = std::move(didp_iter->second);
        domains_.erase(didp_iter);
        return fire_domain(ev.domain_id(), std::move(dom.second));
    }

    // fire_domain
    // - propagates the particle
    // - it removes shells if appropreate

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
    fire_domain(const DomainID& did, Domain dom)
    {
        assert(domains_.count(did) == 0); // the domain is moved out from container
        switch(dom.kind())
        {
            case Domain::DomainKind::SingleSpherical:
            {
                return this->fire_single_spherical(did, std::move(dom.as_single_spherical()));
            }
            case Domain::DomainKind::Multi:
            {
                return this->fire_multi(did, std::move(dom.as_multi()));
            }
            default:
            {
                throw_exception<NotImplemented>("NGFRD::fire_domain: unknown "
                        "domain kind (", static_cast<int>(dom.kind()), ").");
            }
        }
    }

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
    fire_single_spherical(const DomainID& did, SingleSphericalDomain dom);

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
    fire_multi(const DomainID& did, MultiDomain dom);

    // -----------------------------------------------------------------------
    // burst_domain
    // - It propagates domains until world_->t(), regardless of the original dt.
    // - It returns the particles after propagation.
    //   - We will keep it if this simulator is finalized.
    //   - We need to add a new domain otherwise.

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
    burst_domain(const DomainID& did, Domain dom)
    {
        assert(domains_.count(did) == 0); // the domain is moved out from container
        switch(dom.kind())
        {
            case Domain::DomainKind::Multi:
            {
                return this->burst_multi(did, std::move(dom.as_multi()));
            }
            default:
            {
                throw_exception<NotImplemented>("NGFRD::burst_domain: unknown "
                        "domain kind (", static_cast<int>(dom.kind()), ").");
            }
        }
    }

    boost::container::small_vector<std::pair<ParticleID, Particle>, 4>
    burst_multi(const DomainID& did, MultiDomain dom)
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

    Polygon const& polygon() const noexcept {return this->world_->polygon();}

private:

    // -----------------------------------------------------------------------
    // determine_positions
    // - It propagates domains until world_->t(), regardless of the original dt.
    // - It automatically adds a new domain to the resulting particles.
    //   - The original domain (that could be troublesome) is removed, and
    //   - Newly formed (less trouble) domain is added
    template<typename Filter>
    void determine_positions_3D_impl(
            const Real3& pos, const Real radius, Filter filter)
    {
        ECELL4_NGFRD_LOG_FUNCTION();

        // fetch shells that overlaps with this region
        boost::container::small_vector<DomainID, 8> dids;
        for(const auto& shd : shells_.list_shells_within_radius_3D(pos, radius))
        {
            if(const auto did = shd.first.second.domain_id())
            {
                if(filter(*did)) {continue;}

                if(std::find(dids.begin(), dids.end(), *did) != dids.end())
                {
                    continue; // already assigned. It may be Multi.
                }
                dids.push_back(*did);
            }
            else
            {
                throw_exception<IllegalState>("shell ", shd.first.first, ": ",
                        shd.first.second, " does not know its own DomainID.");
            }
        }

        // shrink all the domains that have listed
        for(const auto& did : dids)
        {
            if(domains_.at(did).second.is_multi())
            {
                continue;
            }
            // burst the domain
            // - propagate particle             via burst_domain()
            // - remove shell  from shells_     via burst_domain()
            // - remove domain from domains_    by myself
            // - remove event  from scheduler_  by myself

            const auto evid = this->domains_.at(did).first;
            this->scheduler_.remove(evid);

            auto dom = std::move(this->domains_.at(did).second);
            this->domains_.erase(did);

            const auto results = burst_domain(did, std::move(dom));

            // re-wrap the particles by tight domains
            for(const auto& pidp : results)
            {
                this->form_tight_domain_3D(pidp.first, pidp.second);
            }
        }
        return;
    }

private:

    // ------------------------------------------------------------------------
    // inherited from SimulatorBase
//     std::shared_ptr<world_type> world_;
//     std::shared_ptr<model_type> model_;
//     Integer num_steps_;

    // ------------------------------------------------------------------------
    // parameters
    Real dt_factor_3D_;
    Real dt_factor_2D_;
    Real reaction_length_;
    std::size_t max_retry_;

    // ------------------------------------------------------------------------
    // domains and shells
    SerialIDGenerator<DomainID> didgen_;
    std::unordered_map<DomainID, std::pair<event_id_type, Domain>> domains_;

    SerialIDGenerator<ShellID> sidgen_;
    ShellContainer             shells_;

    // ------------------------------------------------------------------------
    // events
    scheduler_type scheduler_;
    Real dt_;
    bool is_uninitialized_;

    // -------------------------------------------------------------------------
    // statistics

    std::size_t rejected_moves_;
    std::size_t zero_step_count_;

    std::vector<std::pair<ReactionRule, ReactionInfo>> last_reactions_;

//     std::array<std::size_t, SingleDomain::EventKinds> single_step_count_;
//     std::array<std::size_t, PairDomain::EventKinds>   pair_step_count_;
//     std::array<std::size_t, MultiDomain::EventKinds>  multi_step_count_;
//     std::array<std::size_t, NUM_DOMAIN_KINDS>         domain_count_per_type_;
};
#undef CHECK

} // egfrd
} // ecell4
#endif /* EGFRDSIMULATOR_HPP */
