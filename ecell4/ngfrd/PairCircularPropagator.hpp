#ifndef ECELL4_NGFRD_PAIR_CIRCULAR_PROPAGATOR_HPP
#define ECELL4_NGFRD_PAIR_CIRCULAR_PROPAGATOR_HPP
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/ngfrd/ReactionInfo.hpp>
#include <ecell4/ngfrd/PairCircularDomain.hpp>
#include <ecell4/ngfrd/NGFRDWorld.hpp>

#include <boost/container/static_vector.hpp>

namespace ecell4
{
namespace ngfrd
{

class NGFRDSimulator; // forward decl

class PairCircularPropagator
{
  public:
    PairCircularPropagator(const DomainID& self,
        const Model& model, NGFRDWorld& world, NGFRDSimulator& sim,
        RandomNumberGenerator& rng, const std::size_t max_retry,
        std::vector<std::pair<ReactionRule, ReactionInfo>>& last_reactions)
        : model_(model), world_(world), sim_(sim), rng_(rng), self_id_(self),
          max_retry_count_(max_retry), last_reactions_(last_reactions)
    {}

    boost::container::static_vector<ParticleID, 3>
    operator()(const PairCircularDomain& dom)
    {
        switch(dom.eventkind())
        {
            case PairCircularDomain::EventKind::ComEscape:
            {
                return this->com_escape(dom);
            }
            case PairCircularDomain::EventKind::IpvEscape:
            {
                return this->ipv_escape(dom);
            }
            case PairCircularDomain::EventKind::SingleReaction1:
            {
                // it passes particle IDs to determine which particle reacts.
                return this->single_reaction(dom, dom.particle1_id(), dom.particle2_id());
            }
            case PairCircularDomain::EventKind::SingleReaction2:
            {
                return this->single_reaction(dom, dom.particle2_id(), dom.particle1_id());
            }
            case PairCircularDomain::EventKind::PairReaction:
            {
                return this->pair_reaction(dom);
            }
            case PairCircularDomain::EventKind::PairEvent:
            {
                const auto& gf_ipv = dom.gf_ipv();
                const auto event_kind = gf_ipv.drawEventType(rng_.uniform(0.0, 1.0), dom.dt());

                if(event_kind == greens_functions::GreensFunction::IV_ESCAPE)
                {
                    return this->ipv_escape(dom);
                }
                else if(event_kind == greens_functions::GreensFunction::IV_REACTION)
                {
                    return this->pair_reaction(dom);
                }
                else
                {
                    throw_exception<IllegalState>("ecell4::ngfrd::PairCircular"
                        "Propagator: invalid greens_functions::EventKind");
                }
            }
            default:
            {
                throw_exception<IllegalState>("ecell4::ngfrd::PairCircular"
                    "Propagator: invalid eventkind");
            }
        }
    }

    // bursting occurs before any event happens. So here we know that the resulting
    // particles IDs are the same as the pair included in this domain.
    std::array<ParticleID, 2> burst(const PairCircularDomain& dom, const Real time_now)
    {
        return this->propagate(dom, time_now - dom.begin_time());
    }

    std::size_t max_retry_count() const noexcept {return max_retry_count_;}

    std::size_t get_rejected_move_count() const noexcept
    {
        return rejected_move_count_;
    }
    std::vector<std::pair<ReactionRule, ReactionInfo>>& last_reactions() const noexcept
    {
        return last_reactions_;
    }

  private:

    boost::container::static_vector<ParticleID, 3>
    com_escape(const PairCircularDomain& dom);
    boost::container::static_vector<ParticleID, 3>
    ipv_escape(const PairCircularDomain& dom);

    boost::container::static_vector<ParticleID, 3>
    single_reaction(const PairCircularDomain& dom,
                    const ParticleID& reactant, const ParticleID& other);

    boost::container::static_vector<ParticleID, 3>
    pair_reaction(const PairCircularDomain& dom);

    boost::container::static_vector<ParticleID, 3>
    attempt_1to1_reaction(const PairCircularDomain& dom,
            const ParticleID&, const Particle&, const ReactionRule& rule,
            const ParticleID& other);
    boost::container::static_vector<ParticleID, 3>
    attempt_1to2_reaction(const PairCircularDomain& dom,
            const ParticleID&, const Particle&, const ReactionRule& rule,
            const ParticleID& other);

    boost::container::static_vector<ParticleID, 3>
    attempt_2to1_reaction(const PairCircularDomain& dom,
            const ParticleID&, const Particle&, const FaceID&,
            const ParticleID&, const Particle&, const FaceID&,
            const ReactionRule& rule, const std::pair<Real3, FaceID>& com);

    std::array<ParticleID, 2> propagate(const PairCircularDomain& dom, const Real dt);

    ReactionRule const& determine_single_reaction_rule(const Species&, const Real rnd);
    ReactionRule const& determine_pair_reaction_rule(const Species&, const Species&, const Real rnd);

    bool is_inside_of_shell(const PairCircularDomain& dom,
            const std::pair<Real3, FaceID>& pos, const Real radius) const;

    Real3 generate_ipv(const FaceID& fid, const Real3& ipv, const Real r, const Real theta) const;
    Real3 draw_2D_displacement(const FaceID& fid, const Real len);

  private:

    const Model&           model_;
    NGFRDWorld&            world_;
    NGFRDSimulator&        sim_;
    RandomNumberGenerator& rng_;
    DomainID               self_id_;
    std::size_t            max_retry_count_;
    std::size_t            rejected_move_count_;
    std::vector<std::pair<ReactionRule, ReactionInfo>>& last_reactions_;
};

} // ngfrd
} // ecell4
#endif// ECELL4_NGFRD_PAIR_CIRCULAR_PROPAGATOR_HPP
