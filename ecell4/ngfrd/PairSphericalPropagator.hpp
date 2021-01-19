#ifndef ECELL4_NGFRD_PAIR_SPHERICAL_PROPAGATOR_HPP
#define ECELL4_NGFRD_PAIR_SPHERICAL_PROPAGATOR_HPP
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/ngfrd/ReactionInfo.hpp>
#include <ecell4/ngfrd/PairSphericalDomain.hpp>
#include <ecell4/ngfrd/NGFRDWorld.hpp>

#include <boost/container/static_vector.hpp>

namespace ecell4
{
namespace ngfrd
{

class NGFRDSimulator; // forward decl

class PairSphericalPropagator
{
  public:
    PairSphericalPropagator(const DomainID& self,
        const Model& model, NGFRDWorld& world, NGFRDSimulator& sim,
        RandomNumberGenerator& rng, const std::size_t max_retry,
        std::vector<std::pair<ReactionRule, ReactionInfo>>& last_reactions)
        : model_(model), world_(world), sim_(sim), rng_(rng), self_id_(self),
          max_retry_count_(max_retry), last_reactions_(last_reactions)
    {}

    boost::container::static_vector<ParticleID, 3>
    operator()(const PairSphericalDomain& dom)
    {
        switch(dom.eventkind())
        {
            case PairSphericalDomain::EventKind::ComEscape:
            {
                return this->com_escape(dom);
            }
            case PairSphericalDomain::EventKind::IpvEscape:
            {
                return this->ipv_escape(dom);
            }
            case PairSphericalDomain::EventKind::SingleReaction1:
            {
                return this->single_reaction_1(dom);
            }
            case PairSphericalDomain::EventKind::SingleReaction2:
            {
                return this->single_reaction_2(dom);
            }
            case PairSphericalDomain::EventKind::PairReaction:
            {
                return this->pair_reaction(dom);
            }
            default:
            {
                throw_exception<IllegalState>("ecell4::ngfrd::PairSpherical"
                    "Propagator: invalid eventkind");
            }
        }
    }

    // bursting occurs before any event happens. So here we know that the resulting
    // particles IDs are the same as the pair included in this domain.
    std::array<ParticleID, 2> burst(const PairSphericalDomain& dom, const Real time_now)
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
    com_escape(const PairSphericalDomain& dom);
    boost::container::static_vector<ParticleID, 3>
    ipv_escape(const PairSphericalDomain& dom);

    boost::container::static_vector<ParticleID, 3>
    single_reaction_1(const PairSphericalDomain& dom);
    boost::container::static_vector<ParticleID, 3>
    single_reaction_2(const PairSphericalDomain& dom);

    boost::container::static_vector<ParticleID, 3>
    pair_reaction(const PairSphericalDomain& dom);

    boost::container::static_vector<ParticleID, 3>
    attempt_1to1_reaction(const PairSphericalDomain& dom,
            const ParticleID&, const Particle&, const ReactionRule& rule);
    boost::container::static_vector<ParticleID, 3>
    attempt_1to2_reaction(const PairSphericalDomain& dom,
            const ParticleID&, const Particle&, const ReactionRule& rule);

    std::array<ParticleID, 2> propagate(const PairSphericalDomain& dom, const Real dt);

    ReactionRule const& determine_single_reaction_rule(const Species&, const Real rnd);
    ReactionRule const& determine_pair_reaction_rule(const Species&, const Species&, const Real rnd);

    bool is_inside_of_shell(const PairSphericalDomain& dom,
                            const Real3& center, const Real radius);

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
#endif// ECELL4_NGFRD_PAIR_SPHERICAL_PROPAGATOR_HPP
