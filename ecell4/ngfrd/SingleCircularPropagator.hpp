#ifndef ECELL4_NGFRD_SINGLE_CIRCULAR_PROPAGATOR_HPP
#define ECELL4_NGFRD_SINGLE_CIRCULAR_PROPAGATOR_HPP
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/Model.hpp>
#include <ecell4/ngfrd/ReactionInfo.hpp>
#include <ecell4/ngfrd/SingleCircularDomain.hpp>
#include <ecell4/ngfrd/NGFRDWorld.hpp>

#include <boost/container/static_vector.hpp>

namespace ecell4
{
namespace ngfrd
{

class NGFRDSimulator; // forward decl

class SingleCircularPropagator
{
  public:
    SingleCircularPropagator(const DomainID& self,
        const Model& model, NGFRDWorld& world, NGFRDSimulator& sim,
        RandomNumberGenerator& rng, const std::size_t max_retry,
        std::vector<std::pair<ReactionRule, ReactionInfo>>& last_reactions)
        : model_(model), world_(world), sim_(sim), rng_(rng), self_id_(self),
          max_retry_count_(max_retry), last_reactions_(last_reactions)
    {}

    boost::container::static_vector<ParticleID, 2>
    operator()(const SingleCircularDomain& dom)
    {
        switch(dom.eventkind())
        {
            case SingleCircularDomain::EventKind::Escape:
            {
                return this->escape(dom);
            }
            case SingleCircularDomain::EventKind::Reaction:
            {
                return this->reaction(dom);
            }
            default:
            {
                throw_exception<IllegalState>("ecell4::ngfrd::SingleCircular"
                    "Propagator: invalid eventkind");
            }
        }
    }

    ParticleID burst(const SingleCircularDomain& dom, const Real time_now)
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

    boost::container::static_vector<ParticleID, 2>
    escape(const SingleCircularDomain& dom);

    boost::container::static_vector<ParticleID, 2>
    reaction(const SingleCircularDomain& dom);

    boost::container::static_vector<ParticleID, 2>
    attempt_1to1_reaction(const SingleCircularDomain& dom,
            const ParticleID&, const Particle&, const FaceID& fid,
            const ReactionRule& rule);
    boost::container::static_vector<ParticleID, 2>
    attempt_1to2_reaction(const SingleCircularDomain& dom,
            const ParticleID&, const Particle&, const FaceID& fid,
            const ReactionRule& rule);

    ParticleID propagate(const SingleCircularDomain& dom, const Real dt);

    ReactionRule const& determine_reaction_rule(const Species& sp, const Real rnd);

    bool is_inside_of_shell(const SingleCircularDomain& dom,
            const std::pair<Real3, FaceID>& pos, const Real radius);

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
#endif// ECELL4_NGFRD_SINGLE_CIRCULAR_PROPAGATOR_HPP