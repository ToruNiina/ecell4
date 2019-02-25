#include "BDSimulator.hpp"

#include <boost/scoped_array.hpp>

#include <cstring>

namespace ecell4
{

namespace bd
{

void BDSimulator::step()
{
    last_reactions_.clear();

    {
        BDPropagator propagator(*model_, *world_, *rng(), dt(), last_reactions_);
        while (propagator())
        {
            ; // do nothing here
        }
    }

    set_t(t() + dt());
    num_steps_++;
}

bool BDSimulator::step(const Real& upto)
{
    const Real t0(t()), dt0(dt()), tnext(next_time());

    if (upto <= t0)
    {
        return false;
    }

    if (upto >= tnext)
    {
        step();
        return true;
    }
    else
    {
        dt_ = upto - t0;
        step();
        dt_ = dt0;
        return false;
    }
}

} // bd

} // ecell4
