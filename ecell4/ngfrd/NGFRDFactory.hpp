#ifndef ECELL4_NGFRD_NGFRD_FACTORY_HPP
#define ECELL4_NGFRD_NGFRD_FACTORY_HPP

#include <ecell4/core/SimulatorFactory.hpp>
#include <ecell4/core/RandomNumberGenerator.hpp>
#include <ecell4/core/STLFileIO.hpp>

#include "NGFRDWorld.hpp"
#include "NGFRDSimulator.hpp"


namespace ecell4
{

namespace ngfrd
{

class NGFRDFactory:
    public SimulatorFactory<NGFRDWorld, NGFRDSimulator>
{
public:

    typedef SimulatorFactory<NGFRDWorld, NGFRDSimulator> base_type;
    typedef base_type::world_type world_type;
    typedef base_type::simulator_type simulator_type;
    typedef NGFRDFactory this_type;

public:

    NGFRDFactory(const Integer3& matrix_sizes = default_matrix_sizes(),
                 Real bd_dt_factor_3D = default_bd_dt_factor_3D(),
                 Real bd_dt_factor_2D = default_bd_dt_factor_2D(),
                 Real bd_reaction_length_factor = default_bd_reaction_length_factor(),
                 Real tree_margin = default_tree_margin(),
                 Integer max_retry_count = default_max_retry_count())
        : base_type(), rng_(), matrix_sizes_(matrix_sizes),
          max_retry_count_(max_retry_count),
          bd_dt_factor_3D_(bd_dt_factor_3D),
          bd_dt_factor_2D_(bd_dt_factor_2D),
          bd_reaction_length_factor_(bd_reaction_length_factor),
          tree_margin_(tree_margin),
          polygon_(nullptr), polygon_file_("", STLFormat::Ascii)
    {
        ; // do nothing
    }

    virtual ~NGFRDFactory()
    {
        ; // do nothing
    }

    static inline Integer3 default_matrix_sizes()
    {
        return Integer3(3, 3, 3);
    }
    static inline Integer default_max_retry_count()
    {
        return 4; // ?
    }

    static inline Real default_bd_dt_factor_3D()
    {
        return 1e-5;
    }
    static inline Real default_bd_dt_factor_2D()
    {
        return 1e-3;
    }
    static inline Real default_bd_reaction_length_factor()
    {
        return 0.1; // 0.05 ~ 0.1
    }
    static inline Real default_tree_margin()
    {
        return 0.1; // tuned later
    }

    this_type& rng(const std::shared_ptr<RandomNumberGenerator>& rng)
    {
        rng_ = rng;
        return (*this);
    }

    inline this_type* rng_ptr(const std::shared_ptr<RandomNumberGenerator>& rng)
    {
        return &(this->rng(rng));  //XXX: == this
    }

    // -----------------------------------------------------------------------
    // get polygon from a list of triangles

    this_type& polygon(const Real3& el, const std::vector<Triangle>& ts)
    {
        polygon_file_.first = ""; // XXX clear polygon file

        this->polygon_ = std::make_shared<Polygon>(el, ts);
        return (*this);
    }
    this_type* polygon_ptr(const Real3& el, const std::vector<Triangle>& ts)
    {
        return std::addressof(this->polygon(el, ts));
    }

    // -----------------------------------------------------------------------
    // read polygon from .STL file

    this_type& polygon(const std::string& fname, const STLFormat fmt)
    {
        polygon_ = nullptr; // XXX clear polygon structure

        this->polygon_file_ = std::make_pair(fname, fmt);
        return (*this);
    }
    this_type* polygon_ptr(const std::string& fname, const STLFormat fmt)
    {
        return std::addressof(this->polygon(fname, fmt));
    }

protected:

    virtual world_type* create_world(const Real3& edge_lengths) const override
    {
        std::shared_ptr<RandomNumberGenerator> rng(nullptr);
        if (rng_)
        {
            rng = rng_;
        }
        else
        {
            rng = std::shared_ptr<RandomNumberGenerator>(new GSLRandomNumberGenerator());
            rng->seed();
        }

        if(this->polygon_)
        {
            return new world_type(edge_lengths, matrix_sizes_, tree_margin_,
                                  std::move(rng), this->polygon_);
        }
        else if(!this->polygon_file_.first.empty())
        {
            return new world_type(edge_lengths, matrix_sizes_, tree_margin_,
                    std::move(rng), std::make_shared<Polygon>(read_polygon(
                        polygon_file_.first, polygon_file_.second, edge_lengths)));
        }
        else
        {
            return new world_type(edge_lengths, matrix_sizes_, tree_margin_, std::move(rng));
        }
    }

    virtual simulator_type* create_simulator(
        const std::shared_ptr<world_type>& w, const std::shared_ptr<Model>& m) const override
    {
        return new simulator_type(w, m, bd_dt_factor_3D_, bd_dt_factor_2D_, bd_reaction_length_factor_, max_retry_count_);
    }


protected:

    std::shared_ptr<RandomNumberGenerator> rng_;
    std::shared_ptr<Polygon>               polygon_;
    std::pair<std::string, STLFormat>      polygon_file_;
    Integer3 matrix_sizes_;
    Integer  max_retry_count_;
    Real bd_dt_factor_3D_;
    Real bd_dt_factor_2D_;
    Real bd_reaction_length_factor_;
    Real tree_margin_;
};

} // ngfrd

} // ecell4

#endif /* ECELL4_NGFRD_NGFRD_FACTORY_HPP */
