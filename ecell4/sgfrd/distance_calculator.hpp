#ifndef ECELL4_SGFRD_DISTANCE_CALCULATOR
#define ECELL4_SGFRD_DISTANCE_CALCULATOR
#include <ecell4/core/Polygon.hpp>
#include <ecell4/core/Circle.hpp>
#include <ecell4/core/Cone.hpp>
#include <boost/variant.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include "distance.hpp"
#include "Shell.hpp"

namespace ecell4
{
namespace sgfrd
{

struct distance_calculator : public boost::static_visitor<Real>
{
    distance_calculator(const Real3& pos): pos_(pos){}

    template<typename shapeT, typename stridT>
    Real operator()(const Shell<shapeT, stridT>& sh) const
    {
        return ecell4::sgfrd::distance(pos_, sh.shape()) - sh.size();
    }

  private:
    Real3 pos_;
};

template<typename strID>
struct distance_calculator_on_surface : public boost::static_visitor<Real>
{
public:
    typedef strID structure_id_type;
    typedef ecell4::Polygon polygon_type;
    typedef typename polygon_type::FaceID   FaceID;
    typedef typename polygon_type::VertexID VertexID;

private:

    template<typename shapeT, typename stridT>
    struct is_2d_shell : boost::mpl::or_<
        boost::mpl::and_<
            boost::is_same<shapeT, ecell4::Circle>,
            boost::is_same<stridT, FaceID>
            >,
        boost::mpl::and_<
            boost::is_same<shapeT, ecell4::ConicalSurface>,
            boost::is_same<stridT, VertexID>
            >
        >
    {};

public:

    distance_calculator_on_surface(
            const std::pair<Real3, structure_id_type>& pos,
            const polygon_type& poly)
        : poly_(poly), pos_(pos)
    {}

    template<typename shapeT, typename stridT>
    typename boost::enable_if<is_2d_shell<shapeT, stridT>, Real>::type
    operator()(const Shell<shapeT, stridT>& sh) const
    {
        return ecell4::polygon::distance(poly_, pos_, sh.get_surface_position()) -
               sh.size();
    }

private:

    polygon_type const& poly_;
    std::pair<Real3, structure_id_type> pos_;
};

} // sgfrd
} // ecell4
#endif// ECELL4_SGFRD_DISTANCE_CALCULATOR
