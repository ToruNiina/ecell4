#ifndef ARRAY_HELPER_HPP
#define ARRAY_HELPER_HPP

#include <boost/array.hpp>
#include <boost/preprocessor/config/limits.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/call_traits.hpp>

namespace ecell4
{
namespace egfrd
{

#define ARRAY_HELPER_INNER_TPL(__z__, __n__, __d__) \
    __d__[__n__] = BOOST_PP_CAT(p, __n__);

#define ARRAY_HELPER_TPL(__z__, __n__, __d__) \
template<typename T_> \
inline ::std::array<T_, __n__> array_gen(\
        BOOST_PP_ENUM_PARAMS(__n__, T_ const& p)) \
{ \
    ::std::array<T_, __n__> retval; \
    BOOST_PP_REPEAT_ ## __z__(__n__, ARRAY_HELPER_INNER_TPL, retval) \
    return retval; \
}

#ifndef WIN32_MSC
BOOST_PP_REPEAT_FROM_TO(0, BOOST_PP_LIMIT_REPEAT, ARRAY_HELPER_TPL, )
#else
BOOST_PP_REPEAT_FROM_TO(0, 128, ARRAY_HELPER_TPL, )
#endif

#undef ARRAY_HELPER_TPL
#undef ARRAY_HELPER_INNER_TPL


} // egfrd
} // ecell4
#endif /* ARRAY_HELPER_HPP */
