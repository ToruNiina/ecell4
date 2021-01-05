#ifndef ECELL4_NGFRD_UTILITY_HPP
#define ECELL4_NGFRD_UTILITY_HPP
#include <algorithm>

namespace ecell4
{
namespace ngfrd
{

template<typename Container, typename Value>
bool unique_push_back(Container& cont, Value v)
{
    using std::begin;
    using std::end;
    if(std::find(begin(cont), end(cont), v) == end(cont))
    {
        cont.push_back(v);
        return true;
    }
    return false;
}

} // ngfrd
} // ecell4
#endif// ECELL4_NGFRD_UTILITY_HPP