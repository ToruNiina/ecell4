#ifndef MATRIX_SPACE_HPP
#define MATRIX_SPACE_HPP

// #include <iostream>
#include <cstddef>
#include <algorithm>
#include <iterator>
#include <boost/multi_array.hpp>
#include <boost/mpl/if.hpp>
#include <boost/range/size.hpp>
#include <boost/range/difference_type.hpp>
// #include "Vector3.hpp"
#include "Real3Type.hpp"
#include "sorted_list.hpp"
#include "utils/array_helper.hpp"
#include "utils/range.hpp"
#include "utils/unassignable_adapter.hpp"

#include <ecell4/core/Integer3.hpp>

namespace ecell4
{
namespace egfrd
{

template<typename Tobj_, typename Tkey_>
class MatrixSpace
{
public:
    typedef typename Tobj_::length_type length_type;
    typedef Tkey_ key_type;
    typedef Tobj_ mapped_type;
    // typedef Vector3<length_type> position_type;
    typedef ecell4::Real3 position_type;

    // typedef std::pair<const key_type, mapped_type> value_type;
    typedef std::pair<key_type, mapped_type> value_type;
    typedef std::vector<value_type> all_values_type;

    typedef sorted_list<std::vector<typename all_values_type::size_type> > cell_type;
    typedef boost::multi_array<cell_type, 3> matrix_type;
    typedef typename cell_type::size_type size_type;
    typedef std::array<typename matrix_type::size_type, 3>
            cell_index_type;
    typedef std::array<typename matrix_type::difference_type, 3>
            cell_offset_type;
    typedef std::unordered_map<key_type, typename all_values_type::size_type>
            key_to_value_mapper_type;

    typedef typename all_values_type::iterator iterator;
    typedef typename all_values_type::const_iterator const_iterator;
    typedef typename all_values_type::reference reference;
    typedef typename all_values_type::const_reference const_reference;

    typedef ecell4::Integer3 matrix_sizes_type;

private:
    typedef std::pair<key_type, mapped_type> nonconst_value_type;

public:

    MatrixSpace(
        const position_type& edge_lengths, const matrix_sizes_type& matrix_sizes)
        : edge_lengths_(edge_lengths),
          cell_sizes_(
            edge_lengths[0] / matrix_sizes[0],
            edge_lengths[1] / matrix_sizes[1],
            edge_lengths[2] / matrix_sizes[2]),
          matrix_(
            boost::extents[matrix_sizes[0]][matrix_sizes[1]][matrix_sizes[2]])
    {
        ;
    }

    // ~MatrixSpace()
    // {
    //     std::cerr << "MatrixSpace was released." << std::endl; //XXX: DEBUG
    // }

    inline cell_index_type index(const position_type& pos,
            double t = 1e-10) const
    {
        return array_gen<typename matrix_type::size_type>(
            static_cast<typename matrix_type::size_type>(
                pos[0] / cell_sizes_[0]) % matrix_.shape()[0],
            static_cast<typename matrix_type::size_type>(
                pos[1] / cell_sizes_[1]) % matrix_.shape()[1],
            static_cast<typename matrix_type::size_type>(
                pos[2] / cell_sizes_[2]) % matrix_.shape()[2]);
    }

    inline bool offset_index(
            cell_index_type& i,
            const cell_offset_type& o) const
    {
        if ((o[0] < 0 && static_cast<size_type>(-o[0]) > i[0])
                || (matrix_.shape()[0] - o[0] <= i[0])
                || (o[1] < 0 && static_cast<size_type>(-o[1]) > i[1])
                || (matrix_.shape()[1] - o[1] <= i[1])
                || (o[2] < 0 && static_cast<size_type>(-o[2]) > i[2])
                || (matrix_.shape()[2] - o[2] <= i[2]))
        {
            return false;
        }
        i[0] += o[0];
        i[1] += o[1];
        i[2] += o[2];
        return true;
    }

    inline position_type offset_index_cyclic(cell_index_type& i,
                                             const cell_offset_type& o) const
    {
        position_type retval;

        if (o[0] < 0 &&
            static_cast<typename matrix_type::size_type>(-o[0]) > i[0])
        {
            typename matrix_type::size_type t(
                (i[0] + matrix_.shape()[0] - (-o[0] % matrix_.shape()[0])) %
                matrix_.shape()[0]);
            retval[0] 
                = (o[0] - 
                   static_cast<typename matrix_type::difference_type>
                   (t - i[0])) * cell_sizes_[0];
            i[0] = t;
        }
        else if (matrix_.shape()[0] - o[0] <= i[0])
        {
            typename matrix_type::size_type t(
                    (i[0] + (o[0] % matrix_.shape()[0])) % matrix_.shape()[0]);
            retval[0] 
                = (o[0] - 
                   static_cast<typename matrix_type::difference_type>
                   (t - i[0])) * cell_sizes_[0];
            i[0] = t;
        }
        else
        {
            i[0] += o[0];
        }

        if (o[1] < 0 &&
                static_cast<typename matrix_type::size_type>(-o[1]) > i[1])
        {
            typename matrix_type::size_type t(
                    (i[1] + matrix_.shape()[1] - (-o[1] % matrix_.shape()[1])) %
                        matrix_.shape()[1]);
            retval[1] = (o[1] - static_cast<typename matrix_type::difference_type>(t - i[1])) * cell_sizes_[1];
            i[1] = t;
        }
        else if (matrix_.shape()[1] - o[1] <= i[1])
        {
            typename matrix_type::size_type t(
                    (i[1] + (o[1] % matrix_.shape()[1])) % matrix_.shape()[1]);
            retval[1] = (o[1] - static_cast<typename matrix_type::difference_type>(t - i[1])) * cell_sizes_[1];
            i[1] = t;
        }
        else
        {
            i[1] += o[1];
        }

        if (o[2] < 0 &&
                static_cast<typename matrix_type::size_type>(-o[2]) > i[2])
        {
            typename matrix_type::size_type t(
                    (i[2] + matrix_.shape()[2] - (-o[2] % matrix_.shape()[2])) %
                        matrix_.shape()[2]);
            retval[2] = (o[2] - static_cast<typename matrix_type::difference_type>(t - i[2])) * cell_sizes_[2];
            i[2] = t;
        }
        else if (matrix_.shape()[2] - o[2] <= i[2])
        {
            typename matrix_type::size_type t(
                    (i[2] + (o[2] % matrix_.shape()[2])) % matrix_.shape()[2]);
            retval[2] = (o[2] - static_cast<typename matrix_type::difference_type>(t - i[2])) * cell_sizes_[2];
            i[2] = t;
        }
        else
        {
            i[2] += o[2];
        }

        return retval;
    }

    inline const cell_type& cell(const cell_index_type& i) const
    {
        return matrix_[i[0]][i[1]][i[2]];
    }

    inline cell_type& cell(const cell_index_type& i)
    {
        return matrix_[i[0]][i[1]][i[2]];
    }

    inline const position_type& edge_lengths() const
    {
        return edge_lengths_;
    }

    inline const position_type& cell_sizes() const
    {
        return cell_sizes_;
    }

    inline const matrix_sizes_type matrix_sizes() const
    {
        typedef typename matrix_type::size_type matrix_size_type;
        const matrix_size_type* sizes(matrix_.shape());
        return matrix_sizes_type(sizes[0], sizes[1], sizes[2]);
    }

    inline size_type size() const
    {
        return values_.size();
    }

    inline iterator update(iterator const& old_value, const value_type& v)
    {
        cell_type* new_cell(&cell(index(v.second.position())));
        cell_type* old_cell(0);

        if (old_value != values_.end())
            old_cell = &cell(index((*old_value).second.position()));

        if (new_cell == old_cell)
        {
            reinterpret_cast<nonconst_value_type&>(*old_value) = v;
            return old_value;
        }
        else
        {
            typename all_values_type::size_type index(0);

            if (old_cell)
            {
                reinterpret_cast<nonconst_value_type&>(*old_value) = v;

                typename cell_type::iterator i(
                        old_cell->find(old_value - values_.begin()));
                index = *i;
                old_cell->erase(i);
                new_cell->push(index);
            }
            else
            {
                index = values_.size();
                values_.push_back(v);
                new_cell->push(index);
                rmap_[v.first] = index;
            }
            return values_.begin() + index;
        }
    }

    inline std::pair<iterator, bool> update(const value_type& v)
    {
        cell_type* new_cell(&cell(index(v.second.position())));
        typename all_values_type::iterator old_value(values_.end());
        cell_type* old_cell(0);

        {
            typename key_to_value_mapper_type::const_iterator i(rmap_.find(v.first));
            if (i != rmap_.end())
            {
                old_value = values_.begin() + (*i).second;
                old_cell = &cell(index(old_value->second.position()));
            }
        }

        if (new_cell == old_cell)
        {
            reinterpret_cast<nonconst_value_type&>(*old_value) = v;
            return std::pair<iterator, bool>(old_value, false);
        }
        else
        {
            typename all_values_type::size_type index(0);

            if (old_cell)
            {
                reinterpret_cast<nonconst_value_type&>(*old_value) = v;

                typename cell_type::iterator i(
                        old_cell->find(old_value - values_.begin()));
                index = *i;
                old_cell->erase(i);
                new_cell->push(index);
                return std::pair<iterator, bool>(values_.begin() + index, false);
            }
            else
            {
                index = values_.size();
                values_.push_back(v);
                new_cell->push(index);
                rmap_[v.first] = index;
                return std::pair<iterator, bool>(values_.begin() + index, true);
            }
        }
    }

    inline bool erase(iterator const& i)
    {
        if (end() == i)
        {
            return false;
        }

        typename all_values_type::size_type const old_index(i - values_.begin());

        bool is_succeeded(cell(index((*i).second.position())).erase(old_index));
        BOOST_ASSERT(is_succeeded);
        // BOOST_ASSERT(cell(index((*i).second.position())).erase(old_index));
        rmap_.erase((*i).first);

        typename all_values_type::size_type const last_index(values_.size() - 1);

        if (old_index < last_index)
        {
            value_type const& last(values_[last_index]);
            cell_type& old_c(cell(index(last.second.position())));
            is_succeeded = old_c.erase(last_index);
            BOOST_ASSERT(is_succeeded);
            // BOOST_ASSERT(old_c.erase(last_index));
            old_c.push(old_index);
            rmap_[last.first] = old_index;
            reinterpret_cast<nonconst_value_type&>(*i) = last; 
        }
        values_.pop_back();
        return true;
    }

    inline bool erase(const key_type& k)
    {
        typename key_to_value_mapper_type::const_iterator p(rmap_.find(k));
        if (rmap_.end() == p)
        {
            return false;
        }
        return erase(values_.begin() + (*p).second);
    }

    inline void clear()
    {
        for (typename matrix_type::element *p(matrix_.data()),
                                           *e(matrix_.data()
                                              + matrix_.num_elements());
             p != e; ++p)
        {
            (*p).clear();
        }
        rmap_.clear();
    }

    inline iterator begin()
    {
        return values_.begin();
    }

    inline const_iterator begin() const
    {
        return values_.begin();
    }

    inline iterator end()
    {
        return values_.end();
    }

    inline const_iterator end() const
    {
        return values_.end();
    }

    inline iterator find(const key_type& k)
    {
        typename key_to_value_mapper_type::const_iterator p(rmap_.find(k));
        if (rmap_.end() == p)
        {
            return values_.end();
        }
        return values_.begin() + (*p).second;
    }

    inline const_iterator find(const key_type& k) const
    {
        typename key_to_value_mapper_type::const_iterator p(rmap_.find(k));
        if (rmap_.end() == p)
        {
            return values_.end();
        }
        return values_.begin() + (*p).second;
    }

    template<typename Tcollect_>
    inline void each_neighbor(const cell_index_type& idx, Tcollect_& collector)
    {
        if(size() == 0)
        {
            return;
        }
        each_neighbor_loops<Tcollect_>(idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor(const cell_index_type& idx, Tcollect_ const& collector)
    {
        if(size() == 0)
        {
            return;
        }
        each_neighbor_loops<Tcollect_ const>(idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor(const cell_index_type& idx, Tcollect_& collector) const
    {
        if(size() == 0)
        {
            return;
        }
        each_neighbor_loops<Tcollect_>(idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor(const cell_index_type& idx, Tcollect_ const& collector) const
    {
        if(size() == 0)
        {
            return;
        }
        each_neighbor_loops<Tcollect_ const>(idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor_cyclic(const cell_index_type& idx,
            Tcollect_& collector)
    {
        if(size() == 0)
        {
            return;
        }
        each_neighbor_cyclic_loops<Tcollect_>(idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor_cyclic(const cell_index_type& idx,
            Tcollect_ const& collector)
    {
        if(size() == 0)
        {
            return;
        }
        each_neighbor_cyclic_loops<Tcollect_ const>(idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor_cyclic(const cell_index_type& idx,
            Tcollect_& collector) const
    {
        if(size() == 0)
        {
            return;
        }
        each_neighbor_cyclic_loops<Tcollect_>(idx, collector);
    }

    template<typename Tcollect_>
    inline void each_neighbor_cyclic(const cell_index_type& idx,
            Tcollect_ const& collector) const
    {
        if(size() == 0)
        {
            return;
        }
        each_neighbor_cyclic_loops<Tcollect_ const>(idx, collector);
    }

private:
    std::pair<cell_type*, cell_type*> cell_range()
    {
        return std::make_pair(
            matrix_.origin(), matrix_.origin() + matrix_.num_elements());
    }

    std::pair<cell_type const*, cell_type const*> cell_range() const
    {
        return std::make_pair(
            matrix_.origin(), matrix_.origin() + matrix_.num_elements());
    }

    template<typename Tcollect_>
    inline void each_neighbor_loops(const cell_index_type& idx,
                                    Tcollect_& collector) const
    {
        cell_offset_type off;

        for (off[2] = -1; off[2] <= 1; ++off[2])
        {
            for (off[1] = -1; off[1] <= 1; ++off[1])
            {
                for (off[0] = -1; off[0] <= 1; ++off[0])
                {
                    cell_index_type _idx(idx);
                    if (!offset_index(_idx, off)) {
                        continue;
                    }
                    cell_type const& c(cell(_idx));
                    for (typename cell_type::const_iterator i(c.begin()); i != c.end(); ++i) 
                    {
                        collector(values_.begin() + *i, position_type());
                    }
                }
            }
        }
    }

    template<typename Tcollect_>
    inline void each_neighbor_loops(const cell_index_type& idx,
                                    Tcollect_& collector)
    {
        cell_offset_type off;

        for (off[2] = -1; off[2] <= 1; ++off[2])
        {
            for (off[1] = -1; off[1] <= 1; ++off[1])
            {
                for (off[0] = -1; off[0] <= 1; ++off[0])
                {
                    cell_index_type _idx(idx);
                    if (!offset_index(_idx, off)) {
                        continue;
                    }
                    cell_type const& c(cell(_idx));
                    for (typename cell_type::const_iterator i(c.begin()); i != c.end(); ++i) 
                    {
                        collector(values_.begin() + *i, position_type());
                    }
                }
            }
        }
    }

    template<typename Tcollect_>
    inline void each_neighbor_cyclic_loops(const cell_index_type& idx,
                                           Tcollect_& collector) const
    {
        cell_offset_type off;

        for (off[2] = -1; off[2] <= 1; ++off[2])
        {
            for (off[1] = -1; off[1] <= 1; ++off[1])
            {
                for (off[0] = -1; off[0] <= 1; ++off[0])
                {
                    cell_index_type _idx(idx);
                    const position_type pos_off(offset_index_cyclic(_idx, off));
                    cell_type const& c(cell(_idx));
                    for (typename cell_type::const_iterator i(c.begin()); i != c.end(); ++i) 
                    {
                        collector(values_.begin() + *i, pos_off);
                    }
                }
            }
        }
    }

    template<typename Tcollect_>
    inline void each_neighbor_cyclic_loops(const cell_index_type& idx,
                                           Tcollect_& collector)
    {
        cell_offset_type off;

        for (off[2] = -1; off[2] <= 1; ++off[2])
        {
            for (off[1] = -1; off[1] <= 1; ++off[1])
            {
                for (off[0] = -1; off[0] <= 1; ++off[0])
                {
                    cell_index_type _idx(idx);
                    const position_type pos_off(offset_index_cyclic(_idx, off));
                    cell_type const& c(cell(_idx));
                    for (typename cell_type::const_iterator i(c.begin()); i != c.end(); ++i) 
                    {
                        collector(values_.begin() + *i, pos_off);
                    }
                }
            }
        }
    }

private:
    const position_type edge_lengths_;
    const position_type cell_sizes_;
    matrix_type matrix_;
    key_to_value_mapper_type rmap_;
    all_values_type values_;
};

template<typename T_, typename Tkey_>
static inline typename MatrixSpace<T_, Tkey_>::cell_index_type&
operator+=(
       typename MatrixSpace<T_,
                Tkey_>::cell_index_type& lhs,
       const typename MatrixSpace<T_,
                Tkey_>::cell_offset_type& rhs)
{
    rhs[0] += lhs[0];
    rhs[1] += lhs[1];
    rhs[2] += lhs[2];
    return rhs;
}

template<typename T_, typename Tkey_>
struct is_sized<MatrixSpace<T_, Tkey_> >: boost::mpl::true_ {};

template<typename T_, typename Tkey_>
struct range_size<MatrixSpace<T_, Tkey_> >
{
    typedef typename MatrixSpace<T_, Tkey_>::size_type type;
};

template<typename T_, typename Tkey_>
struct range_size_retriever<MatrixSpace<T_, Tkey_> >
{
    typedef MatrixSpace<T_, Tkey_> argument_type;
    typedef typename range_size<argument_type>::type result_type;

    result_type operator()(argument_type const& range) const
    {
        return range.size();
    }
};

} // egfrd
} // ecell4
#endif /* MATRIX_SPACE_HPP */
