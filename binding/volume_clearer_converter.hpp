#ifndef BINDING_VOLUME_CLEARER_CONVERTER_HPP
#define BINDING_VOLUME_CLEARER_CONVERTER_HPP

#include <boost/python.hpp>
#include <boost/utility/in_place_factory.hpp>
#include <boost/scoped_ptr.hpp>
#include "peer/utils.hpp"

namespace binding {

template<typename Tbase_>
class VolumeClearerWrapper: public Tbase_
{
public:
    typedef Tbase_ wrapped_type;
    typedef typename wrapped_type::particle_shape_type particle_shape_type;
    typedef typename wrapped_type::particle_id_type particle_id_type;

public:
    virtual ~VolumeClearerWrapper() {}

    virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore)
    {
        return callable_(shape, ignore);
    }

    virtual bool operator()(particle_shape_type const& shape, particle_id_type const& ignore0, particle_id_type const& ignore1)
    {
        return callable_(shape, ignore0, ignore1);
    }

    VolumeClearerWrapper(boost::python::object callable): callable_(callable) {}

private:
    boost::python::object callable_;
};

template<typename Tbase_>
struct volume_clearer_converter
{
    typedef VolumeClearerWrapper<Tbase_> native_type;

    static PyTypeObject const* expected_pytype()
    {
        return &PyBaseObject_Type;
    }

    static void* convert(PyObject* pyo)
    {
        if (!pyo || pyo == Py_None)
        {
            return 0;
        }

        if (!PyCallable_Check(pyo))
        {
            boost::python::throw_error_already_set();
        }

        native_type* retval(
            new native_type(
                boost::python::object(boost::python::borrowed(pyo))));
        peer::util::install_instance_holder<boost::scoped_ptr<native_type> >(pyo, boost::in_place(retval));
        return retval;
    }
};

template<typename Timpl_>
void register_volume_clearer_converter()
{
    using namespace boost::python;
    typedef Timpl_ impl_type;

    peer::util::to_native_lvalue_converter<impl_type, volume_clearer_converter<impl_type> >();
}

} // namespace binding

#endif /* BINDING_VOLUME_CLEARER_CONVERTER_HPP */
