#include "MantidKernel/FilteredTimeSeriesProperty.h"

#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <boost/python/return_value_policy.hpp>

using Mantid::Kernel::TimeSeriesProperty;
using Mantid::Kernel::FilteredTimeSeriesProperty;
using namespace boost::python;

namespace
{
  /// Macro to reduce copy-and-paste
  #define EXPORT_FILTEREDTIMESERIES_PROP(TYPE, Prefix)\
    register_ptr_to_python<FilteredTimeSeriesProperty<TYPE>*>();\
    register_ptr_to_python<const FilteredTimeSeriesProperty<TYPE>*>();\
    implicitly_convertible<FilteredTimeSeriesProperty<TYPE>*,const FilteredTimeSeriesProperty<TYPE>*>();\
    \
    class_<FilteredTimeSeriesProperty<TYPE>, bases<TimeSeriesProperty<TYPE> >, boost::noncopyable>(#Prefix"FilteredTimeSeriesProperty", no_init)\
      .def("unfiltered", &FilteredTimeSeriesProperty<TYPE>::unfiltered, return_value_policy<return_by_value>(),\
           "Returns a time series containing the unfiltered data") \
      ;
}

void export_FilteredTimeSeriesProperty()
{
  EXPORT_FILTEREDTIMESERIES_PROP(double, Float);
  EXPORT_FILTEREDTIMESERIES_PROP(bool, Bool);
  EXPORT_FILTEREDTIMESERIES_PROP(int32_t, Int32);
  EXPORT_FILTEREDTIMESERIES_PROP(int64_t, Int64);
  EXPORT_FILTEREDTIMESERIES_PROP(std::string, String);
}

