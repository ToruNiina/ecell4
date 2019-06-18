file(GLOB CPP_FILES *.cpp)

add_library(ecell4-egfrd STATIC ${CPP_FILES})
target_link_libraries(ecell4-egfrd INTERFACE ecell4-core)
target_link_libraries(ecell4-egfrd PRIVATE GSL::gsl GSL::gslcblas greens_functions)

# add_subdirectory(tests)
add_subdirectory(samples)
