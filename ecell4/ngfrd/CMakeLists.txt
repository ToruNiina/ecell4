# file(GLOB CPP_FILES *.cpp)
set(CPP_FILES
    BDMath.cpp
    BDPropagator.cpp
    MultiDomain.cpp
    NGFRDSimulator.cpp
    SingleSphericalPropagator.cpp
    SingleCircularPropagator.cpp
    SingleConicalPropagator.cpp
    PairSphericalPropagator.cpp
    PairCircularPropagator.cpp
)

add_definitions(-DECELL4_NGFRD_LOG_DEBUG)
add_library(ecell4-ngfrd STATIC ${CPP_FILES})
target_link_libraries(ecell4-ngfrd INTERFACE ecell4-core greens_functions ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})

# add_subdirectory(tests)
add_subdirectory(samples)
