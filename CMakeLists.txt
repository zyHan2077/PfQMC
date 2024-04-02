cmake_minimum_required(VERSION 3.13)

PROJECT(MQMC)

include_directories(
${PROJECT_SOURCE_DIR}/inc /opt/Eigen/eigen-3.4.0
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(SRC_LIST src/skewMatUtils.cpp src/operator.cpp)

set(MKL_INTERFACE lp64)

enable_testing()
find_package(MKL CONFIG REQUIRED)
find_package(GTest CONFIG REQUIRED)
#message(STATUS "${MKL_IMPORTED_TARGETS}") #Provides available list of targets based on input

# add_library( libyaml STATIC IMPORTED )

# set_target_properties( libyaml PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/lib/libyaml-cpp.a )

set(CMAKE_CXX_COMPILER "icpx")
set(CMAKE_C_LINK)
set(FP_MODEL_FLAGS "-fp-model=precise -fimf-arch-consistency=true -no-fma")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${FP_MODEL_FLAGS}")

# add_executable(main ./main.cpp ${SRC_LIST})
# target_compile_options(main PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
# target_include_directories(main PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
# target_link_libraries(main PUBLIC $<LINK_ONLY:MKL::MKL>)

# target_link_libraries(main_test GTest::gtest_main)
add_executable(main_test ./main_test.cpp ${SRC_LIST})
target_compile_options(main_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(main_test PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(main_test PUBLIC $<LINK_ONLY:MKL::MKL> GTest::gtest_main)

include(GoogleTest)
gtest_discover_tests(main_test)

# add_test(NAME mytest COMMAND myapp)
# if(MKL_ENV)
#   set_tests_properties(mytest PROPERTIES ENVIRONMENT "${MKL_ENV}")
# endif()
# add_executable(rapidtest src/rapidjsontest.cpp)
# target_link_libraries(main libyaml)