cmake_minimum_required(VERSION 3.13)

PROJECT(MQMC)

include_directories(
${PROJECT_SOURCE_DIR}/include
)

find_package(MKL CONFIG REQUIRED)
#message(STATUS "${MKL_IMPORTED_TARGETS}") #Provides available list of targets based on input

# add_library( libyaml STATIC IMPORTED )

# set_target_properties( libyaml PROPERTIES IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/lib/libyaml-cpp.a )

set(CMAKE_CXX_COMPILER "icpx")
set(CMAKE_C_LINK)

add_executable(main src/main.cpp)

target_compile_options(main PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(main PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(main PUBLIC $<LINK_ONLY:MKL::MKL>)

# add_test(NAME mytest COMMAND myapp)
# if(MKL_ENV)
#   set_tests_properties(mytest PROPERTIES ENVIRONMENT "${MKL_ENV}")
# endif()
# add_executable(rapidtest src/rapidjsontest.cpp)
# target_link_libraries(main libyaml)