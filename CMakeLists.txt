cmake_minimum_required(VERSION 3.17)
project(bose_polaron_gpe VERSION 0.3.1 LANGUAGES CXX)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(ZLIB REQUIRED)
find_package(LibArchive REQUIRED)
add_library(gpe_solver src/solver.cpp src/continuation.cpp src/driven.cpp)
add_library(gpe::solver ALIAS gpe_solver)
set_target_properties(gpe_solver PROPERTIES EXPORT_NAME solver CXX_EXTENSIONS OFF)
target_compile_features(gpe_solver PUBLIC cxx_std_17)
target_include_directories(gpe_solver PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(gpe_solver PUBLIC Eigen3::Eigen PRIVATE ZLIB::ZLIB LibArchive::LibArchive)
add_executable(gpe src/main.cpp)
target_compile_features(gpe PRIVATE cxx_std_17)
target_link_libraries(gpe PRIVATE gpe::solver)
set_target_properties(gpe PROPERTIES CXX_EXTENSIONS OFF)
include(GNUInstallDirs)
install(TARGETS gpe gpe_solver EXPORT gpeTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT gpeTargets NAMESPACE gpe:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gpe)
include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/gpeConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/gpeConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gpe)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/gpeConfigVersion.cmake
    VERSION ${PROJECT_VERSION} COMPATIBILITY SameMajorVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/gpeConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/gpeConfigVersion.cmake DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/gpe)
include(CTest)
if(BUILD_TESTING)
    add_executable(test_numerics tests/numerics.cpp)
    target_compile_features(test_numerics PRIVATE cxx_std_17)
    target_link_libraries(test_numerics PRIVATE Eigen3::Eigen)
    set_target_properties(test_numerics PROPERTIES CXX_EXTENSIONS OFF)
    add_test(NAME numerics COMMAND test_numerics)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        add_test(NAME legacy_picard COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/legacy_picard.py ${CMAKE_CXX_COMPILER}
            "$<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>")
    endif()
    add_executable(test_picard tests/picard.cpp)
    target_link_libraries(test_picard PRIVATE gpe::solver)
    add_test(NAME picard COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/picard.py $<TARGET_FILE:test_picard>)
    add_executable(test_driven tests/driven.cpp)
    target_link_libraries(test_driven PRIVATE gpe::solver)
    add_test(NAME driven COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/driven.py $<TARGET_FILE:test_driven>)
    add_executable(test_library tests/library.cpp)
    target_link_libraries(test_library PRIVATE gpe::solver)
    add_test(NAME library COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/library.py $<TARGET_FILE:test_library> $<TARGET_FILE:gpe>)
    add_test(NAME integration COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration.py $<TARGET_FILE:gpe>)
    add_executable(test_continuation tests/continuation.cpp)
    target_link_libraries(test_continuation PRIVATE gpe::solver)
    add_test(NAME continuation COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/continuation.py $<TARGET_FILE:test_continuation>)
endif()
