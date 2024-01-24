cmake_minimum_required (VERSION 3.5.1)

# Get project version from `VERSION` file
file(READ "${CMAKE_SOURCE_DIR}/VERSION" version_file)
# Parse string
if (${version_file} MATCHES "(([0-9]+)\\.([0-9]+)\\.([0-9]+)([+-]([A-Za-z0-9_-]+))?)")
  set(SMESH_VERSION ${CMAKE_MATCH_1})
  set(SMESH_VERSION_MAJOR ${CMAKE_MATCH_2})
  set(SMESH_VERSION_MINOR ${CMAKE_MATCH_3})
  set(SMESH_VERSION_PATCH ${CMAKE_MATCH_4})
  set(SMESH_VERSION_PRERELEASE ${CMAKE_MATCH_6})
else()
  message(FATAL_ERROR "Unable to parse version from `VERSION` file")
endif()

project(smesh
        VERSION ${SMESH_VERSION_MAJOR}.${SMESH_VERSION_MINOR}.${SMESH_VERSION_PATCH}
        LANGUAGES Fortran)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

add_library(smesh SHARED
  ${CMAKE_SOURCE_DIR}/src/smesh.f90
)

add_library(input_output SHARED
  ${CMAKE_SOURCE_DIR}/src/input_output.f90
)

add_executable(test_smesh
  ${CMAKE_SOURCE_DIR}/src/test_smesh.f90
)
target_link_libraries(test_smesh PRIVATE smesh input_output)

install(TARGETS smesh input_output test_smesh)
