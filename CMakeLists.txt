cmake_minimum_required (VERSION 3.13.0)

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

add_library(smesh_io SHARED
  ${CMAKE_SOURCE_DIR}/src/smesh_io.f90
)

add_executable(smesh_run
  ${CMAKE_SOURCE_DIR}/src/smesh_run.f90
)
target_link_libraries(smesh_run PRIVATE smesh smesh_io)


# Set appropriate compile flags
# Note: Require Fortran 2018 standard due to use of C interoperability functionality
target_compile_options(smesh PUBLIC "-fPIC")
target_compile_options(smesh_io PUBLIC "-fPIC")
target_compile_options(smesh PRIVATE -Wall -Wextra -Werror -std=f2018)
target_compile_options(smesh_io PRIVATE -Wall -Wextra -Werror -std=f2018)
target_compile_options(smesh_run PRIVATE -Wall -Wextra -Werror -std=f2018)

install(TARGETS smesh smesh_io smesh_run)

# Handle code coverage generation
option(WITH_COVERAGE "Generate code coverage" OFF)
if (WITH_COVERAGE)
    target_compile_options(smesh     PRIVATE -g --coverage -O0)
    target_compile_options(smesh_io  PRIVATE -g --coverage -O0)
    target_compile_options(smesh_run PRIVATE -g --coverage -O0)
    target_link_options(smesh     PRIVATE --coverage)
    target_link_options(smesh_io  PRIVATE --coverage)
    target_link_options(smesh_run PRIVATE --coverage)

    add_custom_target(reset-coverage
      COMMAND lcov --directory . --zerocounters
    )
    add_custom_target(process-coverage
      COMMAND lcov --directory . --capture --output-file lcov.info
    )
endif()
