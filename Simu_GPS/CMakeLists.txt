# (1)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(attenuationInBlock)

# (2)
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# (3)
include(${Geant4_USE_FILE})

find_package(ROOT REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIR})

# (4)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# (5)
add_executable(attenuationInBlock attenuationInBlock.cc ${sources} ${headers})
target_link_libraries(attenuationInBlock ${Geant4_LIBRARIES})
target_link_libraries(attenuationInBlock ${ROOT_LIBRARIES})

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# (6)
install(TARGETS attenuationInBlock DESTINATION bin)

# (7)
configure_file(vis.mac vis.mac COPYONLY)
