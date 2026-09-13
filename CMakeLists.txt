cmake_minimum_required(VERSION 3.16...3.27)
project(dragon_g4_pilot CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(dragon_g4_pilot ${sources} ${headers})
target_include_directories(dragon_g4_pilot PRIVATE include)
target_link_libraries(dragon_g4_pilot PRIVATE ${Geant4_LIBRARIES})

foreach(_script vis.mac init_vis.mac run_chain.mac)
  configure_file(${PROJECT_SOURCE_DIR}/macros/${_script} ${PROJECT_BINARY_DIR}/${_script} COPYONLY)
endforeach()

# Reaction config files (see ReactionConfig.hh) -- copied next to the
# executable so --track-reaction/--vis Reaction find the default
# (o15ag_19ne.reaction) by a plain relative path, same convention as the
# macros above.
file(GLOB _reaction_files ${PROJECT_SOURCE_DIR}/reactions/*.reaction)
foreach(_reaction_file ${_reaction_files})
  get_filename_component(_reaction_name ${_reaction_file} NAME)
  configure_file(${_reaction_file} ${PROJECT_BINARY_DIR}/${_reaction_name} COPYONLY)
endforeach()
