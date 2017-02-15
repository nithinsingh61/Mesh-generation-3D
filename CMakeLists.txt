# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( 3D_Code1_singvertices_ )

cmake_minimum_required(VERSION 2.8.10)

find_package(CGAL QUIET COMPONENTS Core )
add_compile_options(-std=c++11 -g)
if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  create_single_source_cgal_program( "bfs.cpp" )
  create_single_source_cgal_program( "filter.cpp" )
  create_single_source_cgal_program( "Improve_DT.cpp" )

else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()

