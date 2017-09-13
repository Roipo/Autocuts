cmake_minimum_required(VERSION 3.0)
SET(CMAKE_C_COMPILER /usr/local/Cellar/llvm/5.0.0/bin/clang)
SET(CMAKE_CXX_COMPILER /usr/local/Cellar/llvm/5.0.0/bin/clang++)

project(Autocuts)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(LIBIGL_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libigl/include)
set(LIBIGL_INCLUDE_DIRS ${LIBIGL_INCLUDE_DIR}  ${LIBIGL_INCLUDE_DIR}/../external/nanogui/include)

# Compilation flags: adapt to your needs 
if(MSVC)
  # Enable parallel compilation
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj") 
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
  # Libigl requires a modern C++ compiler that supports c++11
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "." )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")
endif()



# libigl options: choose between header only and compiled static library
# Header-only is preferred for small projects. For larger projects the static build
# considerably reduces the compilation times
option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" OFF)

# add a customizable menu bar
option(LIBIGL_WITH_NANOGUI     "Use Nanogui menu"   ON)
# libigl options: choose your dependencies (by default everything is OFF except opengl) 
option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_OPENGL      "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW"           ON)
option(LIBIGL_WITH_BBW         "Use BBW"            OFF)
option(LIBIGL_WITH_EMBREE      "Use Embree"         OFF)
option(LIBIGL_WITH_PNG         "Use PNG"            ON)
option(LIBIGL_WITH_TETGEN      "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE    "Use Triangle"       OFF)
option(LIBIGL_WITH_XML         "Use XML"            OFF)
option(LIBIGL_WITH_LIM         "Use LIM"            OFF)
option(LIBIGL_WITH_COMISO      "Use CoMiso"         OFF)
option(LIBIGL_WITH_MATLAB      "Use Matlab"         OFF) # This option is not supported yet
option(LIBIGL_WITH_MOSEK       "Use MOSEK"          OFF) # This option is not supported yet
option(LIBIGL_WITH_CGAL        "Use CGAL"           OFF)
option(LIBIGL_WITH_PYTHON        "Use Python"           OFF)

## definition for schuller's fork
option(LIBIGL_VIEWER_WITH_NANOGUI               "Use Nanogui menu"             ON)
option(LIBIGL_VIEWER_WITH_NANOGUI_MULTIMESH     "Show mesh selection"          ON)
option(LIBIGL_VIEWER_WITH_NANOGUI_IO            "Show IO menu"                 ON)
option(LIBIGL_VIEWER_WITH_NANOGUI_SERIALIZATION "Show serialization menu"      ON)
# Adding libigl: choose the path to your local copy libigl
# This is going to compile everything you requested 
#message(FATAL_ERROR "${PROJECT_SOURCE_DIR}/../libigl/cmake")
add_subdirectory("${LIBIGL_INCLUDE_DIR}/../shared/cmake" "libigl")
# libigl information
message("libigl includes: ${LIBIGL_INCLUDE_DIRS}")
message("libigl libraries: ${LIBIGL_LIBRARIES}")
message("libigl extra sources: ${LIBIGL_EXTRA_SOURCES}")
message("libigl extra libraries: ${LIBIGL_EXTRA_LIBRARIES}")
message("libigl definitions: ${LIBIGL_DEFINITIONS}")

OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
ENDIF()

set(SRCFILES
   src/ApplicationPipeline.cpp
   src/Newton.cpp
   src/SolverPlugin.cpp
   src/BBox.cpp
   src/SolverWrapper.cpp
   src/DataObject.cpp
   src/Position.cpp
   src/Utils.cpp
   src/Energy.cpp
   src/Separation.cpp
   src/main.cpp
   src/EnergySymDir.cpp
   src/Solver.cpp)

OPTION (USE_PARDISO "Use Pardiso" OFF)
IF(USE_PARDISO)
    message("Use PARDISO")
    set(LIBIGL_EXTRA_LIBRARIES ${LIBIGL_EXTRA_LIBRARIES} ${CMAKE_SOURCE_DIR}/libpardiso500-MACOS-X86-64.dylib)
    set(SRCFILES ${SRCFILES} PardisoSolver.cpp)
ENDIF()
# Prepare the build environment
include_directories(${LIBIGL_INCLUDE_DIRS}
    ${LIBIGL_INCLUDE_DIR}/../external/nanogui/include/
    ${LIBIGL_INCLUDE_DIR}/../external/nanogui/ext/nanovg/src)
add_definitions(${LIBIGL_DEFINITIONS} -DENABLE_SERIALIZATION_CORE_ONLY -D_CRT_SECURE_NO_WARNINGS -DNOMINMAX)


add_executable(${PROJECT_NAME}_bin ${SRCFILES} ${LIBIGL_EXTRA_SOURCES})
target_link_libraries(${PROJECT_NAME}_bin ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})
