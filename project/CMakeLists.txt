cmake_minimum_required(VERSION 3.12...3.26)

project(ultimaille-slim LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# include(libigl)

if(TARGET igl::core)
    return()
endif()

include(FetchContent)
FetchContent_Declare(
    libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(libigl)

FetchContent_Declare(
  param_parser
  GIT_REPOSITORY https://github.com/ultimaille/param-parser
  GIT_TAG        v1.0.2
)
FetchContent_MakeAvailable(param_parser)

include_directories(
    /usr/local/include/eigen3
    ${param_parser_SOURCE_DIR}
)
igl_include(glfw)

add_executable(SLIM_intern SLIM_intern.cpp)

target_link_libraries(SLIM_intern ${CMAKE_DL_LIBS} ultimaille $<$<BOOL:${OpenMP_CXX_FOUND}>:OpenMP::OpenMP_CXX>)
target_link_libraries(SLIM_intern igl::glfw)
target_link_libraries(SLIM_intern ${CMAKE_DL_LIBS} param-parser)

IF (NOT WIN32)
    target_link_libraries(SLIM_intern m)
ENDIF()

# Copy asset files to build directory
file(COPY mesh_test DESTINATION .)