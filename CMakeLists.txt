cmake_minimum_required(VERSION 3.12...3.26)

project(sandbox LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

get_property(is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT is_multi_config AND NOT (CMAKE_BUILD_TYPE OR DEFINED ENV{CMAKE_BUILD_TYPE}))
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Release default")
endif()

find_package(OpenMP COMPONENTS CXX)

if (MSVC)
    # warning level 4 (and all warnings as errors, /WX)
    add_compile_options(/W4)
else()
    # lots of warnings and all warnings as errors
    add_compile_options(-Wall -Wextra -pedantic)
endif()

include(FetchContent)
FetchContent_Declare(
  ultimaille
  GIT_REPOSITORY https://github.com/ssloy/ultimaille
  GIT_TAG        v2.0.0
)
FetchContent_MakeAvailable(ultimaille)
include_directories(${ultimaille_SOURCE_DIR})

# Define repository information
set(GITHUB_REPO "BrunoLevy/GraphiteThree")
set(ZIP_FILE "${CMAKE_SOURCE_DIR}/graphite.zip")
set(EXTRACT_DIR "${CMAKE_BINARY_DIR}")

# Get the latest release data using curl
execute_process(
    COMMAND curl -s https://api.github.com/repos/${GITHUB_REPO}/releases/latest
    OUTPUT_VARIABLE RELEASE_DATA
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Check if the release data was retrieved successfully
if(NOT RELEASE_DATA)
    message(FATAL_ERROR "Failed to get the latest release data.")
endif()

# Debug output to check the fetched JSON
message(STATUS "Release data: ${RELEASE_DATA}")

# Initialize the variable for the latest release URL
set(LATEST_RELEASE_URL "")

# Extract the download URL for the Windows 64-bit zip file
string(REGEX MATCH "\"browser_download_url\": \"[^\"]*win64.zip\"" URL_MATCH "${RELEASE_DATA}")
if(URL_MATCH)
    string(REGEX REPLACE ".*\"browser_download_url\": \"([^\"]*)\".*" "\\1" LATEST_RELEASE_URL "${URL_MATCH}")
endif()

# Check if the URL was retrieved successfully
if(NOT LATEST_RELEASE_URL)
    message(FATAL_ERROR "Failed to get the latest release URL.")
endif()

# Debug output to check the extracted URL
message(STATUS "Latest release URL: ${LATEST_RELEASE_URL}")

# Download the latest release zip file
file(DOWNLOAD ${LATEST_RELEASE_URL} ${ZIP_FILE})

# Unzip the downloaded file
execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xzf ${ZIP_FILE}
    WORKING_DIRECTORY ${EXTRACT_DIR}
)

# Remove the downloaded zip file
file(REMOVE_RECURSE ${ZIP_FILE})

add_subdirectory(project)