cmake_minimum_required(VERSION 3.26)

set(KALLES_PROJECT_NAME "kalles-fraktaler" CACHE STRING "Project name")
project("${KALLES_PROJECT_NAME}" LANGUAGES CXX)

add_subdirectory(fraktal_sft)
add_subdirectory(common)
