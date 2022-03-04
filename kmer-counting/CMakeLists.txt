cmake_minimum_required(VERSION 3.1)
project(Polyassembly)

set (CMAKE_CXX_STANDARD 17)
set (-D_GLIBCXX_USE_CXX11_ABI=0)
#add for debugging information
#set(CMAKE_BUILD_TYPE RelWithDebInfo)

enable_testing()

find_package(PkgConfig REQUIRED)
pkg_check_modules(JELLYFISH REQUIRED jellyfish-2.0)

link_directories(${JELLYFISH_LIBRARY_DIRS})

include_directories(${JELLYFISH_INCLUDE_DIRS} ${CMAKE_BINARY_DIR}/src)
add_compile_options(${JELLYFISH_CFLAGS_OTHER})

add_subdirectory(src)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
	message(STATUS "${dir}")
endforeach()
