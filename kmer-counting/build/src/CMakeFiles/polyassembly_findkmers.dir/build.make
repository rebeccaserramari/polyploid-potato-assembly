# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/rebecca/work/hifi-potato/paper-code/kmer-counting

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/rebecca/work/hifi-potato/paper-code/kmer-counting/build

# Include any dependencies generated for this target.
include src/CMakeFiles/polyassembly_findkmers.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/polyassembly_findkmers.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/polyassembly_findkmers.dir/flags.make

src/CMakeFiles/polyassembly_findkmers.dir/polyassembly_findkmers.cpp.o: src/CMakeFiles/polyassembly_findkmers.dir/flags.make
src/CMakeFiles/polyassembly_findkmers.dir/polyassembly_findkmers.cpp.o: ../src/polyassembly_findkmers.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/rebecca/work/hifi-potato/paper-code/kmer-counting/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/CMakeFiles/polyassembly_findkmers.dir/polyassembly_findkmers.cpp.o"
	cd /home/rebecca/work/hifi-potato/paper-code/kmer-counting/build/src && /home/rebecca/miniconda3/bin/x86_64-conda_cos6-linux-gnu-c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/polyassembly_findkmers.dir/polyassembly_findkmers.cpp.o -c /home/rebecca/work/hifi-potato/paper-code/kmer-counting/src/polyassembly_findkmers.cpp

src/CMakeFiles/polyassembly_findkmers.dir/polyassembly_findkmers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/polyassembly_findkmers.dir/polyassembly_findkmers.cpp.i"
	cd /home/rebecca/work/hifi-potato/paper-code/kmer-counting/build/src && /home/rebecca/miniconda3/bin/x86_64-conda_cos6-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/rebecca/work/hifi-potato/paper-code/kmer-counting/src/polyassembly_findkmers.cpp > CMakeFiles/polyassembly_findkmers.dir/polyassembly_findkmers.cpp.i

src/CMakeFiles/polyassembly_findkmers.dir/polyassembly_findkmers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/polyassembly_findkmers.dir/polyassembly_findkmers.cpp.s"
	cd /home/rebecca/work/hifi-potato/paper-code/kmer-counting/build/src && /home/rebecca/miniconda3/bin/x86_64-conda_cos6-linux-gnu-c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/rebecca/work/hifi-potato/paper-code/kmer-counting/src/polyassembly_findkmers.cpp -o CMakeFiles/polyassembly_findkmers.dir/polyassembly_findkmers.cpp.s

# Object files for target polyassembly_findkmers
polyassembly_findkmers_OBJECTS = \
"CMakeFiles/polyassembly_findkmers.dir/polyassembly_findkmers.cpp.o"

# External object files for target polyassembly_findkmers
polyassembly_findkmers_EXTERNAL_OBJECTS =

src/polyassembly_findkmers: src/CMakeFiles/polyassembly_findkmers.dir/polyassembly_findkmers.cpp.o
src/polyassembly_findkmers: src/CMakeFiles/polyassembly_findkmers.dir/build.make
src/polyassembly_findkmers: src/libPolyassemblyLib.so
src/polyassembly_findkmers: src/CMakeFiles/polyassembly_findkmers.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/rebecca/work/hifi-potato/paper-code/kmer-counting/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable polyassembly_findkmers"
	cd /home/rebecca/work/hifi-potato/paper-code/kmer-counting/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/polyassembly_findkmers.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/polyassembly_findkmers.dir/build: src/polyassembly_findkmers

.PHONY : src/CMakeFiles/polyassembly_findkmers.dir/build

src/CMakeFiles/polyassembly_findkmers.dir/clean:
	cd /home/rebecca/work/hifi-potato/paper-code/kmer-counting/build/src && $(CMAKE_COMMAND) -P CMakeFiles/polyassembly_findkmers.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/polyassembly_findkmers.dir/clean

src/CMakeFiles/polyassembly_findkmers.dir/depend:
	cd /home/rebecca/work/hifi-potato/paper-code/kmer-counting/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/rebecca/work/hifi-potato/paper-code/kmer-counting /home/rebecca/work/hifi-potato/paper-code/kmer-counting/src /home/rebecca/work/hifi-potato/paper-code/kmer-counting/build /home/rebecca/work/hifi-potato/paper-code/kmer-counting/build/src /home/rebecca/work/hifi-potato/paper-code/kmer-counting/build/src/CMakeFiles/polyassembly_findkmers.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/polyassembly_findkmers.dir/depend
