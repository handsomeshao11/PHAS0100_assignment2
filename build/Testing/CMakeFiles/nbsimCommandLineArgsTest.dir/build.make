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
CMAKE_SOURCE_DIR = /workspaces/PHAS0100Assignment2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /workspaces/PHAS0100Assignment2/build

# Include any dependencies generated for this target.
include Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/depend.make

# Include the progress variables for this target.
include Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/progress.make

# Include the compile flags for this target's objects.
include Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/flags.make

Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCommandLineArgsTest.cpp.o: Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/flags.make
Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCommandLineArgsTest.cpp.o: ../Testing/nbsimCommandLineArgsTest.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/PHAS0100Assignment2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCommandLineArgsTest.cpp.o"
	cd /workspaces/PHAS0100Assignment2/build/Testing && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCommandLineArgsTest.cpp.o -c /workspaces/PHAS0100Assignment2/Testing/nbsimCommandLineArgsTest.cpp

Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCommandLineArgsTest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCommandLineArgsTest.cpp.i"
	cd /workspaces/PHAS0100Assignment2/build/Testing && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/PHAS0100Assignment2/Testing/nbsimCommandLineArgsTest.cpp > CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCommandLineArgsTest.cpp.i

Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCommandLineArgsTest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCommandLineArgsTest.cpp.s"
	cd /workspaces/PHAS0100Assignment2/build/Testing && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/PHAS0100Assignment2/Testing/nbsimCommandLineArgsTest.cpp -o CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCommandLineArgsTest.cpp.s

Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCatchMain.cpp.o: Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/flags.make
Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCatchMain.cpp.o: ../Testing/nbsimCatchMain.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/workspaces/PHAS0100Assignment2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCatchMain.cpp.o"
	cd /workspaces/PHAS0100Assignment2/build/Testing && /usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCatchMain.cpp.o -c /workspaces/PHAS0100Assignment2/Testing/nbsimCatchMain.cpp

Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCatchMain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCatchMain.cpp.i"
	cd /workspaces/PHAS0100Assignment2/build/Testing && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /workspaces/PHAS0100Assignment2/Testing/nbsimCatchMain.cpp > CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCatchMain.cpp.i

Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCatchMain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCatchMain.cpp.s"
	cd /workspaces/PHAS0100Assignment2/build/Testing && /usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /workspaces/PHAS0100Assignment2/Testing/nbsimCatchMain.cpp -o CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCatchMain.cpp.s

# Object files for target nbsimCommandLineArgsTest
nbsimCommandLineArgsTest_OBJECTS = \
"CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCommandLineArgsTest.cpp.o" \
"CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCatchMain.cpp.o"

# External object files for target nbsimCommandLineArgsTest
nbsimCommandLineArgsTest_EXTERNAL_OBJECTS =

bin/nbsimCommandLineArgsTest: Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCommandLineArgsTest.cpp.o
bin/nbsimCommandLineArgsTest: Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/nbsimCatchMain.cpp.o
bin/nbsimCommandLineArgsTest: Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/build.make
bin/nbsimCommandLineArgsTest: bin/libphas0100assignment2.a
bin/nbsimCommandLineArgsTest: Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/workspaces/PHAS0100Assignment2/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable ../bin/nbsimCommandLineArgsTest"
	cd /workspaces/PHAS0100Assignment2/build/Testing && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/nbsimCommandLineArgsTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/build: bin/nbsimCommandLineArgsTest

.PHONY : Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/build

Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/clean:
	cd /workspaces/PHAS0100Assignment2/build/Testing && $(CMAKE_COMMAND) -P CMakeFiles/nbsimCommandLineArgsTest.dir/cmake_clean.cmake
.PHONY : Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/clean

Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/depend:
	cd /workspaces/PHAS0100Assignment2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /workspaces/PHAS0100Assignment2 /workspaces/PHAS0100Assignment2/Testing /workspaces/PHAS0100Assignment2/build /workspaces/PHAS0100Assignment2/build/Testing /workspaces/PHAS0100Assignment2/build/Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Testing/CMakeFiles/nbsimCommandLineArgsTest.dir/depend

