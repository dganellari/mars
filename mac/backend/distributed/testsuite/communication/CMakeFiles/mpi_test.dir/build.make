# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/gandanie/scratch/santis/mars

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/gandanie/scratch/santis/mars/mac

# Include any dependencies generated for this target.
include backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/compiler_depend.make

# Include the progress variables for this target.
include backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/progress.make

# Include the compile flags for this target's objects.
include backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/flags.make

backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.o: backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/flags.make
backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.o: /Users/gandanie/scratch/santis/mars/backend/distributed/testsuite/communication/mars_test_mpi.cpp
backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.o: backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/gandanie/scratch/santis/mars/mac/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.o"
	cd /Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.o -MF CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.o.d -o CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.o -c /Users/gandanie/scratch/santis/mars/backend/distributed/testsuite/communication/mars_test_mpi.cpp

backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.i"
	cd /Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/gandanie/scratch/santis/mars/backend/distributed/testsuite/communication/mars_test_mpi.cpp > CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.i

backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.s"
	cd /Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication && /Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/gandanie/scratch/santis/mars/backend/distributed/testsuite/communication/mars_test_mpi.cpp -o CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.s

# Object files for target mpi_test
mpi_test_OBJECTS = \
"CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.o"

# External object files for target mpi_test
mpi_test_EXTERNAL_OBJECTS =

backend/distributed/testsuite/communication/mpi_test: backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/mars_test_mpi.cpp.o
backend/distributed/testsuite/communication/mpi_test: backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/build.make
backend/distributed/testsuite/communication/mpi_test: libmars.a
backend/distributed/testsuite/communication/mpi_test: lib/libgtest.a
backend/distributed/testsuite/communication/mpi_test: lib/libgmock.a
backend/distributed/testsuite/communication/mpi_test: /Users/gandanie/scratch/view_mars/lib/libmpi.dylib
backend/distributed/testsuite/communication/mpi_test: /Users/gandanie/scratch/view_mars/lib/libkokkoskernels.dylib
backend/distributed/testsuite/communication/mpi_test: /Users/gandanie/scratch/view_mars/lib/libkokkoscontainers.4.0.0.dylib
backend/distributed/testsuite/communication/mpi_test: /Users/gandanie/scratch/view_mars/lib/libkokkoscore.4.0.0.dylib
backend/distributed/testsuite/communication/mpi_test: /Users/gandanie/scratch/view_mars/lib/libkokkossimd.4.0.0.dylib
backend/distributed/testsuite/communication/mpi_test: lib/libgtest.a
backend/distributed/testsuite/communication/mpi_test: backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/gandanie/scratch/santis/mars/mac/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable mpi_test"
	cd /Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/mpi_test.dir/link.txt --verbose=$(VERBOSE)
	cd /Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication && /opt/homebrew/Cellar/cmake/3.28.3/bin/cmake -D TEST_TARGET=mpi_test -D TEST_EXECUTABLE=/Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication/mpi_test -D TEST_EXECUTOR= -D TEST_WORKING_DIR=/Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication -D TEST_EXTRA_ARGS= -D TEST_PROPERTIES=VS_DEBUGGER_WORKING_DIRECTORY -D TEST_PREFIX= -D TEST_SUFFIX= -D TEST_FILTER= -D NO_PRETTY_TYPES=FALSE -D NO_PRETTY_VALUES=FALSE -D TEST_LIST=mpi_test_TESTS -D CTEST_FILE=/Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication/mpi_test[1]_tests.cmake -D TEST_DISCOVERY_TIMEOUT=5 -D TEST_XML_OUTPUT_DIR= -P /opt/homebrew/Cellar/cmake/3.28.3/share/cmake/Modules/GoogleTestAddTests.cmake

# Rule to build all files generated by this target.
backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/build: backend/distributed/testsuite/communication/mpi_test
.PHONY : backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/build

backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/clean:
	cd /Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication && $(CMAKE_COMMAND) -P CMakeFiles/mpi_test.dir/cmake_clean.cmake
.PHONY : backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/clean

backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/depend:
	cd /Users/gandanie/scratch/santis/mars/mac && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/gandanie/scratch/santis/mars /Users/gandanie/scratch/santis/mars/backend/distributed/testsuite/communication /Users/gandanie/scratch/santis/mars/mac /Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication /Users/gandanie/scratch/santis/mars/mac/backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : backend/distributed/testsuite/communication/CMakeFiles/mpi_test.dir/depend

