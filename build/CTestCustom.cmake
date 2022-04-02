#/*============================================================================
#
#  PHAS0100ASSIGNMENT2: PHAS0100 Assignment 2 Gravitational N-body Simulation
#
#  Copyright (c) University College London (UCL). All rights reserved.
#
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.
#
#  See LICENSE.txt in the top level directory for details.
#
#============================================================================*/

set (CTEST_EXTRA_COVERAGE_GLOB ${CTEST_EXTRA_COVERAGE_GLOB})
foreach (E IN ITEMS .c .cpp .h .txx)
  list (APPEND CTEST_EXTRA_COVERAGE_GLOB
    "Code/*/*${E}"
)
endforeach ()

set(CTEST_CUSTOM_COVERAGE_EXCLUDE
  ${CTEST_CUSTOM_COVERAGE_EXCLUDE}

  # Exclude try_compile sources from coverage results:
  ".*/CMakeFiles/CMakeTmp/.*"

  # Exclude files generated by the moc pre-compiler
  ".*/moc_.*"

  # Exclude files generated by the uic pre-compiler
  ".*/ui_.*"
  )

# The following tests should not be run under valgrind
set(CTEST_CUSTOM_MEMCHECK_IGNORE
  
  )

set(CTEST_CUSTOM_ERROR_MATCH
  ${CTEST_CUSTOM_ERROR_MATCH}
  "CMake Error[ :]"
  )

set(CTEST_CUSTOM_WARNING_MATCH
  ${CTEST_CUSTOM_WARNING_MATCH}
  "CMake Warning[ :]"
  )

set(CTEST_CUSTOM_WARNING_EXCEPTION
  ${CTEST_CUSTOM_WARNING_EXCEPTION}
  
  # kwstyle suppressions
  "[Kk][Ww][Ss]tyle.*kws.*cxx"
  "[Kk][Ww][Ss]tyle.*kws.*h"
  "[Kk][Ww][Ss]tyle.*metaCommand.*cxx"
  )
