# Foldseek Test Suite
# Include directories for accessing headers
include_directories(${CMAKE_SOURCE_DIR}/src/commons)
include_directories(${CMAKE_SOURCE_DIR}/lib/mmseqs/src/commons)

# Define test executable
add_executable(test_simple3alphsw TestSimple3AlphSW.cpp)

# Link against foldseek-framework (contains Simple3AlphSW)
# and mmseqs-framework (contains FileUtil, Util, etc.)
target_link_libraries(test_simple3alphsw foldseek-framework mmseqs-framework)

# Add dependency on foldseek-framework
add_dependencies(test_simple3alphsw foldseek-framework)

# Optionally install test (commented out by default)
# install(TARGETS test_simple3alphsw DESTINATION bin)
