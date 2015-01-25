INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RegionMarkerID1D
  SOURCES basic_test1D.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RegionMarkerID2D
  SOURCES basic_test2D.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  RegionMarkerID3D
  SOURCES basic_test3D.cpp
  ARGS -c
  NUM_MPI_PROCS 1
  COMM serial mpi
#  STANDARD_PASS_OUTPUT
  )


TRIBITS_COPY_FILES_TO_BINARY_DIR(data_region_marker_id
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(data_region_marker_id
  SOURCE_FILES square20x20.msh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/core/data/mesh/freefem/
)

