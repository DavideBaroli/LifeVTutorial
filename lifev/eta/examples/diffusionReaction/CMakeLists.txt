INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DiffusionReactionExample
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ScalarVectorial
  SOURCES mainSV.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  )



TRIBITS_COPY_FILES_TO_BINARY_DIR(dataETA_DiffusionReaction
  CREATE_SYMLINK
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverParamListAztecOO_xml_DiffusionReaction
  SOURCE_FILES SolverParamListAztecOO.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)
TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverParamListBelos_xml_DiffusionReaction
  SOURCE_FILES SolverParamListBelos.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(cube.mesh_test
  SOURCE_FILES cube.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/eta/data/mesh/inria/
)

