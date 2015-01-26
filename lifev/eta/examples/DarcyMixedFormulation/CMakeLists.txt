INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DarcyMixedFormulation
  SOURCES main.cpp
  ARGS -c
  NUM_MPI_PROCS 2
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_DarcyMixedFormulation
  SOURCE_FILES data
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_DarcyMixedFormulation
  SOURCE_FILES cube.mesh
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/lifev/eta/data/mesh/inria/
)
TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverParamListAztecOO_xml_DiffusionReaction
  SOURCE_FILES SolverParamListAztecOO.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)
TRIBITS_COPY_FILES_TO_BINARY_DIR(SolverParamListBelos_xml_DiffusionReaction
  SOURCE_FILES SolverParamListBelos.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
)


