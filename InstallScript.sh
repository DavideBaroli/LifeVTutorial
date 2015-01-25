
LifeV_DIR=$HOME
LifeV_SRC=$LifeV_DIR/LifeV-src
LifeV_BUILD=$LifeV_DIR/LifeV-build
LifeV_LIB=$LifeV_DIR/LifeV-install
mkdir  $LifeV_SRC
mkdir  $LifeV_BUILD
mkdir  $LifeV_LIB

parmetisLib="${mkParmetisLib};${mkMetisLib}"
parmetisInc="${mkParmetisInc};${mkMetisInc}"

HOMESCRIPT=$PWD

cd $LifeV_DIR
${mkGitBin}/git clone https://github.com/DavideBaroli/LifeVTutorial.git LifeV-src

cd $LifeV_BUILD

     rm  -rf   MakeCache.txt
     rm -rf CMakeFiles

     $mkCmakeBin/cmake \
    -D BUILD_SHARED_LIBS:BOOL=ON \
    -D CMAKE_BUILD_TYPE:STRING=RELEASE \
    -D CMAKE_INSTALL_PREFIX:PATH=${LifeV_LIB}\
    -D CMAKE_C_COMPILER:STRING="$mkOpenmpiBin/mpicc" \
    -D CMAKE_CXX_COMPILER:STRING="$mkOpenmpiBin/mpicxx" \
    -D CMAKE_CXX_FLAGS:STRING="-O3 -msse3 -ansi" \
    -D CMAKE_C_FLAGS:STRING="-O3 -msse3 -ansi" \
    -D CMAKE_Fortran_COMPILER:STRING="$mkOpenmpiBin/mpif90" \
    -D CMAKE_Fortran_FLAGS:STRING="-Og -g" \
    -D CMAKE_AR:STRING="/usr/bin/ar" \
    -D CMAKE_MAKE_PROGRAM:STRING="/usr/bin/gmake" \
       \
    -D  TPL_ENABLE_AMD=ON \
     -D      AMD_INCLUDE_DIRS=${mkSuitesparseInc} \
     -D       AMD_LIBRARY_DIRS=${mkSuitesparseLib} \
     -D       AMD_LIBRARY_NAMES=amd \
     -D    TPL_ENABLE_BLAS=ON \
     -D       BLAS_INCLUDE_DIRS=${mkBlasInc} \
     -D       BLAS_LIBRARY_DIRS=${mkBlasLib} \
     -D       BLAS_LIBRARY_NAMES=${mkBlas} \
     -D    TPL_ENABLE_Boost=ON \
     -D       Boost_INCLUDE_DIRS=${mkBoostInc} \
     -D   TPL_ENABLE_BoostLib=ON \
     -D Boost_NO_BOOST_CMAKE:BOOL=ON \
     -D       BoostLib_INCLUDE_DIRS=${mkBoostInc} \
     -D       BoostLib_LIBRARY_DIRS=${mkBoostLib} \
     -D   TPL_ENABLE_HDF5=ON \
     -D       HDF5_INCLUDE_DIRS=${mkHdf5Inc} \
     -D       HDF5_LIBRARY_DIRS=${mkHdf5Lib} \
     -D   TPL_ENABLE_LAPACK=ON \
     -D       LAPACK_INCLUDE_DIRS=${mkLapackInc} \
     -D       LAPACK_LIBRARY_DIRS=${mkLapackLib} \
     -D       LAPACK_LIBRARY_NAMES=${mkLapack} \
     -D   TPL_ENABLE_MPI=ON \
     -D       MPI_BASE_DIR:PATH=$mkOpenmpiHome \
     -D       ParMETIS_INCLUDE_DIRS=$parmetisInc \
     -D       ParMETIS_LIBRARY_DIRS=$parmetisLib \
     -D   TPL_ENABLE_UMFPACK=ON \
     -D       UMFPACK_INCLUDE_DIRS=${mkSuitesparseInc} \
     -D      UMFPACK_LIBRARY_DIRS=${mkSuitesparseLib} \
     -D       UMFPACK_LIBRARY_NAMES=umfpack \
     -D   Trilinos_INCLUDE_DIRS:PATH="${mkTrilinosInc}" \
     -D   Trilinos_LIBRARY_DIRS:PATH="${mkTrilinosLib}" \
     -D   LifeV_ENABLE_TESTS:BOOL=ON \
     -D   LifeV_ENABLE_ALL_PACKAGES:BOOL=ON \
     -D   LifeV_ENABLE_STRONG_CXX_COMPILE_WARNINGS:BOOL=ON \
     -D   LifeV_ENABLE_CPP11:BOOL=ON\
       ${LifeV_SRC}
cd $LifeV_BUILD
make


cd $HOMESCRIPT
