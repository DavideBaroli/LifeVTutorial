//@HEADER
/*
*******************************************************************************

    Copyright (C) 2004, 2005, 2007 EPFL, Politecnico di Milano, INRIA
    Copyright (C) 2010 EPFL, Politecnico di Milano, Emory University

    This file is part of LifeV.

    LifeV is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    LifeV is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with LifeV.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************
*/
//@HEADER

/*!
    @file
    @brief Diffusion-Reaction Scalar and Vectorial Problem

    @author Davide Baroli <davide.baroli@polimi.it>
    @date 2015

    ETA stands Expression Template Assembly, in reference
    to the metaprogramming technique used.

 */

// ---------------------------------------------------------------
// We include here the MPI headers for the parallel computations.
// ---------------------------------------------------------------


#include <Epetra_ConfigDefs.h>
// ---------------------------------------------------------------
// Include Teuchos XML reader and Trilinos Smart Pointer Definition
//
// ---------------------------------------------------------------
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>



#ifdef EPETRA_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

// -----------------------------------
// STL header
// ----------------------------------
#include <memory>
#include <iostream>

// ---------------------------------------------------------------
// We include then the required headers from LifeV. First of all,
// the definition file and mesh related files. We also include
// the MatrixEpetra since this is the kind of object that we want
// to assemble.
// ---------------------------------------------------------------

#include <lifev/core/LifeV.hpp>

// --------------------------------------------------------------
// LifeV Mesh Utility, Partitioner, Hypercube structured generator
// --------------------------------------------------------------

#include <lifev/core/mesh/MeshPartitioner.hpp>
#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/RegionMesh.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshData.hpp>

// ===============================================
//  LifeV GetPot - Reader from dataFile
//=================================================
#include <lifev/core/filter/GetPot.hpp>
// ===============================================
//  LifeV Chrono
//=================================================

#include <lifev/core/util/LifeChrono.hpp>


// -------------------------------------------------------
// LifeV Linear Algebra
// -------------------------------------------------------

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>
// -------------------------------------------------------
// LifeV Preconditioner(ML, Ifpack)  and Solver (Belos)
// --------------------------------------------------------


#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
//#include <lifev/core/algorithm/SolverAztecOO.hpp>
#include <lifev/core/algorithm/LinearSolver.hpp>

// --------------------------------------------------
// LifeV Exporter
// ----------------------------------------------

#include <lifev/core/filter/ExporterEnsight.hpp>
#include <lifev/core/filter/ExporterHDF5.hpp>


// ---------------------------------------------------------------
// In order to use the ETA framework, a special version of the
// FESpace structure must be used. It is called ETFESpace and
// has basically the same role as the FESpace.
// ---------------------------------------------------------------

#include <lifev/eta/fem/ETFESpace.hpp>
#include <lifev/core/fem/FESpace.hpp>
/// We can Choose Between Scalar FE and Vectorial FE

//#define SCALARFE 0
#define VECTORIALFE 0

//--------------------------------------------------
//LifeV  Boundary Conditions
//------------------------------------------------
#include <lifev/core/fem/BCHandler.hpp>
#include <lifev/core/fem/BCManage.hpp>

// -----------------------------------------------
// Post-Processing Utility
// ----------------------------------------------

#include <lifev/core/fem/GradientRecovery.hpp>


// ---------------------------------------------------------------
// The most important file to include is the Integrate.hpp file
// which contains all the definitions required to perform the
// different integrations.
// ---------------------------------------------------------------

#include <lifev/eta/expression/Integrate.hpp>
#include <lifev/eta/fem/QuadratureBoundary.hpp>

// ---------------------------------------------------------------
// Finally, we include shared pointer from boost since we use
// them explicitly in this tutorial.
// ---------------------------------------------------------------

#include <boost/shared_ptr.hpp>


// ---------------------------------------------------------------
// As usual, we work in the LifeV namespace. For clarity, we also
// make two typedefs for the mesh type and matrix type.
// ---------------------------------------------------------------

using namespace LifeV;
/* Register the preconditioner */
namespace
{
static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );
}

typedef RegionMesh<LinearTetra> mesh_Type;
typedef MatrixEpetra<Real> matrix_Type;
typedef boost::shared_ptr<matrix_Type> matrixPtr_Type;

typedef VectorEpetra vector_Type;
typedef boost::shared_ptr<vector_Type> vectorPtr_Type;

typedef LifeV::Preconditioner             basePrec_Type;
typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;

typedef LifeV::PreconditionerIfpack       prec_Type;
typedef boost::shared_ptr<prec_Type>      precPtr_Type;


const Real pi = 3.141592653589793;

// --------------------------------------------------------------
// We define the permeability coefficients which vary spatially
// ---------------------------------------------------------------
// Change Permeability...

Real scalarpermeability ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 1.; // std::sin (2* pi / y ) * std::cos ( 2*pi / x ) * std::exp ( z ) ;
}

Real scalarpermeability11 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 2; //  std::sin (2* pi / y ) * std::cos ( 2*pi / x ) * std::exp ( z ) ;
}

Real scalarpermeability12 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 0.0;
}

Real scalarpermeability13 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 0.0;
}
Real scalarpermeability21   ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 0.0;
}
Real scalarpermeability22   ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 1.; // std::sin (2* pi / y ) * std::cos ( 2*pi / x ) * std::exp ( z ) ;
}
Real scalarpermeability23   ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 0.0;
}

Real scalarpermeability31   ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 0.0;
}
Real scalarpermeability32   ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 0.0 ;
}
Real scalarpermeability33   ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 1.; // std::sin (2* pi / y ) * std::cos ( 2*pi / x ) * std::exp ( z );
}

// Source Rhs
Real ScalarRhs ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return  1* (std::sqrt( std::pow(x,2)+std::pow(y,2)+std::pow(z,2) )<1 ) ;
}

// Component Rhs
Real ScalarRhs0 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return  1* (std::sqrt( std::pow(x,2)+std::pow(y,2)+std::pow(z,2) )<1 ) ; //  x;
}
// Component Rhs
Real ScalarRhs1 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return  1* (std::sqrt( std::pow(x,2)+std::pow(y,2)+std::pow(z,2) )<1 ) ; // x;

}

// Component Rhs
Real ScalarRhs2 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return  1* (std::sqrt( std::pow(x,2)+std::pow(y,2)+std::pow(z,2) )<1 ) ; // x;
}

Real scalarReaction( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
return x*y-0.5*z ;
}

// Boundary Conditions

Real dirichlet ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
return  x;
}

Real neumann1 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{

    VectorSmall<3> gradient;
    gradient[0]= -1. * ( 4. * x * y * y + 2. * x * x * y + 12. - 2. * x * x * x - 2. * y );
    gradient[1]=-1. * ( 2. * y * x * x + 2. * x * y * y + 6. - 2. * y - x * x * x );
    gradient[2]=-1. * ( 5. - 4. * z );

    VectorSmall<3> normalz;
    normalz[0] = 0;
    normalz[1] = 0;
    normalz[2] = 1;
return gradient.dot(normalz);
}


Real neumann2 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{

    VectorSmall<3> gradient;
    gradient[0]= -1. * ( 4. * x * y * y + 2. * x * x * y + 12. - 2. * x * x * x - 2. * y );
    gradient[1]=-1. * ( 2. * y * x * x + 2. * x * y * y + 6. - 2. * y - x * x * x );
    gradient[2]=-1. * ( 5. - 4. * z );

    VectorSmall<3> normalx;
    normalx[0] = 1;
    normalx[1] = 0;
    normalx[2] = 0;
return -1.*gradient.dot(normalx);
}


Real robin ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
return  neumann1(0,x,y,z,0)-neumann2(0,x,y,z,0) ;
}


// Functor

class ScalarPermeabilityFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
        return scalarpermeability (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  ) ;
    }

    ScalarPermeabilityFunctor() {}
   ScalarPermeabilityFunctor (const ScalarPermeabilityFunctor&) {}
     ~ScalarPermeabilityFunctor() {}
};

class MatrixPermeabilityFunctor
{
public:
    typedef MatrixSmall<3,3>  return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
         MatrixSmall<3,3> Perm;
         Perm[0][0]= scalarpermeability11 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[0][1]= scalarpermeability12 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[0][2]= scalarpermeability13 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[1][0]= scalarpermeability21 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[1][1]= scalarpermeability22 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[1][2]= scalarpermeability23 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[2][0]= scalarpermeability31 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[2][1]= scalarpermeability32 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[2][2]= scalarpermeability33 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );

       return Perm;
   }

    MatrixPermeabilityFunctor() {}
    MatrixPermeabilityFunctor (const MatrixPermeabilityFunctor&) {}
     ~MatrixPermeabilityFunctor() {}
};

/* LaplacianRhs */
class DarcyRhsScalarFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
        return ScalarRhs( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  ) ;
    }

    DarcyRhsScalarFunctor() {}
    DarcyRhsScalarFunctor (const DarcyRhsScalarFunctor&) {}
    ~DarcyRhsScalarFunctor() {}
};


/* LaplacianRhs */
class DarcyRhsVectorialFunctor
{
public:
    typedef VectorSmall<3>  return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
       return_Type RhsVector;
       RhsVector[0]=  ScalarRhs0( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  ) ;
       RhsVector[1]=  ScalarRhs1( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  ) ;
       RhsVector[2]=  ScalarRhs2( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  ) ;

       return RhsVector;
    }

    DarcyRhsVectorialFunctor() {}
    DarcyRhsVectorialFunctor (const DarcyRhsVectorialFunctor&) {}
    ~DarcyRhsVectorialFunctor() {}
};

/* LaplacianNeumann1 */
class DarcyNeumann1VectorialFunctor
{
public:
    typedef VectorSmall<3>  return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {

       return_Type RhsNeumannVector;

       RhsNeumannVector[0]=0.  ;
       RhsNeumannVector[1]=0.  ;
       RhsNeumannVector[2]=neumann1( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
       return RhsNeumannVector;
    }

    DarcyNeumann1VectorialFunctor() {}
    DarcyNeumann1VectorialFunctor (const DarcyNeumann1VectorialFunctor&) {}
    ~DarcyNeumann1VectorialFunctor() {}
};
/* LaplacianNeumann2 */
class DarcyNeumann2VectorialFunctor
{
public:
    typedef VectorSmall<3>  return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {

       return_Type RhsNeumannVector;
       RhsNeumannVector[0]=neumann2( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
        RhsNeumannVector[1]=0.  ;
       RhsNeumannVector[2]=0.  ;
      return RhsNeumannVector;
    }

    DarcyNeumann2VectorialFunctor() {}
    DarcyNeumann2VectorialFunctor (const DarcyNeumann2VectorialFunctor&) {}
    ~DarcyNeumann2VectorialFunctor() {}
};



/* LaplacianNeumann1 */
class DarcyNeumann1ScalarFunctor
{
public:
    typedef Real  return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {

      return neumann1( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
     }

    DarcyNeumann1ScalarFunctor() {}
    DarcyNeumann1ScalarFunctor (const DarcyNeumann1ScalarFunctor&) {}
    ~DarcyNeumann1ScalarFunctor() {}
};
/* LaplacianNeumann2 */
class DarcyNeumann2ScalarFunctor
{
public:
    typedef Real  return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {

      return neumann2( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
     }

    DarcyNeumann2ScalarFunctor() {}
    DarcyNeumann2ScalarFunctor (const DarcyNeumann2ScalarFunctor&) {}
    ~DarcyNeumann2ScalarFunctor() {}
};



class scalarReactionFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
        return scalarReaction( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  ) ;
    }

    scalarReactionFunctor() {}
     scalarReactionFunctor (const  scalarReactionFunctor&) {}
    ~scalarReactionFunctor() {}
};


class robinFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
        return robin ( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  ) ;
    }

    robinFunctor() {}
    robinFunctor (const robinFunctor&) {}
    ~robinFunctor() {}
};

class robinVectFunctor
{
public:
    typedef VectorSmall<3> return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
    VectorSmall<3> robinVect;
    robinVect[0]= robin ( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  ) ;
    robinVect[1]= robin ( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 1  ) ;
    robinVect[2]= robin ( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 2  ) ;

    return robinVect;
    }

    robinVectFunctor() {}
    robinVectFunctor (const robinVectFunctor&) {}
    ~robinVectFunctor() {}
};



// Robin Parameter
Real beta ( 1.0 );
Real alpha ( 1.0 );


// =======================================================
//               Boundary Data ID
// ======================================================
namespace BCFlags
{
// Flags for structured meshes


// Walls
const Int LEFTWALL   = 4;
const Int RIGHTWALL  = 2;
const Int FRONTWALL  = 1;
const Int BACKWALL   = 3;
const Int TOPWALL    = 6;
const Int BOTTOMWALL = 5;
// Edges
const Int BOTTOMEDGE1 =  7;
const Int BOTTOMEDGE2 =  8;
const Int BOTTOMEDGE3 =  9;
const Int BOTTOMEDGE4 = 10;
const Int SIDEEDGE1   = 11;
const Int SIDEEDGE2   = 12;
const Int SIDEEDGE3   = 13;
const Int SIDEEDGE4   = 14;
const Int TOPEDGE1    = 15;
const Int TOPEDGE2    = 16;
const Int TOPEDGE3    = 17;
const Int TOPEDGE4    = 18;
// Corners
const Int BOTTOMCORNER1 = 19;
const Int BOTTOMCORNER2 = 20;
const Int BOTTOMCORNER3 = 21;
const Int BOTTOMCORNER4 = 22;
const Int TOPCORNER1    = 23;
const Int TOPCORNER2    = 24;
const Int TOPCORNER3    = 25;
const Int TOPCORNER4    = 26;
}



// ---------------------------------------------------------------
// We start the programm by the definition of the communicator
// (as usual) depending on whether MPI is available or not. We
// also define a boolean to allow only one process to display
// messages.
// ---------------------------------------------------------------

int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
    // Only local Processor
    //   boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_SELF) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    const bool verbose (Comm->MyPID() == 0);


    if (verbose)
    {
        std::cout << " Solve  Div( Grad( K u)) - g_reaction*u  = - f(u)    ... " << std::flush;
    }


    // ==============================================================
    // Read from DataFile
    // =============================================================
    // Open and read the data file
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile ( data_file_name );




    // ---------------------------------------------------------------
    // The next step is to build the mesh. We use here a structured
    // cartesian mesh over the  domain (startx,startx+lengthMeshx)x(starty,starty+lengthMeshy)x(startz,startz+lengthMeshz).
    // The mesh is the partitioned for the parallel computations and
    // the original mesh is deleted.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building and partitioning the mesh ... " << std::flush;
    }

    const UInt Nelementsx=dataFile("mesh/numElementalonx "  ,10);
    const UInt Nelementsy=dataFile("mesh/numElementalony "  ,10);
    const UInt Nelementsz=dataFile("mesh/numElementalonz "  ,10);
    const Real l_x =dataFile("mesh/lengthMeshx "  , 2.);
    const Real l_y =dataFile("mesh/lengthMeshy "  , 2.);
    const Real l_z =dataFile("mesh/lengthMeshz "  , 2.);

    const Real t_x =dataFile("mesh/startx "  , 2.);
    const Real t_y =dataFile("mesh/starty "  , 2.);
    const Real t_z =dataFile("mesh/startz "  , 2.);



    boost::shared_ptr< mesh_Type > fullMeshPtr (new mesh_Type ( Comm ) );

    regularMesh3D ( *fullMeshPtr, 1, Nelementsx, Nelementsy, Nelementsz, false,
                l_x ,   l_y,   l_z,
                t_x ,   t_y,      t_z  );

    // ==============================================================
    // Read from Mesh ( conforming with LifeV marker , see testsuite/mesh )
    // =============================================================
    // Load the mesh
 /*
    MeshData dataMesh;
    dataMesh.setup (dataFile, "mesh");
    boost::shared_ptr < mesh_Type > fullMeshPtr (new mesh_Type);
    readMesh (*fullMeshPtr, dataMesh);
 */


    boost::shared_ptr< mesh_Type > meshPtr;
    {
        MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, Comm);
        meshPtr = meshPart.meshPartition();
    }

    fullMeshPtr.reset();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We define now the ETFESpace that we need for the assembly.
    // Remark that we use a shared pointer because other structures
    // will require this ETFESpace to be alive. We can also observe
    // that the ETFESpace has more template parameters than the
    // classical FESpace (this is the main difference). The 3
    // indicates that the problem is in 3D while the 1 indicate that
    // the unknown is scalar.
    //
    // After having constructed the ETFESpace, we display the number
    // of degrees of freedom of the problem.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Building ETFESpaces ... " << std::flush;
    }


  std::string  uOrder=dataFile("fespace/order "  ,"P1");

   #ifdef SCALARFE
   boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uFESpace ( new FESpace< mesh_Type, MapEpetra > (meshPtr, uOrder, 1, Comm) );

    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > ETuFESpace
  ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPtr, &(uFESpace->refFE()), Comm) );
   #endif
   #ifdef VECTORIALFE
   boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uFESpace ( new FESpace< mesh_Type, MapEpetra > (meshPtr, uOrder, 3, Comm) );

  boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > ETuFESpace
    ( new ETFESpace< mesh_Type, MapEpetra, 3, 3 > (meshPtr,
      & (uFESpace->refFE()), Comm) );
   #endif

if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }
    if (verbose)
    {
        std::cout << " ---> Dofs: " << ETuFESpace->dof().numTotalDof() << std::endl;
    }


    // ---------------------------------------------------------------
    // The matrix is then defined using the map of the FE space.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Defining the matrix ... " << std::flush;
    }

    QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria4pt) );


    boost::shared_ptr<matrix_Type> SystemMatrix (new matrix_Type ( ETuFESpace->map() ) );


    *SystemMatrix *= 0.0;

    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }


    // ---------------------------------------------------------------
    // We start now the assembly of the matrix.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Assembling the Diffusion-Reaction Equation   ... " << std::flush;
    }


    // ---------------------------------------------------------------
    // To use the ETA framework, it is mandatory to use a special
    // namespace, called ExpressionAssembly. This namespace is useful
    // to avoid collisions with keywords used for the assembly. A
    // special scope is opened to keep only that part of the code
    // in the ExpressionAssembly namespace.
    // ---------------------------------------------------------------
    boost::shared_ptr<ScalarPermeabilityFunctor> ScalarPermeabilityFct ( new ScalarPermeabilityFunctor );
    boost::shared_ptr<scalarReactionFunctor> ScalarReactionFct ( new scalarReactionFunctor );



    boost::shared_ptr<MatrixPermeabilityFunctor> MatrixPermeabilityFct ( new MatrixPermeabilityFunctor );

    {
        using namespace ExpressionAssembly;

        // ---------------------------------------------------------------
        // We can now proceed with assembly. The next instruction
        // assembles the stiffness and mass matrix.
        //
        // The first argument of the integrate function indicates that the
        // integration is done on the elements of the mesh located in the
        // ETFESpace defined earlier.
        //
        // The second argument is simply the quadrature rule to be used.
        //
        // The third argument is the finite element space of the test
        // functions.
        //
        // The fourth argument is the finite element space of the trial
        // functions (those used to represent the solution).
        //
        // The last argument is the expression to be integrated, i.e.
        // that represents the weak formulation of the problem. The
        // keyword phi_i stands for a generic test function and phi_j
        // a generic trial function. The function grad applied to them
        // indicates that the gradient is considered and the dot function
        // indicates a dot product between the two gradients. The
        // expression to be integrated is then the dot product between
        // the gradient of the test function and the gradient of the trial
        // function. This corresponds to the left hand side of the weak
        // formulation of the Laplace problem.
        //
        // Finally, the operator >> indicates that the result of the
        // integration must be added to the SystemMatrix.
        // ---------------------------------------------------------------

    #ifdef SCALARFE

integrate (  elements (ETuFESpace->mesh() ),
                     quadRuleTetra4pt,
                     ETuFESpace,
                     ETuFESpace,
    eval( ScalarPermeabilityFct, X) * dot ( grad (phi_j) , grad (phi_i) )
    - eval( ScalarReactionFct, X) *   phi_j * phi_i

                  )
                >> *SystemMatrix;

// - \int _\GammaN \grad u cdot n =- \int_GammaN gN v + \int_\GammaN  alpha u cdot v

// Assembly of LHS of Robin Boundary Condition

integrate (  boundary (ETuFESpace->mesh(), BCFlags::BOTTOMWALL ),
                     myBDQR,
                     ETuFESpace,
                     ETuFESpace,
            value ( alpha ) * phi_j * phi_i

 )
                >> *SystemMatrix;




    #endif

    #ifdef VECTORIALFE
integrate (  elements (ETuFESpace->mesh() ),
                     quadRuleTetra4pt,
                     ETuFESpace,
                     ETuFESpace,
          dot( eval(MatrixPermeabilityFct, X)* grad (phi_j) , grad (phi_i) )
          -    eval( ScalarReactionFct, X) * dot(   phi_j , phi_i)


                  )
                >> *SystemMatrix;



integrate (  boundary (ETuFESpace->mesh(), BCFlags::BOTTOMWALL ),
                     myBDQR,
                     ETuFESpace,
                     ETuFESpace,
            value ( alpha )* dot( phi_j , phi_i)

 )
                >> *SystemMatrix;

#endif
   }

    if (verbose)
    {
        std::cout << " Assembling the rhs ... " << std::flush;
    }
   boost::shared_ptr<DarcyRhsScalarFunctor> ScalarFctRhs ( new DarcyRhsScalarFunctor );
   boost::shared_ptr<DarcyNeumann1ScalarFunctor> ScalarFctN1 ( new DarcyNeumann1ScalarFunctor );
   boost::shared_ptr<DarcyNeumann2ScalarFunctor> ScalarFctN2 ( new DarcyNeumann2ScalarFunctor );

   boost::shared_ptr<robinFunctor> ScalarFctR ( new robinFunctor );


     boost::shared_ptr<DarcyRhsVectorialFunctor> VectorialFctRhs ( new DarcyRhsVectorialFunctor );

   boost::shared_ptr<DarcyNeumann1VectorialFunctor> VectorialFctN1 ( new DarcyNeumann1VectorialFunctor );
   boost::shared_ptr<DarcyNeumann2VectorialFunctor> VectorialFctN2 ( new DarcyNeumann2VectorialFunctor );

   boost::shared_ptr<robinVectFunctor> VectorialFctR ( new robinVectFunctor );


   vector_Type uRhs ( ETuFESpace->map() , Repeated );
    uRhs *= 0.0;

    {
        using namespace ExpressionAssembly;
    #ifdef SCALARFE
        integrate ( elements (ETuFESpace->mesh()), // Mesh


                   uFESpace->qr(), // QR

                    ETuFESpace,

                     eval ( ScalarFctRhs, X ) * phi_i

                  )
                  >> uRhs;
        integrate ( boundary (ETuFESpace->mesh(),BCFlags::TOPWALL), // Mesh

                    myBDQR, // QR

                    ETuFESpace,

                     eval ( ScalarFctN1 , X ) * phi_i

                  )
                  >>uRhs;
        integrate ( boundary (ETuFESpace->mesh(),BCFlags::FRONTWALL), // Mesh

                    myBDQR, // QR

                    ETuFESpace,

                     eval ( ScalarFctN2 , X ) * phi_i

                  )
                  >>uRhs;

        integrate ( boundary (ETuFESpace->mesh(),BCFlags::BOTTOMWALL), // Mesh

                    myBDQR, // QR

                    ETuFESpace,

                     eval ( ScalarFctR , X ) * phi_i

                  )
                  >>uRhs;

    #endif
    #ifdef VECTORIALFE
        integrate ( elements (ETuFESpace->mesh() ), // Mesh

                    uFESpace->qr(), // QR

                    ETuFESpace,

                   dot(  eval ( VectorialFctRhs, X ) , phi_i)

                  )
                >> uRhs;

        integrate ( boundary (ETuFESpace->mesh(),BCFlags::TOPWALL), // Mesh

                    myBDQR, // QR

                    ETuFESpace,

                     dot(  eval ( VectorialFctN1 , X ) , phi_i)

                  )
                  >>uRhs;

     integrate ( boundary (ETuFESpace->mesh(),BCFlags::FRONTWALL), // Mesh

                    myBDQR, // QR

                    ETuFESpace,

                     dot(  eval ( VectorialFctN2 , X ) , phi_i)

                  )
                  >>uRhs;

     integrate ( boundary (ETuFESpace->mesh(),BCFlags::BOTTOMWALL), // Mesh

                    myBDQR, // QR

                    ETuFESpace,

                     dot(  eval ( VectorialFctR , X ) , phi_i)

                  )
                  >>uRhs;


  #endif
   }



    if (verbose)
    {
        std::cout << " done! " << std::endl;
    }
    // ---------------------------------------------------------------
    // As we are already done with the assembly of the matrix, we
    // finalize it to be able to work on it, e.g. to solve a linear
    // system.
    // ---------------------------------------------------------------

    if (verbose)
    {
        std::cout << " -- Closing the matrix ... " << std::flush;
    }
    vectorPtr_Type uRhsUnique ( new vector_Type ( uRhs, Unique ) );

    SystemMatrix->globalAssemble();

    if (verbose)
    {
        std::cout << " done ! " << std::endl;
    }

    if (verbose)
    {
        std::cout << " Applying Dirichlet boundary conditions ... " << std::flush;
    }

   BCHandler bcHandler;

    BCFunctionBase dirichletBCFct ( dirichlet );

//======================================================
// Essential BC are identified by  coordinates x; y; z
//=====================================================

   #ifdef SCALARFE

    bcHandler.addBC ("Left", BCFlags::LEFTWALL, Essential, Full, dirichletBCFct, 1);
    bcHandler.addBC ("Right", BCFlags::RIGHTWALL, Essential, Full , dirichletBCFct, 1);
    bcHandler.addBC ("Back", BCFlags::BACKWALL, Essential, Full, dirichletBCFct, 1);

    #endif

    #ifdef VECTORIALFE
    bcHandler.addBC ("Left", BCFlags::LEFTWALL, Essential, Full , dirichletBCFct, 3);
    bcHandler.addBC ("Right", BCFlags::RIGHTWALL, Essential, Full , dirichletBCFct, 3);
    bcHandler.addBC ("Back", BCFlags::BACKWALL, Essential, Full , dirichletBCFct, 3);
    #endif

    bcHandler.bcUpdate ( *meshPtr, uFESpace->feBd(), uFESpace->dof() );
    bcManage (*SystemMatrix, *uRhsUnique,
              *uFESpace->mesh(), uFESpace->dof(),
              bcHandler, uFESpace->feBd(), 1.0, Real (0.0) );




    vectorPtr_Type uSolution ( new vector_Type ( ETuFESpace->map() , Unique) );


   if ( verbose )
    {
        std::cout << "Setting up LinearSolver (Belos)... " << std::flush;
    }

    prec_Type* precRawPtr;
    basePrecPtr_Type precPtr;
    precRawPtr = new prec_Type;
    precRawPtr->setDataFromGetPot ( dataFile, "prec" );
    precPtr.reset ( precRawPtr );

    Teuchos::RCP< Teuchos::ParameterList > belosList = Teuchos::rcp ( new Teuchos::ParameterList );
    belosList = Teuchos::getParametersFromXmlFile ( "SolverParamList.xml" );

    LinearSolver linearSolver;
    linearSolver.setCommunicator ( Comm );
    linearSolver.setParameters ( *belosList );
    linearSolver.setPreconditioner ( precPtr );
    if ( verbose )
    {
        std::cout << "done" << std::endl;
    }

   if ( verbose )
    {
        std::cout << " Solving the system " << std::endl;
    }

    linearSolver.setOperator ( SystemMatrix );
    linearSolver.setRightHandSide ( uRhsUnique );
    linearSolver.solve ( uSolution );




    if (verbose)
    {
        std::cout << " Building the exporter " << std::endl;
    }

    std::string const exporterFileName    =  dataFile ( "exporter/filename", "cube");
    ExporterHDF5<mesh_Type> exporter ( dataFile, meshPtr, exporterFileName, Comm->MyPID() );
    exporter.setMultimesh (false);

    boost::shared_ptr<vector_Type> uExported ( new vector_Type (ETuFESpace->map(), exporter.mapType() ) );
#ifdef SCALARFE
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "u", uFESpace, uExported, UInt (0) );

#endif
#ifdef VECTORIALFE
    exporter.addVariable ( ExporterData<mesh_Type>::VectorField, "u", uFESpace, uExported, UInt (0) );
#endif



    if ( verbose )
    {
        std::cout << "done" << std::endl;
    }

// Assign to uExporter the Solution

    *uExported = *uSolution;

    exporter.postProcess ( 1.0 );


// Close Exporter
    exporter.closeFile();


   bool success ( true );

   if (!success)
    {
        if (verbose)
        {
            std::cout << "End Result: TEST NOT PASSED" << std::endl;
        }
    }
    else
    {
        if (verbose)
        {
            std::cout << "End Result: TEST PASSED" << std::endl;
        }
    }





    // ---------------------------------------------------------------
    // We finalize the MPI session if MPI was used
    // ---------------------------------------------------------------

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    if ( !success )
    {
        return ( EXIT_FAILURE );
    }
    else
    {
        return ( EXIT_SUCCESS );
    }


}


