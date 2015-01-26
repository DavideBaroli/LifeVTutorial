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
    @brief Mixed Darcy Formulation
   @author dbaroli
   @date 2015
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
#include <lifev/core/util/LifeChronoManager.hpp>


// -------------------------------------------------------
// LifeV Linear (Block) Algebra
// -------------------------------------------------------

#include <lifev/core/array/MatrixEpetra.hpp>
#include <lifev/core/array/MapEpetra.hpp>


#include <lifev/core/array/MatrixEpetraStructured.hpp>
#include <lifev/core/array/VectorBlockMonolithicEpetra.hpp>



// -------------------------------------------------------
// LifeV Preconditioner(ML, Ifpack)  and Solver (Belos)
// --------------------------------------------------------


#include <lifev/core/algorithm/Preconditioner.hpp>
#include <lifev/core/algorithm/PreconditionerIfpack.hpp>
#include <lifev/core/algorithm/PreconditionerML.hpp>
#include <lifev/core/algorithm/SolverAztecOO.hpp>
//#include <lifev/core/algorithm/LinearSolver.hpp>

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


Real scalardiffusion11 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 2; //  std::sin (2* pi / y ) * std::cos ( 2*pi / x ) * std::exp ( z ) ;
}

Real scalardiffusion12 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 0.0;
}

Real scalardiffusion13 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 0.0;
}
Real scalardiffusion21   ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 0.0;
}
Real scalardiffusion22   ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 1.; // std::sin (2* pi / y ) * std::cos ( 2*pi / x ) * std::exp ( z ) ;
}
Real scalardiffusion23   ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 0.0;
}

Real scalardiffusion31   ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 0.0;
}
Real scalardiffusion32   ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 0.0 ;
}
Real scalardiffusion33   ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
  return 1.; // std::sin (2* pi / y ) * std::cos ( 2*pi / x ) * std::exp ( z );
}
// Source Rhs
Real scalarRhs ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return  1* (std::sqrt( std::pow(x,2)+std::pow(y,2)+std::pow(z,2) )<1 ) ;
}



// Component Rhs
Real scalarRhs0 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return  1* (std::sqrt( std::pow(x,2)+std::pow(y,2)+std::pow(z,2) )<1 ) ; //  x;
}
// Component Rhs
Real scalarRhs1 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return  1* (std::sqrt( std::pow(x,2)+std::pow(y,2)+std::pow(z,2) )<1 ) ; // x;

}

// Component Rhs
Real scalarRhs2 ( const Real& /*t*/, const Real& x , const Real& y, const Real& z , const ID& /*i*/)
{
    return  1* (std::sqrt( std::pow(x,2)+std::pow(y,2)+std::pow(z,2) )<1 ) ; // x;
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

class MatrixdiffusionFunctor
{
public:
    typedef MatrixSmall<3,3>  return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
         MatrixSmall<3,3> Perm;
         Perm[0][0]= scalardiffusion11 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[0][1]= scalardiffusion12 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[0][2]= scalardiffusion13 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[1][0]= scalardiffusion21 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[1][1]= scalardiffusion22 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[1][2]= scalardiffusion23 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[2][0]= scalardiffusion31 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[2][1]= scalardiffusion32 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );
         Perm[2][2]= scalardiffusion33 (0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  );

       return Perm;
   }

    MatrixdiffusionFunctor() {}
    MatrixdiffusionFunctor (const MatrixdiffusionFunctor&) {}
     ~MatrixdiffusionFunctor() {}
};


class rhsScalarFunctor
{
public:
    typedef Real return_Type;

    return_Type operator() ( const VectorSmall<3> spaceCoordinates )
    {
        return scalarRhs( 0, spaceCoordinates[0], spaceCoordinates[1], spaceCoordinates[2] , 0  ) ;
    }

    rhsScalarFunctor() {}
    rhsScalarFunctor (const rhsScalarFunctor&) {}
    ~rhsScalarFunctor() {}
};



/* LaplacianNeumann1 */
class Neumann1VectorialFunctor
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

    Neumann1VectorialFunctor() {}
    Neumann1VectorialFunctor (const Neumann1VectorialFunctor&) {}
    ~Neumann1VectorialFunctor() {}
};
/* LaplacianNeumann2 */
class Neumann2VectorialFunctor
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

    Neumann2VectorialFunctor() {}
    Neumann2VectorialFunctor (const Neumann2VectorialFunctor&) {}
    ~Neumann2VectorialFunctor() {}
};





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

static bool regIF = (PRECFactory::instance().registerProduct ( "Ifpack", &createIfpack ) );
static bool regML = (PRECFactory::instance().registerProduct ( "ML", &createML ) );



typedef RegionMesh<LinearTetra> mesh_Type;
typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

typedef MatrixEpetraStructured<Real> matrix_block_type;
typedef VectorBlockMonolithicEpetra vector_block_type;
typedef MatrixEpetra<Real> matrix_type;
typedef VectorEpetra vector_type;
typedef LifeV::Preconditioner             basePrec_Type;
typedef boost::shared_ptr<basePrec_Type>  basePrecPtr_Type;

typedef LifeV::PreconditionerIfpack       prec_Type;
typedef boost::shared_ptr<prec_Type>      precPtr_Type;


int main ( int argc, char** argv )
{

#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_MpiComm (MPI_COMM_WORLD) );
#else
    boost::shared_ptr<Epetra_Comm> Comm (new Epetra_SerialComm);
#endif

    // a flag to see who's the leader for output purposes
    bool verbose = (Comm->MyPID() == 0);

        LifeChronoManager<> chronoMgr ( Comm );

        LifeChrono initTime;
        chronoMgr.add ( "Initialization Time", &initTime );
        initTime.start();

    // Open and read the data file
    GetPot command_line (argc, argv);
    string data_file_name = command_line.follow ("data", 2, "-f", "--file");
    //GetPot dataFile( data_file_name );
    GetPot dataFile ( "data" );

     initTime.stop();

        // Build and partition the mesh

        LifeChrono meshTime;
        chronoMgr.add ( "Mesh reading/creation Time", &initTime );
        meshTime.start();

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


       if ( verbose )
        {
            std::cout << " done ! " << std::endl;
        }
        if ( verbose )
        {
            std::cout << "mesh elements = " << fullMeshPtr->numElements() << "\n"
                      << "mesh points   = " << fullMeshPtr->numPoints() << std::endl;
        }

        meshTime.stop();

        if ( verbose )
        {
            std::cout << " -- Partitioning the mesh ... " << std::flush;
        }

        LifeChrono partTime;
        chronoMgr.add ( "Partition Time", &partTime );
        partTime.start();
        boost::shared_ptr< mesh_Type > meshPtr;
    {
        MeshPartitioner< mesh_Type >   meshPart (fullMeshPtr, Comm);
        meshPtr = meshPart.meshPartition();
    }

   partTime.stop();

        Int localMeshNum[ 2 ];
        localMeshNum[ 0 ] = meshPtr->numElements();
        localMeshNum[ 1 ] = meshPtr->numPoints();
        Int maxMeshNum[ 2 ] = { 0, 0 };
        Comm->MaxAll ( localMeshNum, maxMeshNum, 2 );

        if ( verbose )
        {
            std::cout << "part mesh elements = " << maxMeshNum[ 0 ] << "\n"
                      << "part mesh points   = " << maxMeshNum[ 1 ] << std::endl;
        }

        if ( verbose )
        {
            std::cout << " -- Freeing the global mesh ... " << std::flush;
        }

        fullMeshPtr.reset();
#ifdef HAVE_LIFEV_DEBUG
        ASSERT ( fullMeshPtr.use_count() == 0, "full mesh not properly freed." );
#endif
        if ( verbose )
        {
            std::cout << " done ! " << std::endl;
        }

        // Build the FESpaces

        if ( verbose )
        {
            std::cout << " -- Building FESpaces ... " << std::flush;
        }
        LifeChrono feSpaceTime;
        chronoMgr.add ( "FESpace creation Time", &feSpaceTime );
        feSpaceTime.start();

        std::string uOrder ("P2");
        std::string pOrder ("P1");

    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > uFESpace ( new FESpace< mesh_Type, MapEpetra > (meshPtr, uOrder, 3, Comm) );
    boost::shared_ptr<FESpace< mesh_Type, MapEpetra > > pFESpace ( new FESpace< mesh_Type, MapEpetra > (meshPtr, pOrder, 1, Comm) );

    if (verbose)
    {
        std::cout << std::endl << " ### Dof Summary ###: " <<  std::endl;
    }
    if (verbose)
    {
        std::cout << " Velocity  : " << uFESpace->map().map (Unique)->NumGlobalElements() << std::endl;
    }
    if (verbose)
    {
        std::cout << " Pressure  : " << pFESpace->map().map (Unique)->NumGlobalElements() << std::endl;
    }

    if (verbose)
    {
        std::cout << " Building EA FESpaces  " << std::endl;
    }


    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 3 > > ETuFESpace ( new ETFESpace< mesh_Type, MapEpetra, 3, 3 > (meshPtr, & (uFESpace->refFE() ), Comm) );
    boost::shared_ptr<ETFESpace< mesh_Type, MapEpetra, 3, 1 > > ETpFESpace ( new ETFESpace< mesh_Type, MapEpetra, 3, 1 > (meshPtr, & (pFESpace->refFE() ), Comm) );

       feSpaceTime.stop();

        if ( verbose )
        {
            std::cout << " done ! " << std::endl;
        }
    if (verbose)
    {
        std::cout << " Create vector " << std::endl;
    }


    vector_block_type DarcySolution (ETuFESpace->map() | ETpFESpace->map(), Unique);
    DarcySolution *= 0.0;
    vector_block_type Rhs (ETuFESpace->map() | ETpFESpace->map(), Repeated);
    Rhs *= 0.0;


    if (verbose)
    {
        std::cout << " Building the solvers " << std::endl;
    }
//    LinearSolver linearSolver;
    SolverAztecOO linearSolver;
    linearSolver.setCommunicator ( Comm );
    linearSolver.setDataFromGetPot (dataFile, "solver");
    linearSolver.setupPreconditioner (dataFile, "prec");


    if ( verbose )
    {
        std::cout << "done" << std::endl;
    }

    if (verbose)
    {
        std::cout << " Building the exporter " << std::endl;
    }


    ExporterHDF5<mesh_Type> exporter ( dataFile, meshPtr, "solution", Comm->MyPID() );
    exporter.setMultimesh (false);

    boost::shared_ptr<vector_type> DarcyExported ( new vector_type (DarcySolution, Repeated) );


    const UInt PressureOffset ( 3 * uFESpace->dof().numTotalDof() );

    exporter.addVariable ( ExporterData<mesh_Type>::VectorField, "velocity", uFESpace, DarcyExported, UInt (0) );
    exporter.addVariable ( ExporterData<mesh_Type>::ScalarField, "pressure", pFESpace, DarcyExported, PressureOffset);

        if ( verbose )
        {
            std::cout << " -- Defining and filling the matrix ... " << std::flush;
        }

        LifeChrono matTime;
        chronoMgr.add ( "Matrix creation Time", &matTime );
        matTime.start();

        boost::shared_ptr<matrix_block_type> DarcyMatrix (new matrix_block_type ( ETuFESpace->map() | ETpFESpace->map() ) );
        *DarcyMatrix *= 0.0;

       matTime.stop();
       QuadratureBoundary myBDQR (buildTetraBDQR (quadRuleTria4pt) );


        boost::shared_ptr<MatrixdiffusionFunctor> MatrixdiffusionFct ( new MatrixdiffusionFunctor );


        LifeChrono assemblyTime;
        chronoMgr.add ( "Assembly Time", &assemblyTime );
        assemblyTime.start();
        {
            using namespace ExpressionAssembly;

            integrate ( elements (ETuFESpace->mesh() ),
                        quadRuleTetra4pt,
                        ETuFESpace,
                        ETuFESpace,
          dot( eval(MatrixdiffusionFct, X) *  phi_j , phi_i )
                      )
                    >> *DarcyMatrix->block (0, 0);

            integrate ( elements (ETuFESpace->mesh() ),
                        quadRuleTetra4pt,
                        ETuFESpace,
                        ETpFESpace,
                       -1.* phi_j * div (phi_i)

                      )
                    >> *DarcyMatrix->block (0, 1);

            integrate ( elements ( ETuFESpace->mesh() ),
                        quadRuleTetra4pt,
                        ETpFESpace,
                        ETuFESpace,
                         -1.* phi_i * div (phi_j)
                      )
                    >> *DarcyMatrix->block (1, 0);
        }
        DarcyMatrix->globalAssemble();
        assemblyTime.stop();

        if ( verbose )
        {
            std::cout << " done! " << std::endl;
        }

 // Assembly RHS

// Neumann BC on velocity
   boost::shared_ptr<Neumann1VectorialFunctor> VectorialFctN1 ( new Neumann1VectorialFunctor );
   boost::shared_ptr<Neumann2VectorialFunctor> VectorialFctN2 ( new Neumann2VectorialFunctor );
// source 2nd term Rhs
   boost::shared_ptr<rhsScalarFunctor> scalarFctRhs ( new rhsScalarFunctor );


   vector_block_type DarcyRhs ( ETuFESpace->map() | ETpFESpace->map() , Repeated );
   DarcyRhs *= 0.0;

    {
            using namespace ExpressionAssembly;
        integrate ( boundary (ETuFESpace->mesh(),BCFlags::TOPWALL), // Mesh

                    myBDQR, // QR

                    ETuFESpace,

                     dot(  eval ( VectorialFctN1 , X ) , phi_i)

                  )
                  >>*DarcyRhs.block(0);

     integrate ( boundary (ETuFESpace->mesh(),BCFlags::FRONTWALL), // Mesh

                    myBDQR, // QR

                    ETuFESpace,

                     dot(  eval ( VectorialFctN2 , X ) , phi_i)

                  )
                  >>*DarcyRhs.block(0);

     integrate ( elements (ETpFESpace->mesh()), // Mesh
                   pFESpace->qr(), // QR

                    ETpFESpace,

                    eval ( scalarFctRhs, X ) * phi_i

                  )
                 >>*DarcyRhs.block(1);


    }

    if (verbose)
    {
        std::cout << " -- Closing the matrix ... " << std::flush;
    }
     DarcyRhs.globalAssemble();
     vector_block_type DarcyRhsUnique ( DarcyRhs, Unique );


LifeChrono applyBC;
        chronoMgr.add ( "Apply BC", &applyBC );
        applyBC.start();
   BCHandler bcHandler;
   BCFunctionBase dirichletBCFct ( dirichlet );
    bcHandler.addBC ("Left", BCFlags::LEFTWALL, Essential, Full , dirichletBCFct, 3);
    bcHandler.addBC ("Right", BCFlags::RIGHTWALL, Essential, Full , dirichletBCFct, 3);
    bcHandler.addBC ("Back", BCFlags::BACKWALL, Essential, Full , dirichletBCFct, 3);
    bcHandler.bcUpdate ( *meshPtr, uFESpace->feBd(), uFESpace->dof() );


        applyBC.stop();


        if (verbose)
        {
            std::cout << "[Darcy] Solving the system " << std::endl;
        }
        LifeChrono solverTime;
        chronoMgr.add ( "Solver Time", &solverTime );
        solverTime.start();

    MapEpetra fullMap(ETuFESpace->map()+ETpFESpace->map());

        boost::shared_ptr<matrix_type> DarcyMatrixNoBlock (new matrix_type (
        fullMap));
        *DarcyMatrixNoBlock +=*DarcyMatrix;
         DarcyMatrixNoBlock->globalAssemble();

         bcManage (*DarcyMatrixNoBlock, DarcyRhsUnique,
              *uFESpace->mesh(), uFESpace->dof(),
              bcHandler, uFESpace->feBd(), 1.0, Real (0.0) );


        linearSolver.setMatrix(*DarcyMatrixNoBlock);
        linearSolver.solveSystem(DarcyRhsUnique,DarcySolution,DarcyMatrixNoBlock);

        solverTime.stop();

    exporter.postProcess (1.0);
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

       // print out times
        chronoMgr.print();


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
