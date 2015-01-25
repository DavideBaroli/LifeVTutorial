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
    @file basic_test.cpp
    @brief Test the consistency of the mesh data structure

    @author Alessio Fumagalli <alessio.fumagalli@polimi.it>
    @contributor Davide Baroli <davide.baroli@polimi.it>
    @maintainer

    @date 2012-09-14

    Colour a mesh with two different colours, count the number
    of elements equal of one of the two.
    Switch On and Off the Macro we can extend the example in 1D,2D and 3D.

 */

// ===================================================
//! Includes
// ===================================================

#include <Epetra_ConfigDefs.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <lifev/core/mesh/RegionMesh3DStructured.hpp>
#include <lifev/core/mesh/MeshUtility.hpp>
#include <lifev/core/mesh/MeshData.hpp>

#include <lifev/core/array/VectorSmall.hpp>
#ifdef HAVE_HDF5
#include <lifev/core/filter/ExporterHDF5.hpp>
#else
#include <lifev/core/filter/ExporterEmpty.hpp>
#endif

//#define DIM dataFile("mesh/Dimension",1);


using namespace LifeV;

UInt colour_fun ( const VectorSmall<3>& barycentre )
{
    if ( barycentre[0] < 0.5 && barycentre[1] < 0.5 )
    {
        return 2;
    }
    return 3;
}







int main (int argc, char* argv[])
{
// Read From dataFile
   GetPot command_line (argc, argv);
    const std::string data_file_name =
    command_line.follow ("data", 2, "-f", "--file");
    GetPot dataFile (data_file_name);

   const std::string discretization_section = "mesh";



#ifdef HAVE_MPI
    MPI_Init (&argc, &argv);
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_MpiComm ( MPI_COMM_WORLD ) );
#else
    boost::shared_ptr<Epetra_Comm> Comm ( new Epetra_SerialComm );
#endif

    // Create the mesh.
    MeshData meshData (dataFile, (discretization_section).c_str() );

 // set the Dimension.

    typedef RegionMesh<LinearTetra> mesh_Type;
    typedef boost::shared_ptr<mesh_Type> meshPtr_Type;

    meshPtr_Type fullMeshPtr (new mesh_Type(Comm) );
   // Select if the mesh is structured or not
    if ( meshData.meshType() != "structured" )
    {
        // Set up the mesh
        readMesh ( *fullMeshPtr, meshData );

    }
    else
    {
      regularMesh3D (*fullMeshPtr, 0,
                   dataFile ( "mesh/nx", 20 ),
                   dataFile ( "mesh/ny", 20 ),
                   dataFile ( "mesh/nz", 20 ),
                    dataFile ( "mesh/verbose", false ),
                   dataFile ( "mesh/lx", 1. ),
                   dataFile ( "mesh/ly", 1. ),
                   dataFile ( "mesh/lz", 1. ),
                   dataFile ( "mesh/tx", 0. ),
                   dataFile ( "mesh/ty", 0. ),
                   dataFile ( "mesh/tz", 0. )
        );

   }





    // Colour the mesh according to a function.
    MeshUtility::assignRegionMarkerID ( *fullMeshPtr, colour_fun );

    // Count the number of elements with colour 2
    const UInt colourElements = fullMeshPtr->elementList().countElementsWithMarkerID ( 2, std::equal_to<markerID_Type>() );

    // Number of elements with colour 2
    const UInt exactNumber = 44;

    {
        // Needed to correctly destroy the exporterHDF5

        // Set the exporter for the mesh region.
    std::string nameoutput;

    if ( meshData.meshType() != "structured" )
    {
       nameoutput="test_External";
    }
    else{


 switch(mesh_Type::S_geoDimensions){
 case 1:
    nameoutput="test_line";
 case 2:
    nameoutput="test_Tria";
 case 3:
    nameoutput="test_Tetra";

   }
  }


#ifdef HAVE_HDF5
        boost::shared_ptr<ExporterHDF5< mesh_Type> > exporter;
        exporter.reset(new ExporterHDF5<mesh_Type>( dataFile,nameoutput));

#else
       boost::shared_ptr<ExporterEmpty< mesh_Type> > exporter;
        exporter.reset(new ExporterEmpty<mesh_Type>( dataFile,nameoutput));
#endif
        exporter->setPostDir("./");
        // Set the mesh.
        exporter->setMeshProcId ( fullMeshPtr, Comm->MyPID() );

        // Export the region marker ID.
        exporter->exportRegionMarkerID ( fullMeshPtr, Comm );

        // Do the export.
        exporter->postProcess ( 0 );

    } // Needed to correctly destroy the exporterHDF5

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    if ( colourElements == exactNumber )
    {
        return ( EXIT_SUCCESS );
    }
    else
    {
        return ( EXIT_FAILURE );
    }
}
