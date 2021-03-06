#define DISABLE_DEPRECATED_METHOD_CHECK 1

#include <config.h>

#include <iostream>
#include <sstream>
#include <string>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/alugrid/dgf.hh>
#include <dune/alugrid/common/fromtogridfactory.hh>

int main (int argc , char **argv)
{
  // this method calls MPI_Init, if MPI is enabled
  Dune::MPIHelper::instance( argc, argv );

  try {
    typedef Dune::ALUGrid< 3, 3, Dune::cube, Dune::nonconforming > GridType;
    std::string filename( "./dgf/simplex-testgrid-3-3.dgf" );
    std::cout << "READING from " << filename << std::endl;
    Dune::GridPtr< GridType > gridPtr(filename);

    Dune::FromToGridFactory< GridType > factory;
    std::shared_ptr< GridType > newPtr = factory.convert( *gridPtr );

    Dune::VTKWriter< GridType::LeafGridView > writer( newPtr->leafGridView() );
    writer.write( "dump.vtu" );

  }
  catch( Dune::Exception &e )
  {
    std::cerr << e << std::endl;
    return 1;
  }
  catch( ... )
  {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;
}
