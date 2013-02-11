/** include config file generated by configure
 *  (i.e., know what grids are present, etc)
 *  this should always be included first */
#include <config.h>
/** standard headers **/
#include <iostream>
/** dune (mpi, field-vector and grid type for dgf) **/
#include <dune/common/mpihelper.hh>     
#include <dune/common/fvector.hh>        
#include <dune/common/timer.hh>        

typedef Dune::GridSelector::GridType Grid;

/** numerical scheme **/
#include "../piecewisefunction.hh"

/** adaptation scheme **/
#include "adaptation.hh"

// method
// ------
void method ( int startLevel, int maxLevel, const char* outpath )
{
  /* Grid construction ... */
  std::string name = "../unitcube3d.dgf" ;
  // create grid pointer and release to free memory of GridPtr
  Grid* gridPtr = Dune::GridPtr<Grid>(name).release() ;

  Grid &grid = *gridPtr;

  LoadBalanceHandle<Grid> ldb(grid);
  typedef Dune::LoadBalanceHandleIF< LoadBalanceHandle<Grid> > DataHandleInterface;
  // grid.loadBalance( (DataHandleInterface&)(ldb) );
  grid.loadBalance();
  const bool verboseRank = grid.comm().rank() == 0 ;

  std::string outPath( outpath );
  const bool writeOutput = ( outPath != "none" ) ;

  /* ... some global refinement steps */
  if( verboseRank ) 
    std::cout << "globalRefine: " << startLevel << std::endl;
  grid.globalRefine( startLevel );

  /* get view to leaf grid */
  typedef Grid::Partition< Dune::Interior_Partition >::LeafGridView GridView;
  GridView gridView = grid.leafView< Dune::Interior_Partition >();

  /* construct data vector for solution */
  typedef PiecewiseFunction< GridView, Dune::FieldVector< double, 1 > > DataType;
  DataType solution( gridView );

  /* create VTK writer for data sequqnce */
  Dune::VTKSequenceWriter< GridView > vtkOut( gridView, "solution", outPath, ".", Dune::VTK::nonconforming );
  if( writeOutput ) 
  {
    VTKData< DataType >::addTo( solution, vtkOut );
    VTKData< DataType >::addPartitioningData( grid.comm().rank(), vtkOut );
  }

  /* create adaptation method */
  typedef LeafAdaptation< Grid > AdaptationType;
  AdaptationType adaptation( grid );

  if( writeOutput ) 
  {
    /* output the initial grid and the solution */
    vtkOut.write( 0.0 );
  }


  /* final time for simulation */
  const double endTime = 1.;
  /* interval for saving data */
  const double saveInterval = 0.02;
  /* first point where data is saved */
  double saveStep = saveInterval;

  /* now do the time stepping */
  unsigned int step = 0;
  double time = 0.0;
  while ( time < endTime ) 
  {
    double dt = saveInterval;

    /* augment time */
    time += dt;
    ++step;

    /* check if data should be written */
    if( time >= saveStep )
    {
      if( writeOutput ) 
      {
        /* visualize with VTK */
        vtkOut.write( time );
      }
      /* set saveStep for next save point */
      saveStep += saveInterval;
    }

    adaptation( solution );

  }           

  if( writeOutput ) 
  {
    /* output final result */
    vtkOut.write( time );
  }

  // delete grid 
  delete gridPtr ;
}
/***************************************************
 ** main program with parameters:                 **
 ** 1) number of problem to use (initial data...) **
 ** 2) number of global refinement steps          **
 ** 3) maximal level to use during refinement     **
 ***************************************************/
int main ( int argc , char **argv )
try
{
  /* initialize MPI, finalize is done automatically on exit */
  Dune::MPIHelper &mpi = Dune::MPIHelper::instance( argc, argv );
  
  if( argc < 1 )
  {
    /* display usage */
    if( mpi.rank() == 0 )
      std::cout << "Usage: " << argv[ 0 ] << " [startLevel] [maxLevel]" << std::endl;
    return 0;
  }

  /* get level to use for computationa */
  const int startLevel = (argc > 1 ? atoi( argv[ 1 ] ) : 0);
  const int maxLevel = (argc > 2 ? atoi( argv[ 2 ] ) : startLevel);

  const char* path = (argc > 3) ? argv[ 3 ] : "./";
  method( startLevel, maxLevel, path );

  /* done */
  return 0;
}
catch( const std::exception &e )
{
  std::cout << "STL ERROR: " << e.what() << std::endl;
  return 1;
}
catch( const Dune::Exception &e )
{
  std::cout << "DUNE ERROR: " << e << std::endl;
  return 1;
}
catch( ... )
{
  std::cout << "Unknown ERROR" << std::endl;
  return 1;
}
