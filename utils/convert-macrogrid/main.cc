#include <config.h>

#include <time.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <list>
#include <vector>

#include <dune/alugrid/impl/serial/serialize.h>


// ElementRawID
// ------------

enum ElementRawID { TETRA_RAW = 4, HEXA_RAW = 8 };


// Vertex
// ------

struct Vertex
{
  int id;
  double x, y, z;
};


// Element
// -------

template< ElementRawID >
struct Element;

template<>
struct Element< TETRA_RAW >
{
  static const int numVertices = 4;
  int vertices[ numVertices ];
};

template<>
struct Element< HEXA_RAW >
{
  static const int numVertices = 8;
  int vertices[ numVertices ];
};


// BndSeg
// ------

template< ElementRawID >
struct BndSeg;

template<>
struct BndSeg< TETRA_RAW >
{
  static const int numVertices = 3;
  int bndid;
  int vertices[ numVertices ];
};

template<>
struct BndSeg< HEXA_RAW >
{
  static const int numVertices = 4;
  int bndid;
  int vertices[ numVertices ];
};


// Periodic
// --------

template< ElementRawID >
struct Periodic;

template<>
struct Periodic< TETRA_RAW >
{
  static const int numVertices = 6;
  int bndid;
  int vertices[ numVertices ];
};

template<>
struct Periodic< HEXA_RAW >
{
  static const int numVertices = 8;
  int bndid;
  int vertices[ numVertices ];
};


// readLegacyFormat
// ----------------

template< ElementRawID rawId >
void readLegacyFormat ( std::istream &input,
                        std::vector< Vertex > &vertices,
                        std::vector< Element< rawId > > &elements,
                        std::vector< BndSeg< rawId > > &bndSegs,
                        std::vector< Periodic< rawId > > &periodics )
{
  // Das alte Format sieht im wesentlichen so aus:
  //
  // <Anzahl der Knoten : int >     /* 1.Zeile der Datei
  // <x-Koordinate : float>  <y-Koo. : float>  <z-Koo. : float>
  // ...            /* f"ur den letzten Knoten
  // <Anzahl der Elemente : int>
  // <KnotenNr. 0: int> ... <KnotenNr. 7: int>  /* f"ur das erste Hexaederelement
  // ...            /* f"ur das letzte Hexaederelement
  // <Anzahl der Randfl"achen : int>
  // <Randtyp>  4  <KnotenNr. 0> ... <KnotenNr. 3>/* erste Randfl"ache
  // ...            /* letzte Randfl"ache
  // <Identifier f"ur den 0. Knoten : int>  /* Identifierliste ist im seriellen
  // ...            /* Verfahren oder beim Aufsetzen aus
  // <Identifier f"ur den letzten Knoten : int> /* einem Gitter optional, sonst muss
  //            /* jeder Vertex eine eigene Nummer haben

  std::cout << "Reading legacy format..." << std::endl;

  int nv = 0;
  input >> nv;
  vertices.resize( nv );
  for( int i = 0; i < nv; ++i )
    input >> vertices[ i ].x >> vertices[ i ].y >> vertices[ i ].z;
  std::cout << "  - read " << nv << " vertices." << std::endl;

  int ne = 0;
  input >> ne;
  elements.resize( ne );
  for( int i = 0; i < ne; ++i )
  {
    for( int vx = 0; vx < Element< rawId >::numVertices; ++vx )
      input >> elements[ i ].vertices[ vx ];
  }
  std::cout << "  - read " << ne << " elements." << std::endl;

  int temp_nb = 0;
  input >> temp_nb;
  bndSegs.reserve( temp_nb );
  periodics.reserve( temp_nb );
  for( int i = 0; i < temp_nb; ++i )
  {
    int n, bndid;
    input >> bndid >> n;
    if( n == BndSeg< rawId >::numVertices )
    {
      BndSeg< rawId > seg;
      seg.bndid = bndid;
      for( int vx = 0; vx < BndSeg< rawId >::numVertices; ++vx )
        input >> seg.vertices[ vx ];
      bndSegs.push_back( seg );
    }
    else if( n == Periodic< rawId >::numVertices )
    {
      Periodic< rawId > seg;
      seg.bndid = bndid;
      for( int vx = 0; vx < Periodic< rawId >::numVertices; ++vx )
        input >> seg.vertices[ vx ];
      periodics.push_back( seg );
    }
    else
    {
      std::cerr << "ERROR (fatal): Invalid number of vertices for boundary object: " << n << "." << std::endl;
      std::exit( 1 );
    }
  }
  std::cout << "  - read " << bndSegs.size() << " boundary segments." << std::endl;
  std::cout << "  - read " << periodics.size() << " periodic boundary segments." << std::endl;

  if( !input )
  {
    std::cerr << "ERROR (fatal): Unexpected end of file." << std::endl;
    std::exit( 1 );
  }

  for( int i = 0; i < nv; ++i )
  {
    int dummy;
    input >> vertices[ i ].id >> dummy;
  }

  if( !input )
  {
    std::cerr << "WARNING (ignored) No parallel identification applied due to incomplete (or non-existent) identifier list." << std::endl;
    for( int i = 0; i < nv; ++ i )
      vertices[ i ].id = i;
  }
}



// readMacroGrid
// -------------

template< class stream_t, ElementRawID rawId >
void readMacroGrid ( stream_t &input,
                     std::vector< Vertex > &vertices,
                     std::vector< Element< rawId > > &elements,
                     std::vector< BndSeg< rawId > > &bndSegs,
                     std::vector< Periodic< rawId > > &periodics )
{
  // ignore firstline for now; it should have been stripped already
  std::string firstline;
  ALUGrid::getline( input, firstline );

  int vertexListSize = 0;
  input >> vertexListSize;
  vertices.resize( vertexListSize );
  for( int i = 0; i < vertexListSize; ++i )
    input >> vertices[ i ].id >> vertices[ i ].x >> vertices[ i ].y >> vertices[ i ].z;

  int elementListSize = 0, type = -1;
  input >> elementListSize >> type;
  if( type != rawId )
  {
    std::cerr << "ERROR (fatal): Invalid number of vertices for element (got " << type << ", expected " << rawId << ")." << std::endl;
    std::exit( 1 );
  }
  elements.resize( elementListSize );
  for( int i = 0; i < elementListSize; ++i )
  {
    for( int j = 0; j < Element< rawId >::numVertices; ++j )
      input >> elements[ i ].vertices[ j ];
  }

  int bndSegListSize = 0;
  int periodicListSize = 0;
  input >> bndSegListSize >> periodicListSize;
  periodics.resize( periodicListSize );
  for( int i = 0; i < periodicListSize; ++i )
  {
    for( int j = 0; j < Periodic< rawId >::numVertices; ++j )
      input >> periodics[ i ].vertices[ j ];
  }
  bndSegs.resize( bndSegListSize );
  for( int i = 0; i < bndSegListSize; ++i )
  {
    input >> bndSegs[ i ].bndid;
    for( int j = 0; j < BndSeg< rawId >::numVertices; ++j )
      input >> bndSegs[ i ].vertices[ j ];
  }
}



// writeMacroGrid
// --------------

template< class stream_t, ElementRawID rawId >
void writeMacroGrid ( stream_t &output,
                      const std::vector< Vertex > &vertices,
                      const std::vector< Element< rawId > > &elements,
                      const std::vector< BndSeg< rawId > > &bndSegs,
                      const std::vector< Periodic< rawId > > &periodics )
{
  const int vertexListSize = vertices.size();
  const int elementListSize = elements.size();
  const int bndSegListSize = bndSegs.size();
  const int periodicListSize = periodics.size();

  ALUGrid::StandardWhiteSpace_t ws;

  // write header line as vector of characters 
  std::ostringstream firstline;
  firstline << (rawId == HEXA_RAW ? "!Hexahedra" : "!Tetrahedra");
  firstline << "  ( noVertices = " << vertexListSize << " | noElements = " << elementListSize << " )" << std::endl;
  output << firstline.str();

  output << std::endl << vertexListSize << std::endl;
  for( int i = 0; i < vertexListSize; ++i )
    output << vertices[ i ].id << ws << vertices[ i ].x << ws << vertices[ i ].y << ws << vertices[ i ].z << std::endl;

  output << std::endl << elementListSize << ws << int( rawId ) << std::endl;
  for( int i = 0; i < elementListSize; ++i )
  {
    output << elements[ i ].vertices[ 0 ];
    for( int j = 1; j < Element< rawId >::numVertices; ++j )
      output << ws << elements[ i ].vertices[ j ];
    output << std::endl;
  }

  output << std::endl << periodicListSize << ws << bndSegListSize << std::endl;
  for( int i = 0; i < periodicListSize; ++i )
  {
    output << periodics[ i ].vertices[ 0 ];
    for( int j = 1; j < Periodic< rawId >::numVertices; ++j )
      output << ws << periodics[ i ].vertices[ j ];
    output << std::endl;
  }
  for( int i = 0; i < bndSegListSize; ++i )
  {
    output << bndSegs[ i ].bndid;
    for( int j = 0; j < BndSeg< rawId >::numVertices; ++j )
      output << ws << bndSegs[ i ].vertices[ j ];
    output << std::endl;
  }
}


// convertLegacyFormat
// -------------------

template< ElementRawID rawId, class stream_t >
void convertLegacyFormat ( std::istream &input, stream_t &output )
{
  const clock_t start = clock();

  std::vector< Vertex > vertices;
  std::vector< Element< rawId > > elements;
  std::vector< BndSeg< rawId > > bndSegs;
  std::vector< Periodic< rawId > > periodics;

  readLegacyFormat( input, vertices, elements, bndSegs, periodics );
  writeMacroGrid( output, vertices, elements, bndSegs, periodics );

  std::cout << "INFO: Conversion of legacy macro grid format used " << (double( clock () - start ) / double( CLOCKS_PER_SEC )) << " s." << std::endl;
}

template< class stream_t >
void convertLegacyFormat ( std::istream &input, stream_t &output )
{
  std::string firstline;
  std::getline( input, firstline );
  if( firstline[ 0 ] == char( '!' ) )
  {
    if( (firstline.find( "Tetrahedra" ) != firstline.npos) || (firstline.find( "Tetraeder" ) != firstline.npos) )
      convertLegacyFormat< TETRA_RAW >( input, output );
    else if( (firstline.find( "Hexahedra" ) != firstline.npos) || (firstline.find( "Hexaeder" ) != firstline.npos) )
      convertLegacyFormat< HEXA_RAW >( input, output );
    else 
    {
      std::cerr << "ERROR: Unknown comment to file format (" << firstline << ")." << std::endl;
      std::exit( 1 );
    }
  }
  else
  {
    std::cerr << "WARNING: No identifier for file format found. Trying to read as hexahedral grid." << std::endl;
    convertLegacyFormat< HEXA_RAW >( input, output );
  }
}


// main
// ----

int main ( int argc, char **argv )
{
  bool writeBinary = false;

  for( int i = 1; i < argc; ++i )
  {
    if( argv[ i ][ 0 ] != '-' )
      continue;

    for( int j = 1; argv[ i ][ j ]; ++j )
    {
      if( argv[ i ][ j ] == 'b' )
        writeBinary = true;
    }

    std::copy( argv + (i+1), argv + argc, argv + i );
    --i; --argc;
  }

  if( argc <= 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " [-b] <input> <output>" << std::endl;
    std::cerr << "Flags: -b : write binary output" << std::endl;
    return 1;
  }

  std::ifstream input( argv[ 1 ] );
  std::ofstream output( argv[ 2 ] );
  if( writeBinary )
  {
    ALUGrid::ObjectStream os;
    convertLegacyFormat( input, os );
    output.write( os.getBuff( 0 ), os.size() );
  }
  else
    convertLegacyFormat( input, output );
}
