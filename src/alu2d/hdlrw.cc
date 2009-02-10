// typdef these stream because this code uses a lot strstream
typedef basic_stringbuf<char>  strstreambuf_t ;

#include "grid.h"
#include "triang.h"
#include "handle.h"

bool Hmesh :: ascireadtriang(const char *filename,
           double& time, unsigned long int& nbr) 
{
#ifndef NDEBUG 
  cerr << "\n  Hmesh_basic::ascireadtriang(?) opens: " ;
  cerr << filename << "\n" << endl ;
#endif
  
  ifstream in;
  in.open(filename, ios::in) ;

  if (!in.good()) {
    in.clear();
    string macro(filename);
    macro+=".macro";
    cerr << "Warning: file " << filename << " not found, trying " << macro << endl;
    in.open(macro.c_str(), ios::in) ;
  }
  assert(in) ;
  return ascireadtriang(in,time,nbr);
}
bool Hmesh :: ascireadtriang(ifstream &in,
           double& time, unsigned long int& nbr)

{
  bool isbackup=false;
  // Wiederaufsetzen?
  int c = in.get () ;
  assert (!in.eof ()) ;
  in.putback (c) ;
  assert (in.good ()) ;
  if (c == int ('!')) {
    // Kommentar gefunden: Die erste Zeile in den strstreambuf buf lesen
    // und auf 'Backup' untersuchen.
    strstreambuf_t buf ;
    in.get () ;         // Das Kommentarzeichen wird entfernt.
    in.get (buf) ;
    int len = in.gcount () ;
    // in.get () ;         // Der folgende Zeilenumbruch wird auch entfernt.
    istream is (& buf) ;
    char * str = new char [len + 1] ;
    assert (str) ;
    is >> str ;         // Das erste Wort nach dem Kommentar steht jetzt in str.
    // Alle weiteren k"onnen noch aus is gelesen werden, das
    // array str ist so lang, wie die gesamte Zeile in 'buf'.
    
    if (0 == strcmp (str, "Backup")) {
      cerr << "Backup-file found" << endl;
      // cerr << "STIRNG: " << buf << endl;
      int rrule;
      is >> time >> nbr;
      is >> _nconfDeg >> rrule; 
      refinement_rule=(Refco::tag_t)rrule;
      assert(_nconfDeg>=0);
      assert(_nconfDeg==0 || refinement_rule==Refco::quart);
      isbackup=true;
    } 
    else if (0 != strcmp (str, "Triangles")) 
    {
      cerr << "Error in Hmesh :: ascireadtriang: "
           << "Wrong macrogrid format: " << endl;
      cerr << "file with !-line but command not recognized" << endl;
      abort();
    }
    delete [] str;
  }
  Hmesh_basic::ascireadtriang(in);
  return isbackup;
}
 
void Hmesh_basic :: ascireadtriang(ifstream &in) {

  Vertex ** v = 0 ;

  int nv = 0 ;

  {
    in >> nv ;
    
#ifndef NDEBUG
    cerr << "    Number of Vertices:           " << nv << endl ;
#endif
    
    v = new Vertex *[nv] ;

    assert(v) ;
    
    for(int i = 0; i < nv ; i ++) {
         
      Vertex * n = new Fullvertex() ;
    
      n->read(in) ;

      v[i] = n ;

      vl.insert(n) ;
      
    }

  }

  {

    int ne = 0 ;
    
    in >> ne ;
  
#ifndef NDEBUG
    cerr << "    Number of MacroElements:      " << ne << endl ;
#endif
    
 
    for(int i = 0; i < ne ; i ++) {
    
      Triang & tr = * new Triang() ;

      tr.read(in, v, nv) ;

      mel.insert( new macroelement_t (tr) ) ;
    
    }

  }

  {

    int nb = 0 ;

    in >> nb ;

#ifndef NDEBUG 
    cerr << "    Number of BoundarySegments:   " << nb << endl ;
#endif

    typedef struct
      {
        int pernb;
        Bndel *b;
      } perbnd_struct;
    perbnd_struct *perbnd_list;
    int perbnd_card=0,perbnd_ok=0;
    //perbnd_list = (perbnd_struct *)malloc(sizeof(perbnd_struct)*nb);
    perbnd_list = new  perbnd_struct[nb];

    typedef struct
      {
        double p0,p1;
        Bndel *b;
      } axis_struct;
    axis_struct *x_axis,*y_axis;
    int y_card=0,x_card=0,x_ok=0,y_ok=0;

    //x_axis = (axis_struct *)malloc(sizeof(axis_struct)*nb);
    //y_axis = (axis_struct *)malloc(sizeof(axis_struct)*nb);
    x_axis = new axis_struct[nb];
    y_axis = new axis_struct[nb];

    for( int i = 0 ; i <  nb ; i ++ ) 
    {

      int lt;
      
      lt=in.peek();

      Bndel::bnd_t t=(Bndel::bnd_t)lt;

      Bndel_triang * b;

      int generalperbnd=0;

      switch (t)
      {
  case Bndel::periodic:
          b=new Bndel_periodic();
    break;
  case Bndel::general_periodic:
    t=Bndel::periodic;
          b=new Bndel_periodic();
    generalperbnd=1;
    break;
        default:
          b=new Bndel_triang(t);
    break;
      }

      b->read(in, v, nv) ;

      mbl.insert( new macrobndel_t (*b) ) ;

      if (t==Bndel::periodic)
      {
  if (generalperbnd) {
    int pernb;
    in >> pernb;
          if (pernb<i) {
      assert(perbnd_list[pernb].pernb==i);
      ((Bndel_periodic*)b)->set_pnb(perbnd_list[pernb].b);
      ((Bndel_periodic*)perbnd_list[pernb].b)->set_pnb(b);
      perbnd_ok++;
    }
    else {
      perbnd_list[i].b=b;
      perbnd_list[i].pernb=pernb;
      perbnd_card++;
    }
  }
        else if (fabs(b->vertex(0)->coord()[0]-b->vertex(1)->coord()[0])<EPS)
  {

          double y0,y1;
          if (b->vertex(0)->coord()[1]<b->vertex(1)->coord()[1])
      { y0=b->vertex(0)->coord()[1];y1=b->vertex(1)->coord()[1]; }
          else
      { y0=b->vertex(1)->coord()[1];y1=b->vertex(0)->coord()[1]; }
          int y;
          for (y=0;y<y_card;y++)
            if (fabs(y_axis[y].p0-y0)+fabs(y_axis[y].p1-y1)<EPS)
      {
              ((Bndel_periodic*)b)->set_pnb(y_axis[y].b);
              ((Bndel_periodic*)y_axis[y].b)->set_pnb(b);
              y_ok++;
              break;
            }
          if (y==y_card)
    {
            y_axis[y_card].p0=y0;
            y_axis[y_card].p1=y1;
            y_axis[y_card].b=b;
            y_card++;
    }
        }
        else 
  {
          assert(fabs(b->vertex(0)->coord()[1]-b->vertex(1)->coord()[1])<EPS);

          double x0,x1;
          if (b->vertex(0)->coord()[0]<b->vertex(1)->coord()[0])
      { x0=b->vertex(0)->coord()[0];x1=b->vertex(1)->coord()[0]; }
          else
      { x0=b->vertex(1)->coord()[0];x1=b->vertex(0)->coord()[0]; }
          int x;
          for (x=0;x<x_card;x++)
            if (fabs(x_axis[x].p0-x0)+fabs(x_axis[x].p1-x1)<EPS)
      {
              ((Bndel_periodic*)b)->set_pnb(x_axis[x].b);
              ((Bndel_periodic*)x_axis[x].b)->set_pnb(b);
              x_ok++;
              break;
            }
          if (x==x_card)
    {
            x_axis[x_card].p0=x0;
            x_axis[x_card].p1=x1;
            x_axis[x_card].b=b;
            x_card++;
    }
        }
      }

    }

    delete [] y_axis;
    delete [] x_axis;
    delete [] perbnd_list;
    //free(y_axis);
    //free(x_axis);

    assert(y_ok==y_card);
    assert(x_ok==x_card);

  }       

  delete[] v ;
  
#ifndef NDEBUG
  cerr << "\n  -------------------------- closed.\n" <<endl ;
#endif

  vl.renumber() ;

  makeneighbours() ;

  {
    Listwalk_impl < macroelement_t > walk(mel) ;
    for (walk.first() ; !walk.done() ; walk.next() )
    {
      Triang& tr=( (Triang&)(*walk.getitem()) );
      for (int l=0;l<3;l++) {
        if (!tr.normaldir(l)) 
        {
          tr.setnormdir(l,1);
          if (tr.neighbour(l)->thinis(Thinelement::element_like))
            tr.nbel(l)->setnormdir(tr.opposite(l),-1); 
        }
  if (tr.neighbour(l)->edge(tr.opposite(l))) {
    tr.edgeconnect(l,tr.neighbour(l)->edge(tr.opposite(l)));
  } else {
    Edge *e=new Edge(this);
    tr.edgeconnect(l,e);
  }
      }
    }
    for (walk.first() ; !walk.done() ; walk.next() )
    {
      walk.getitem()->setorientation();
      walk.getitem()->setrefine();
    }
  }
  {
    Listwalk_impl < macrobndel_t > walk(mbl) ;
    for (walk.first() ; !walk.done() ; walk.next() ) {
      walk.getitem()->setorientation();
      walk.getitem()->edgeconnect(0,walk.getitem()->neighbour(0)->edge(walk.getitem()->opposite(0)));
    }
  }
}

void Hmesh::asciwritetriang(const char *filename,
          double time, unsigned long int nbr) {

#ifndef NDEBUG
  cerr << "\n  Hmesh_basic::asciwritetriang(?) opens: " ;
  cerr << filename << "\n" << endl ;
#endif

  vl.renumber() ;
  
  ofstream out(filename, ios::out|ios::trunc) ;
  
  out.setf(ios::fixed, ios::floatfield) ;
  
  out.precision(14) ;

  out << "!Backup ";
  out << time << " " << nbr << " ";
  out << _nconfDeg << " " << refinement_rule << endl;
  Hmesh_basic::asciwritetriang(out);
}
 
void Hmesh_basic::asciwritetriang(ofstream &out) {
  {
 
    Listwalk_impl < Vertex > walk(vl) ;
  
#ifndef NDEBUG
    cerr << "    Number of Vertices:       " << walk.size() << endl ;
#endif
    
    int nr = 0;
    
    for( walk.first() ; ! walk.done() ; walk.next() ) {
    
      Vertex & v = walk.getitem() ;

      if (v.level()==-1) ++nr;
              
    }

    out << nr << endl;

    for( walk.first() ; ! walk.done() ; walk.next() ) {
    
      Vertex & v = walk.getitem() ;

      if (v.level()==-1) 
        v.write(out) ;
    
    }
    
  }
  
  {

    Listwalk_impl < macroelement_t > walk(mel) ;

    int count = 0 ;
    
#ifndef NDEBUG
    cerr << "    Number of macro Elements:  " << walk.size() << endl ;
#endif
    
    out << walk.size() << endl;
    
    for( walk.first() ; ! walk.done() ; walk.next() ) {
    
      walk.getitem()->write(out) ;

      count += walk.getitem()->count() ;

    }

#ifndef NDEBUG
    cerr << "    Number of Elements:       " << count << endl ;
#endif

  }
  
  {

    Listwalk_impl < macrobndel_t > walk(mbl) ;

   int count = 0 ;
    
#ifndef NDEBUG
    cerr << "    Number of macro boundary Elements:   " << walk.size() << endl ;
#endif
    
    out << walk.size() << endl ;
    
    for( walk.first() ; ! walk.done() ; walk.next() ) {
    
      walk.getitem()->write(out) ;

      count += walk.getitem()->count() ;

    }

#ifndef NDEBUG
    cerr << "    Number of boundary Elements:       " << count << endl ;
#endif

  }
  
#ifndef NDEBUG
  cerr << "\n  -------------------------- closed.\n" <<endl ;
#endif
 
}

void
Hmesh::storeGrid(const char* fbase,
     double time, unsigned long int nbr)
{
  char *filename;

  filename=new char[strlen(fbase)+15];
  // filename=(char*)alloca(strlen(fbase)+15);
  sprintf(filename,"%s",fbase);

  asciwritetriang (filename,time,nbr);

  sprintf(filename,"%s.refine",fbase);

#ifndef NDEBUG
  cerr << "Hmesh::writeRecoverFile(): writing file \""
       << filename << "\" ...";
#endif

  ofstream out(filename, ios::out|ios::trunc);
  assert(out);

  // out.setf(ios::scientific, ios::floatfield);
  // out.precision(16);

  // Zeit und Nummer des Zeitschritts schreiben
  /*
  out << time << " ";
  out << nbr << endl;
  */
  // Status des Gitters sichern
  for( int level = 0 ;; level++ ) {
    Levelwalk<Element> walk(mel, level);
    if( !walk.size() ) 
    {
      break;
    } 
    else 
    {
      for( walk.first() ; !walk.done() ; walk.next() )
      	out.put(walk.getitem().splitrule());
      // out << walk.getitem().splitrule() << " ";
    }
    // out << endl;
  }
  // out << endl;

  /*
  // Daten sichern
  Leafwalk<Element> walk(mel);
  for( walk.first() ; !walk.done() ; walk.next() ) {
    out << walk.getitem();
  }
  */
  storeIndicies(out);

  delete [] filename;
#ifndef NDEBUG
  cout << " done." << endl;
#endif
}
void
Hmesh::storeIndicies(ostream& out) 
{
  // backup index managers 
  for (int i=0;i<numOfIndexManager2d; ++i) 
  {
    indexmanager[i].backupIndexSet(out);
  }

  // backup vertex indices 
  {
    Listwalk_impl < Vertex > walk(vl) ;
    for( walk.first() ; ! walk.done() ; walk.next() ) 
    {
      int idx=walk.getitem().getIndex();
      out.write( ((const char *) &idx ), sizeof(int) ) ;
    }
  }
  
  // backup element and edge indices 
  {
    Levelwalk<Element> walk(mel, 0);
    for( walk.first() ; !walk.done() ; walk.next() ) 
    {
      SubtreeIterator<Element> hier(&(walk.getitem()));
      for (hier.first() ; !hier.done() ; hier.next() ) 
      {
        // element 
        {
	        int idx=hier.getitem().getIndex();
	        out.write( ((const char *) &idx ), sizeof(int) ) ;
        }

        // edges 
	      for (int e=0;e<3; ++e) 
        {
	        int idx=hier.getitem().edge(e)->getIndex();
	        out.write( ((const char *) &idx ), sizeof(int) ) ;
	      }
      }
    }
  }
}

bool
Hmesh::recoverGrid(const char* recoverFile,
                   double& time, unsigned long int &nbr)
{
  int compwarn = 0;

  cout << "Hmesh::recoverGrid(): trying to read file \""
       << recoverFile << "\" ...";

  char filename[1000];
  sprintf(filename,"%s.refine",recoverFile);

  ifstream in(filename);
  if( !in )
  {
    cout << " FAILED." << endl;
    return false;
  }

  // Zeit und Nummer des Zeitschritts lesen
  // in >> time;
  // in >> nbr;
  // Gitter wiederherstellen
  for( int level = 0 ;; level++ ) {
    {
      Levelwalk<Element> walk(mel, level);
      if( !walk.size() )
        break;
      for( walk.first() ; !walk.done() ; walk.next() ) {
        char flag;
        in.get(flag);
        switch (flag) {
	case Thinelement::unsplit:
	  break;
	case Thinelement::triang_bnd:
	  cerr << "ERROR (Hmesh::recoverGrid()): "
	       << "splitrule \"triang_bnd\" is not allowed for elements!"
	       << endl;
	  abort();
	  break;
	case Thinelement::triang_conf2:
	  walk.getitem().mark(Refco::ref_1);
	  break;
	case Thinelement::triang_quarter:          
	  walk.getitem().mark(Refco::quart);
	  break;
	case Thinelement::compatibility:
	  if (!compwarn)
	    {
              cerr << "WARNING (Hmesh::recoverGrid()): "
		   << "using compatibility mode for obsolete file format!"
		   << endl;
              compwarn = 1;
	    }
	  walk.getitem().mark(Refco::ref_1);
	  break;
	default:
	  cerr << "ERROR (Hmesh::recoverGrid()): "
	       << "unknown splitrule!"
	       << endl;
	  abort();
	}
      }
    }
    refine();
  }


  // read indices 
  recoverIndicies(in);
 
  cout << " done." << endl;

  return true;
}

void
Hmesh::recoverIndicies(istream& in) 
{
  // reads maxIndex of Index Manager 
  for (int i=0;i<numOfIndexManager2d; ++i) 
  {
    indexmanager[i].restoreIndexSet(in);
  }
  
  //////////////////////////////////////////
  //  read vertices 
  //////////////////////////////////////////
  {
    IndexManager2dType& vertexManager = indexmanager[IM_Vertices];
    const int idxSize = vertexManager.getMaxIndex();
    // create vector, all entries are marked true 
    vector<bool> isHole (idxSize, true );

    Listwalk_impl < Vertex > walk(vl) ;
    for( walk.first() ; ! walk.done() ; walk.next() ) 
    {
      Vertex& vx = walk.getitem();
      in.read ( ((char *) &(vx.setIndex())), sizeof(int) );
      assert( vx.getIndex() < idxSize );
      isHole[vx.getIndex()] = false;  
    }

    // all remaining indices are reinserted as holes 
    vertexManager.generateHoles( isHole );
  }
  
  //////////////////////////////////////////
  //  read elements and edges 
  //////////////////////////////////////////
  {
    IndexManager2dType& elementManager = indexmanager[IM_Elements];
    const int elSize = elementManager.getMaxIndex();
    
    IndexManager2dType& edgeManager = indexmanager[IM_Edges];
    const int edgeSize = edgeManager.getMaxIndex();
    // create vector, all entries are marked true
    vector<bool> elementIsHole (elSize, true );
    vector<bool> edgeIsHole  (edgeSize, true );
    
    Levelwalk<Element> walk(mel, 0);
    for( walk.first() ; !walk.done() ; walk.next() ) 
    {
      SubtreeIterator<Element> hier(&(walk.getitem()));
      for (hier.first() ; !hier.done() ; hier.next() ) 
      {
        Element& elem = hier.getitem();

        // read element index 
      	in.read ( ((char *) &(elem.setIndex())), sizeof(int) );
        assert( elem.getIndex() < elSize );
        elementIsHole[elem.getIndex()] = false;
        
        // read edges 
       	for (int e=0; e<3; ++e) 
        {
          int edgeNum = -1; 
	        in.read ( ((char *) &(edgeNum)), sizeof(int) );
          assert( edgeNum < edgeSize );
          edgeIsHole[edgeNum] = false;
          // set edge index 
          elem.edge(e)->setIndex() = edgeNum;
	      }
      }
    } 

    // reinsert remaining indices as holes 
    elementManager.generateHoles( elementIsHole );
    edgeManager.generateHoles( edgeIsHole );
  }
}
