
#ifndef __HEADER__GRID
#define __HEADER__GRID

// ***************************************************
// #begin(header)
// #filename:
//   grid.h
// #description:
//   Basisklassen f"ur Punkte, Gitterelemente und 
//   Hierarchien. 
// #classes:
//   template < class A > class Listagent
//   class Restrict_basic
//   class Basic
//   class Vertex : public Listagent < Vertex >, public Basic
//   class Fullvertex : public Vertex
//   class Thinelement : public Basic
//   class Refco
//   class Refco_el : private Refco
//   class Element : public Thinelement, public Refco_el
//   template < class A > class Hier : public A
//   class Bndel : public Thinelement, protected Refco
// #copyright:
//   IAM Freiburg
// #end(header)
// ***************************************************
#define _DUNE_USES_ALU3DGRID_ 

#include <stdio.h>

#include <stdlib.h>
#include <assert.h>
#include <cmath>

#ifdef _ANSI_HEADER
  using namespace std;
  #include <vector>
  #include <fstream>
  #include <iostream>
#else
  #include <vector.h>
  #include <fstream.h>
  #include <iostream.h>
#endif

#include "xdisplay.h"
#include "vtx_btree.h"

class Hmesh_basic;
class Basic {

  private :

    int refcount ;

    Basic(const Basic & ) ;

    Basic & operator=(const Basic &) ;

  protected :

    Basic() : refcount(0), hdl(0), _idx(-1) {}

    virtual ~Basic() { assert(!refcount) ; }

    int _idx;
    Hmesh_basic *hdl;
    int &setIndex () { return _idx ; }
 public :
    virtual void sethdl(Hmesh_basic *phdl) = 0 ; // {hdl=phdl;}
    inline int getIndex () const { return _idx; }

    enum {  ncoord = 2,

            nparts = 4, max_points = 4 } ;

    void attach() { refcount ++ ; }

    void detach() { refcount -- ; }

    int isfree() const { return refcount == 0 ; }


    virtual void write(ofstream &) const = 0 ;

    virtual void read(ifstream &) = 0 ;

} ;

class Vertex;
class Multivertexadapter ;
template < class A > class Listagency ;
template <class A> class Hier;

// ***************************************************
// #begin(class)
// #description:
// #definition: 
template < class A > class Listagent {

  A * next_list ;
  
  A * prev_list ;

  int number_list ;
    
  Listagency < A > * agnc ;
  
  Listagent(const Listagent &) { }

  Listagent & operator=(const Listagent &) { }

  protected :
  
    Listagent() : next_list(0), prev_list(0), agnc(0) { }
     
  public :
  
    virtual ~Listagent() { }
  
    A * next() const { return next_list ; }
    
    A * prev() const { return prev_list ; }

    int number() const { return number_list ; }
          
  friend class Listagency < A > ;
  
} ;
// #end(class)
// ***************************************************


class Element;

class Restrict_basic;
class Prolong_basic;

typedef struct nconf_vtx
{
  Vertex *vtx;
  Hier<Element> *el[2];
  nconf_vtx(Vertex *v,Hier<Element> *e1,Hier<Element> *e2) 
   {vtx=v;el[0]=e1;el[1]=e2;} 
} nconf_vtx_t;

// ***************************************************


// ***************************************************
// #begin(class)
// #description:
//   Interfaceklasse f"ur Punkte
// #definition:
class Vertex : public Listagent < Vertex >, public Basic {

  private :

    Vertex & operator=(const Vertex &) ;
    
    Vertex(const Vertex &) ;

 protected :

    int nr_of_pernbs;

    Vertex *pernb[3];

    int _level;

  protected :
    Vertex() {
      nr_of_pernbs = 0;
      pernb[0] = (Vertex *)0;
      pernb[1] = (Vertex *)0;
      pernb[2] = (Vertex *)0;
      _level = -1;
    }

  public :

    virtual ~Vertex(); 

   
    virtual const double (& coord() const )[ncoord] = 0 ;


    virtual void write(ofstream &) const = 0 ;


    virtual void read(ifstream &) = 0 ;
 

    int get_nr_of_per_nbs()
    {
      return nr_of_pernbs;
    }

    void set_pernb(Vertex *pv) 
    {
      int li, lnew = 1;

      for (li=0;li<nr_of_pernbs;li++)
      {
        if (pv == pernb[li]) lnew = 0;
      }

      if (pv == this) lnew = 0;

      if (lnew)
      {
        assert(nr_of_pernbs < 3);

        pernb[nr_of_pernbs] = pv;
        nr_of_pernbs++;
      }
    }

    int no_pernbs() {
      return nr_of_pernbs;
    }

    Vertex *get_pernb(int pidx)
    {
      Vertex *lret;

      if (pidx < nr_of_pernbs) lret = pernb[pidx];
      else lret = (Vertex *)0;

      return lret; 
    } 
    
    int level() {
      return _level;
    }
} ;
// #end(class)
// ***************************************************


// ***************************************************
// #begin(class)
// #description:
//  Punkte
// #definition:
class Fullvertex : public Vertex {

  private:

    double vcoord[ncoord] ;
  
    Fullvertex(const Fullvertex &) ;

    Fullvertex & operator = (const Fullvertex &) ;

  public :

    virtual void sethdl(Hmesh_basic *phdl);

    Fullvertex() {}

    Fullvertex(double (&)[ncoord],int ) ;

    Fullvertex(double , double, int ) ;

   ~Fullvertex() { }
  
    const double (& coord() const )[ncoord] { return vcoord ; }

    void write(ofstream &) const ;

    void read(ifstream &) ;

    friend ostream& operator<<(ostream& os, const Fullvertex& fv) {
      return os << "(" << fv.vcoord[0] << "," << fv.vcoord[1] << ")";
    }
} ;
// #end(class)
// ***************************************************
class Edge : public Basic {
  virtual void sethdl(Hmesh_basic *phdl);
 public:
    Edge(Hmesh_basic *phdl) {
      sethdl(phdl);
    }
    ~Edge();
    void write(ofstream &) const ;
    void read(ifstream &) ;
};
class Triang;
class Bndel_triang;
// ***************************************************
// #begin(class)
// #description:
//   Basisklasse f"ur Verfeinerungsinformation
// #definition:
class Refco {

  public :
    
  typedef enum { crs = -1 , ref = 1 , 
                 none = 0, crs_1 = -2, crs_2 = -3, crs_3 = -4, crs_4 = -5,

                 ref_1 = 2, ref_2 = 3, ref_3 = 4, ref_4 = 5, notref = 10,

		 quart=9, notcrs = -10  } tag_t ;
      
    void clear(tag_t ) { }
      
    int is(tag_t ) const { return 0 ; }  // ist nie irgendein tag

    void mark(Refco::tag_t t) {  }
      
} ;
// #end(class)
// ***************************************************


// ***************************************************
// #begin(class)
// #description:
//   Interfaceklasse f"ur Elemente
// #definition:
class Thinelement : public Basic {

  public :

    typedef enum { bndel_like, element_like } thintype_t ;

    typedef enum { unsplit = 0, triang_bnd = -2, compatibility = 1,
                   triang_conf2 = 2, triang_quarter = 4 } splitrule_t ;

  private :
  
    Thinelement(const Thinelement & ) ;

    Thinelement & operator=(const Thinelement &) ;

  protected :
  
    splitrule_t mysplit ;

    int nfaces, nedges, nvertices ;
  
    Thinelement() : mysplit(unsplit), nfaces(0), nedges(0), nvertices(0) { }
    
    virtual ~Thinelement() { }
   
  public :
    
    virtual int thinis(thintype_t ) const = 0 ;

    splitrule_t splitrule() const { return mysplit ; }
     
    int numfaces() const { return nfaces ; }

    int numedges() const { return nedges ; }

    int numvertices() const { return nvertices ; }

    virtual double area() const = 0;

    virtual void edge_vtx(int , Vertex * (&) [2]) const = 0 ;


    virtual int numfacevertices(int ) const = 0 ;
    
    virtual int facevertex(int , int ) const = 0 ;


    virtual Vertex * vertex(int ) const = 0 ;
    
    inline Vertex * vertex(int fce, int j) const { return vertex(facevertex(fce, j)) ; }


    virtual Thinelement * neighbour(int ) const = 0 ;

    virtual int opposite(int fce) const = 0 ;

    virtual int edge_idx(int fce) const = 0 ;

    virtual Edge *edge(int fce) const = 0 ;

    virtual void nbconnect(int , Thinelement *, int ) = 0 ;


    virtual void write(ofstream &) const = 0 ;

    virtual void read(ifstream &, Vertex **, const int ) = 0 ;

    virtual int split(void * (&) [nparts], Listagency < Vertex > *,
                      Multivertexadapter &,nconf_vtx_t *,splitrule_t,
		      int,Refco::tag_t,Prolong_basic *pro_el) = 0 ;

    virtual int docoarsen(nconf_vtx_t*,int,Restrict_basic *rest_el) { return 1; };

    virtual void draw(Xdisplay & ) {};

    Triang * nbel(int l) const
    {
      Thinelement *el=neighbour(l);
      if (el->thinis(element_like))
        return (Triang*)el;
      else
        return NULL;
    }   
    Bndel_triang * nbbnd(int l) const
    {
      Thinelement *el=neighbour(l);
      //if (el==NULL)
      //return NULL;
      assert(el!=NULL);
      if (el->thinis(bndel_like))
        return (Bndel_triang*)el;
      else
        return NULL;
    }   

} ;
// ***************************************************

// ***************************************************
// #begin(class)
// #description:
//   Verfeinerungsinformation
// #definition:
class Refco_el : protected Refco {

//public :

//  Refco::tag_t ;

  private :

    Refco::tag_t tag, tag_last ;   
    
  public :

    Refco_el() : tag(none), tag_last(none) { }

    void clear(Refco::tag_t t = none) { tag = (t == tag) ? none : tag ; }
      
    void mark(Refco::tag_t t) { tag_last = tag; tag = t ; }
      
    int is(Refco::tag_t t) const { return tag == t ? 1 : 0 ; }

    // so kann diese Funktion nicht bleiben...
    int was(Refco::tag_t t) const { return tag_last == t ? 1 : 0 ; }    

} ;
// #end(class)
// ***************************************************


// ***************************************************
// #begin(class)
// #description:
//   Elementinterface mit verfeinerungsinfo.
// #definition:

template<class T>
class Hier;

class Element : public Thinelement, public Refco_el {

  private :
             
    Element(const Element &) ;

    Element & operator = (const Element &) ;
       
  public :

    virtual void sethdl(Hmesh_basic *phdl);
    int thinis(thintype_t t) const { return t == element_like ; }

  protected:

    struct c {	//  die vertex und (thin-)element verbindungen

      enum {pv=2,nv=3,nf=3};

      Vertex  * vtx [3] ;

      Thinelement * nb [3] ;

      Edge * edge [3];

      mutable Vtx_btree* hvtx[3];

      short int bck [3] ;
      
      short int normdir [3] ;

      c() ;

     ~c() ;

      void set(Vertex * v, int i) { (vtx[i] = v)->attach() ; }

      void unset(int i) { vtx[i]->detach() ; vtx[i] = 0 ; }

      void write(ofstream &) const ;

      void read(ifstream &, Vertex ** , const int ) ;

      int check();

      friend class Hmesh_basic ;

    } connect ;

    double _area;
    double _minheight;
    double _outernormal[3][2];
    double _sidelength[3];

  public :
   
    Element() : _area(-1.0), _minheight(-1.0)
    {
      int i;

      for (i=0;i<3;i++)
      {
        _outernormal[i][0] =  0.0;
        _outernormal[i][1] =  0.0;
        _sidelength[i]     = -1.0;
      }
    }

    virtual ~Element();

    int numfacevertices(int ) const { return connect.pv ; }
    
    int facevertex(int , int ) const ;
      
    void edge_vtx(int e, Vertex * (& ) [2] ) const ;


    Vertex * vertex(int ) const ;

    Vertex * vertex(int fce, int j) const { return vertex(facevertex(fce, j)) ; }

    
    Thinelement * neighbour(int ) const ;

    int opposite(int ) const ;

    int edge_idx(int ) const ;

    Edge *edge(int ) const ;

    int normaldir(int ) const ;

    void nbconnect(int , Thinelement * , int ) ;

    void edgeconnect(int, Edge *) ;

    void setnormdir(int , int ) ;

    void setrefine(int ) ;
  
    void init() ;

    int setrefine() ;

    int setorientation();


    double sidelength(int pfce) const { return _sidelength[pfce%connect.nf]; }

    double minheight() const { return _minheight; }

    virtual double area() const { return _area; }

    void outernormal(int ,double (& )[ncoord]) const;

    void dirnormal(int ,double (& )[ncoord]) const;

    void tolocal(const double (& )[ncoord],double (& )[3]) const;

    void fromlocal(const double (& )[3],double (& )[ncoord]) const; // Fullvertex nehmen

    void midpoint(int ,double (& )[3]) const;

    void facepoint(int ,double ,double (& )[3]) const;

    int isecpoint(int , const double (& )[ncoord], double (& )[ncoord]) const; 

    int isecpoint(int , const double (& )[ncoord], double (& )[3]) const; 

    int inside(Vertex &) const ;

    int inside(const double (& )[2]) const;

    void addhvtx(Vertex* inv, Thinelement *lnb,Thinelement *rnb,int fce);

    int hashvtx(int fce) const {
      return (connect.hvtx[fce] != 0);
    }

    void removehvtx(int fce,Vertex *vtx) {
      if (connect.hvtx[fce]->count()==1) {
	assert(connect.hvtx[fce]->getHead()==vtx);
	delete connect.hvtx[fce];
	connect.hvtx[fce] = 0;
      } else {
	bool found=connect.hvtx[fce]->remove(vtx);
	assert(found);
      }
    }

    int check() ;

    void draw(Xdisplay & ) ;

    friend ostream&
      operator<< (ostream& out, const Element& elem) {
        return out;
      }

    friend istream&
      operator>> (istream& in, Element& elem) {
	return in;
      }
} ;
// #author:
// #end(class)
// ***************************************************

template <class T>
class SubtreeIterator;

// ***************************************************
// #begin(class)
// #description:
//   Hierachieklasse
// #definition:
template < class A > class Hier : public A {
 
  Hier * dwn ;

  Hier * nxt ;
  
  Hier * up;

  int lvl ;

  int childNr_;

  protected :

    int numchild ;

    Hier() : dwn(0), nxt(0), up(0), lvl(0) { }

    void deletesubtree() { delete dwn ; dwn = 0; 
      //  this->check();
    };  // Wird f"ur konf. Dreiecke gebraucht

  public :

    virtual ~Hier() {

      if(dwn) delete dwn ;

      if(nxt) delete nxt ;

    }

    Hier * down() const { return dwn ; }

    Hier * father() const { return up ; }

    Hier * next() const { return nxt ; }

    int leaf() const { return ! dwn ; }

    int level() const { return lvl ; }

    int childNr() const { return childNr_; }

    int nchild() const { return numchild ; }

    int count() const { return (nxt ? nxt->count() : 0) + (dwn ? dwn->count() : 1) ; }

    int count(int i) const { 

      return (nxt ? nxt->count(i) : 0) + 

             (lvl == i ? 1 : (dwn ? dwn->count(i) : 0) );

    }

    SubtreeIterator<A> stIterator() {
      return SubtreeIterator<A>(this);
    }

    int deepestLevel() {
      SubtreeIterator<A> iter = stIterator();
      int currLevel = iter->level();
      int result = currLevel;
      while( ++iter ) {
        if( iter->level() > result )
        result = iter->level();
      }
      return result - currLevel;
    }

    int coarse(nconf_vtx_t *ncv,int nconfDeg,Restrict_basic *rest_el) {

      if(dwn ? dwn->coarse(ncv,nconfDeg,rest_el) == numchild : 0 )

        if ( this->docoarsen(ncv,nconfDeg,rest_el) )
	{
          this->deletesubtree();
          this->mysplit = this->unsplit;
	}

      int i = (nxt ? nxt->coarse(ncv,nconfDeg,rest_el) : 0 ) + (this->is(Refco::crs) ? 1 : 0 ) ;

      this->clear(Refco::crs) ;

      return i ;

    }

    int refine_leaf(Listagency < Vertex > * a, 
		    Multivertexadapter * b ,nconf_vtx_t *ncv,
		    int nconfDeg,Refco::tag_t default_ref,
		    Prolong_basic *pro_el)  {

      int count = 0;

      assert( leaf() );

      if (this->is(Refco::ref))
	this->mark(default_ref);

      if(this->is(Refco::quart) || 
	 this->is(Refco::ref_1) || this->is(Refco::ref_2) || 
	 this->thinis(this->bndel_like)) {
          
        void * els [this->nparts] ;

        if (this->is(Refco::quart))
        {
          numchild = this->split(els, a, *b, ncv, this->triang_quarter,nconfDeg,default_ref,pro_el);
          this->clear(Refco::quart);
        }
        else if (this->is(Refco::ref_1))
        {
          numchild = split(els, a, *b, ncv, this->triang_conf2,nconfDeg,default_ref,pro_el);
          this->clear(Refco::ref_1);
        }
        else if (this->is(Refco::ref_2))
        {
          numchild = split(els, a, *b, ncv, this->triang_conf2,nconfDeg,default_ref,pro_el);
          this->clear(Refco::ref_2);
        }
        else
	{
          assert(this->thinis(this->bndel_like));
          this->numchild = split(els, a, *b, ncv, this->triang_bnd,nconfDeg,default_ref,pro_el);
	}

        dwn = (Hier *)els[0] ;

        dwn->lvl = lvl + 1 ;

	dwn->up = this;

	dwn->childNr_ = 0;

        for(int i = 1 ; i < numchild ; i ++ ) {

          ((Hier *)els[i])->lvl = lvl + 1 ;

          ((Hier *)els[i])->up = this ;

	  ((Hier *)els[i])->childNr_ = i;

          ((Hier *)els[i-1])->nxt = (Hier *)els[i] ;

        }
	
	//this->check();

        count = numchild;

      }

      return count;

    }

    int refine(Listagency < Vertex > * a, Multivertexadapter * b,
	       nconf_vtx_t *ncv,
	       int nconfDeg,Refco::tag_t default_ref,Prolong_basic *pro_el) {

      int count =  nxt ? nxt->refine(a, b,ncv, nconfDeg,default_ref,pro_el) : 0 ;
      
      if(dwn) count += dwn->refine(a, b,ncv,nconfDeg,default_ref,pro_el) ;

      else {

 // Neue Behandlung der Bl"atter:
 // Wegen rek. Aufbau der Dreiecksverf. ist eine Funktion n"otig, die Verf.
 // aber nicht u"ber den Baum l"auft. 
 // Weitere "Anderung: count+= statt count=, falls n"amlich von der nxt-Rek.
 // etwas in count steht.
 // Bei Rekursivem Verf. stimmt R"uckgabe sowieso nicht

        count += refine_leaf(a,b,ncv,nconfDeg,default_ref,pro_el) ; 

      }

      return count ;

    }


    void write(ofstream & ) const { }

    void read(ifstream & ) { }

} ;
// #end(class)
// ***************************************************

// ***************************************************
// #begin(class)
// #description:
//   Interfaceklasse f"ur Randelemente
// #definition:
class Bndel : public Thinelement, public Refco {

  protected :

  struct c {

    enum {nf=1,nv=2};

    Vertex * vtx[max_points] ;

    Thinelement * nb ;

    short int bck ;

    Edge *edge;

    c() ;

   ~c() ;

    void set(Vertex * a, int i) { (vtx[i] = a)->attach() ; }

    void unset(int i) { vtx[i]->detach() ; vtx[i] = 0 ; }

    void write(ofstream &) const ;

    void read(ifstream &, Vertex ** , const int ) ;

  } connect ;

  public :
  
  enum {periodic = 1111,general_periodic = 2222};
    typedef int bnd_t ;

    // typedef int bnd_part_t[max_bndnr+offset_bndnr+1];

    virtual void sethdl(Hmesh_basic *phdl);
  private :
      
    Bndel(const Bndel &) ;
    
    Bndel & operator = (const Bndel &) ;
 
  protected :

    Bndel(bnd_t t = none) : typ(t) { }

    bnd_t typ ;

    double (*lf)(double);

    double (*lDf)(double);
 
  public :

    void check() {}
    
    virtual ~Bndel();

    bnd_t type() const { return typ ; }

    int thinis(thintype_t t) const { return t == bndel_like ; }

    int facevertex(int , int ) const ;

    int numfacevertices(int ) const { return connect.nv ; }

 
    void edge_vtx(int e, Vertex * (& ) [2] ) const ;


    Vertex * vertex(int ) const ;
 
    Vertex * vertex(int , int j) const { return vertex(j) ; }


    Thinelement * neighbour(int ) const { return connect.nb ; }

    int neighbours(int ) const { return 1; }

    int opposite(int ) const { return connect.bck ; }

    int edge_idx(int ) const { return connect.edge->getIndex(); }

    Edge *edge(int ) const { return connect.edge; }

    void nbconnect(int , Thinelement * , int ) ;
    void edgeconnect(int , Edge *) ;

    void set_bndfunctions(double (*pf)(double), double (*pDf)(double))
    {
      lf  = pf;
      lDf = pDf;
    }

    int get_splitpoint(double (& ) [2]) ;

    void mirror_opposite(double (& )[2]) const;

    void mirror(const double (& )[2], double (& )[2]) const ;

    virtual double area() const ;

    virtual Bndel *create(Vertex * , Vertex *,bnd_t) const = 0;

    int setorientation() 
    {
      int ret = (connect.vtx[0] != connect.nb->vertex(connect.bck+1));
      if(ret)
      {
        Vertex *tmpv=connect.vtx[0];
        connect.vtx[0]=connect.vtx[1];
        connect.vtx[1]=tmpv;
      }
      return ret;
    }
 
    void draw(Xdisplay & ) ; 
} ;
// #end(class)
// ***************************************************



class Restrict_basic
{
  public: 
  virtual void operator ()(Element *parent,Element **child,
                           int anzahl)
  {
  }
};

class Prolong_basic
{
  public:
  virtual void operator ()(const Element *parent,Element **child,
                           int anzahl)
  {
  }
};

#endif


