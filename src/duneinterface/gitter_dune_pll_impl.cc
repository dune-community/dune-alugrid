#ifndef GITTER_DUNE_PLL_IMPL_CC_INCLUDED
#define GITTER_DUNE_PLL_IMPL_CC_INCLUDED

#include "gitter_dune_pll_impl.h"
#include "gitter_dune_pll_mgb.cc"

IteratorSTI < Gitter :: helement_STI > * GitterDunePll :: 
leafIterator (const helement_STI *)
{
  return new Insert < PureElementAccessIterator < Gitter :: helement_STI > :: Handle,
    TreeIterator < Gitter :: helement_STI, is_leaf < Gitter :: helement_STI> > > (container ()) ;
}

IteratorSTI < Gitter :: helement_STI > * GitterDunePll ::
leafIterator (const IteratorSTI < helement_STI > * p)
{
  return new Insert < PureElementAccessIterator < Gitter :: helement_STI > :: Handle,
    TreeIterator < Gitter :: helement_STI, is_leaf < Gitter :: helement_STI> > >
    (*(const Insert < PureElementAccessIterator < Gitter :: helement_STI > :: Handle,
       TreeIterator < Gitter :: helement_STI, is_leaf < Gitter :: helement_STI> > > *) p) ;
}


bool GitterDunePll :: duneNotifyNewGrid ()
{
  assert (debugOption (20) ? (cout << "**GitterDunePll :: duneNotifyNewGrid () " << endl, 1) : 1) ;
  const int np = mpAccess ().psize () ;
  LoadBalancer :: DataBase db ;
  {
    AccessIterator < hface_STI > :: Handle w (containerPll ()) ;
    for (w.first () ; ! w.done () ; w.next ()) w.item ().accessPllX ().ldbUpdateGraphEdge (db) ;
  }
  {
    AccessIterator < helement_STI > :: Handle w (containerPll ()) ;
    for (w.first () ; ! w.done () ; w.next ()) w.item ().accessPllX ().ldbUpdateGraphVertex (db) ;
  }
  bool neu = false ;
  {
    // Kriterium, wann eine Lastneuverteilung vorzunehmen ist:
    // 
    // load  - eigene ElementLast
    // mean  - mittlere ElementLast
    // nload - Lastverh"altnis

    double load = db.accVertexLoad () ;
    vector < double > v (mpAccess ().gcollect (load)) ;
    double mean = accumulate (v.begin (), v.end (), 0.0) / double (np) ;

    for (vector < double > :: iterator i = v.begin () ; i != v.end () ; i ++)
      neu |= (*i > mean ? (*i > (_ldbOver * mean) ? true : false) : (*i < (_ldbUnder * mean) ? true : false)) ;
  }
  return neu;
}



void GitterDunePll :: duneNotifyGridChanges ()
{
  Gitter :: notifyGridChanges () ;
  duneExchangeDynamicState () ;
  return ;
}


// done call notify and loadBalancer  
bool GitterDunePll :: dAdapt ()   
{
  __STATIC_myrank = mpAccess ().myrank () ;
  __STATIC_turn ++ ;
  assert (debugOption (20) ? (cout << "**INFO GitterDunePll :: dAdapt ()" << endl, 1) : 1) ;
  assert (! iterators_attached ()) ;
  int start = clock () ;
  //bool refined = 
  this->refine() ;
  int lap = clock () ;
  this->coarse ();
  int end = clock () ;
  if (debugOption (1))
    {
      float u1 = (float)(lap - start)/(float)(CLOCKS_PER_SEC) ;
      float u2 = (float)(end - lap)/(float)(CLOCKS_PER_SEC) ;
      float u3 = (float)(end - start)/(float)(CLOCKS_PER_SEC) ;
      cout << "**INFO GitterDunePll :: adapt () [ref (loops)|cse|all] " << u1 << " ("
     << _refineLoops << ") " << u2 << " " << u3 << endl ;
    }
  duneNotifyGridChanges () ;
  //balanceGrid_ = duneNotifyNewGrid();
  
  return true;
  //return refined;
}

// done call notify and loadBalancer  
bool GitterDunePll :: duneAdapt (AdaptRestrictProlongType & arp)   
{
  this->setAdaptRestrictProlongOp(arp);
  bool refined = this->dAdapt();
  this->removeAdaptRestrictProlongOp ();
  return refined;
}

bool GitterDunePll :: duneLoadBalance () 
{
  loadBalancerGridChangesNotify () ;
  return true;
}

// returns true if grid was repartitioned 
bool GitterDunePll :: duneLoadBalance (GatherScatterType & gs, AdaptRestrictProlongType & arp) {
  
  this->setAdaptRestrictProlongOp(arp);
  assert (debugOption (20) ? (cout << "**GitterDunePll :: duneLoadBalance () " << endl, 1) : 1) ;
  const int np = mpAccess ().psize () ;
  LoadBalancer :: DataBase db ;
  {
    AccessIterator < hface_STI > :: Handle w (containerPll ()) ;
    for (w.first () ; ! w.done () ; w.next ()) w.item ().accessPllX ().ldbUpdateGraphEdge (db) ;
  }
  {
    AccessIterator < helement_STI > :: Handle w (containerPll ()) ;
    for (w.first () ; ! w.done () ; w.next ()) w.item ().accessPllX ().ldbUpdateGraphVertex (db) ;
  }
  bool neu = false ;
  {
    // Kriterium, wann eine Lastneuverteilung vorzunehmen ist:
    // 
    // load  - eigene ElementLast
    // mean  - mittlere ElementLast
    // nload - Lastverh"altnis
  
    double load = db.accVertexLoad () ;
    vector < double > v (mpAccess ().gcollect (load)) ;
    double mean = accumulate (v.begin (), v.end (), 0.0) / double (np) ;

    for (vector < double > :: iterator i = v.begin () ; i != v.end () ; i ++)
      neu |= (*i > mean ? (*i > (_ldbOver * mean) ? true : false) : (*i < (_ldbUnder * mean) ? true : false));
  }
  if (neu) 
    {
      if (mpAccess ().gmax (_ldbMethod)) 
  {
    assert (debugOption (5) ? (cout << "**GitterDunePll :: repartitioning macro grid! " << endl, 1) : 1) ;
    duneRepartitionMacroGrid (db, gs) ;
    notifyMacroGridChanges () ;
  }
    }
  this->removeAdaptRestrictProlongOp ();
  return neu;
}

void GitterDunePll :: duneExchangeDynamicState () 
{
  // Die Methode wird jedesmal aufgerufen, wenn sich der dynamische
  // Zustand des Gitters ge"andert hat: Verfeinerung und alle Situationen
  // die einer "Anderung des statischen Zustands entsprechen. Sie wird in
  // diesem Fall NACH dem Update des statischen Zustands aufgerufen, und
  // kann demnach von einem korrekten statischen Zustand ausgehen. F"ur
  // Methoden die noch h"aufigere Updates erfordern m"ussen diese in der
  // Regel hier eingeschleift werden.
  {

    const int nl = mpAccess ().nlinks () ;
  
#ifndef NDEBUG 
    // if debug mode, then count time 
    const int start = clock () ;
#endif
  
    try 
      {
  typedef Insert < AccessIteratorTT < hface_STI > :: InnerHandle,
    TreeIterator < hface_STI, is_def_true < hface_STI > > > InnerIteratorType;
  typedef Insert < AccessIteratorTT < hface_STI > :: OuterHandle, 
    TreeIterator < hface_STI, is_def_true < hface_STI > > > OuterIteratorType;
                
  vector < ObjectStream > osv (nl) ;
  {
    for (int l = 0 ; l < nl ; l ++) 
      {
        {
    AccessIteratorTT < hface_STI > :: InnerHandle mif (this->containerPll (),l) ;
    AccessIteratorTT < hface_STI > :: OuterHandle mof (this->containerPll (),l) ;

    InnerIteratorType wi (mif);
    for (wi.first () ; ! wi.done () ; wi.next ()) 
      {
        pair < ElementPllXIF_t *, int > p = wi.item ().accessPllX ().accessInnerPllX () ;
        p.first->writeDynamicState (osv [l], p.second) ;
      }
        
    OuterIteratorType wo (mof);
    for (wo.first () ; ! wo.done () ; wo.next ()) 
      {
        pair < ElementPllXIF_t *, int > p = wo.item ().accessPllX ().accessInnerPllX () ;
        p.first->writeDynamicState (osv [l], p.second) ;
      }
        }
      } 
  }
    
  osv = mpAccess ().exchange (osv) ;
    
  { 
    for (int l = 0 ; l < nl ; l ++ ) 
      {
        {
    AccessIteratorTT < hface_STI > :: OuterHandle mof (this->containerPll (),l) ;
    AccessIteratorTT < hface_STI > :: InnerHandle mif (this->containerPll (),l) ;
        
    OuterIteratorType wo (mof) ;
    for (wo.first () ; ! wo.done () ; wo.next ()) 
      {
        pair < ElementPllXIF_t *, int > p = wo.item ().accessPllX ().accessOuterPllX () ;
        p.first->readDynamicState (osv [l], p.second) ;
      }
        
    InnerIteratorType wi (mif);
    for (wi.first () ; ! wi.done () ; wi.next ()) 
      {
        pair < ElementPllXIF_t *, int > p = wi.item ().accessPllX ().accessOuterPllX () ;
        p.first->readDynamicState (osv [l], p.second) ;
      }
        }
      } 
  }
      } 
    catch (Parallel ::  AccessPllException) 
      {
  cerr << "  FEHLER Parallel :: AccessPllException entstanden in: " << __FILE__ << " " << __LINE__ << endl ;
      }
    assert (debugOption (20) ? (cout << "**INFO GitterDunePll :: exchangeDynamicState () used " << (float)(clock () - start)/(float)(CLOCKS_PER_SEC) << " sec. " << endl, 1) : 1 ) ;
  }
}

#if 0
// reun only over leaf Level 
void GitterDunePll :: duneExchangeDataLeaf (GatherScatterType & gs) 
{
  // Die Methode wird jedesmal aufgerufen, wenn sich der dynamische
  // Zustand des Gitters ge"andert hat: Verfeinerung und alle Situationen
  // die einer "Anderung des statischen Zustands entsprechen. Sie wird in
  // diesem Fall NACH dem Update des statischen Zustands aufgerufen, und
  // kann demnach von einem korrekten statischen Zustand ausgehen. F"ur
  // Methoden die noch h"aufigere Updates erfordern m"ussen diese in der
  // Regel hier eingeschleift werden.
  {
    const int nl = mpAccess ().nlinks () ;
#ifndef NDEBUG
    const int start = clock () ;
#endif
    try 
      {
  typedef LeafIteratorTT< hface_STI > InnerIteratorType;
  typedef LeafIteratorTT< hface_STI > OuterIteratorType;
                
  vector < ObjectStream > osv (nl) ;
  {
    for (int l = 0 ; l < nl ; l ++) 
      {
        {
    LeafIteratorTT < hface_STI > w (*this,l);
    for (w.inner() .first () ; ! w.inner().done () ; w.inner().next ()) 
      {
        pair < ElementPllXIF_t *, int > p = w.inner().item ().accessPllX ().accessInnerPllX () ;
        p.first->writeDynamicState (osv [l], p.second) ;
        p.first->writeDynamicState (osv [l], gs) ;
      }
        
    for (w.outer().first () ; ! w.outer().done () ; w.outer().next ()) 
      {
        pair < ElementPllXIF_t *, int > p = w.outer().item ().accessPllX ().accessInnerPllX () ;
        p.first->writeDynamicState (osv [l], p.second) ;
        p.first->writeDynamicState (osv [l], gs) ;
      }
        }
      } 
  }
    
  osv = mpAccess ().exchange (osv) ;
    
  { 
    for (int l = 0 ; l < nl ; l ++ ) 
      {
        {
    LeafIteratorTT< hface_STI > w (*this,l) ;
    for (w.outer().first () ; ! w.outer().done () ; w.outer().next ()) 
      {
        pair < ElementPllXIF_t *, int > p = w.outer().item ().accessPllX ().accessOuterPllX () ;
        p.first->readDynamicState (osv [l], p.second) ;
        p.first->readDynamicState (osv [l], gs) ;
      }
        
    for (w.inner().first () ; ! w.inner().done () ; w.inner().next ()) 
      {
        pair < ElementPllXIF_t *, int > p = w.inner().item ().accessPllX ().accessOuterPllX () ;
        p.first->readDynamicState (osv [l], p.second) ;
        p.first->readDynamicState (osv [l], gs ) ;
      }
        }
      } 
  }
      } 
    catch (Parallel ::  AccessPllException) 
      {
  cerr << "  FEHLER Parallel :: AccessPllException entstanden in: " << __FILE__ << " " << __LINE__ << endl ;
      }
    assert (debugOption (20) ? (cout << "**INFO GitterDunePll :: duneExchangeData () used " << (float)(clock () - start)/(float)(CLOCKS_PER_SEC) << " sec. " << endl, 1) : 1 ) ;
  }
}
#endif

// go over all levels 
void GitterDunePll :: duneExchangeDataAll (GatherScatterType & gs) 
{
  // Die Methode wird jedesmal aufgerufen, wenn sich der dynamische
  // Zustand des Gitters ge"andert hat: Verfeinerung und alle Situationen
  // die einer "Anderung des statischen Zustands entsprechen. Sie wird in
  // diesem Fall NACH dem Update des statischen Zustands aufgerufen, und
  // kann demnach von einem korrekten statischen Zustand ausgehen. F"ur
  // Methoden die noch h"aufigere Updates erfordern m"ussen diese in der
  // Regel hier eingeschleift werden.
  {
    const int nl = mpAccess ().nlinks () ;
#ifndef NDEBUG
    const int start = clock () ;
#endif
    try 
      {
  typedef Insert < AccessIteratorTT < hface_STI > :: InnerHandle,
    TreeIterator < hface_STI, is_def_true < hface_STI > > > InnerIteratorType;
  typedef Insert < AccessIteratorTT < hface_STI > :: OuterHandle, 
    TreeIterator < hface_STI, is_def_true < hface_STI > > > OuterIteratorType;
                
  vector < ObjectStream > osv (nl) ;
  {
    for (int l = 0 ; l < nl ; l ++) 
      {
        {
    AccessIteratorTT < hface_STI > :: InnerHandle mif (this->containerPll (),l) ;
    AccessIteratorTT < hface_STI > :: OuterHandle mof (this->containerPll (),l) ;

    InnerIteratorType wi (mif);
    for (wi.first () ; ! wi.done () ; wi.next ()) 
      {
        pair < ElementPllXIF_t *, int > p = wi.item ().accessPllX ().accessInnerPllX () ;
        p.first->writeDynamicState (osv [l], p.second) ;
        p.first->writeDynamicState (osv [l], gs) ;
      }
        
    OuterIteratorType wo (mof);
    for (wo.first () ; ! wo.done () ; wo.next ()) 
      {
        pair < ElementPllXIF_t *, int > p = wo.item ().accessPllX ().accessInnerPllX () ;
        p.first->writeDynamicState (osv [l], p.second) ;
        p.first->writeDynamicState (osv [l], gs) ;
      }
        }
      } 
  }
    
  osv = mpAccess ().exchange (osv) ;
    
  { 
    for (int l = 0 ; l < nl ; l ++ ) 
      {
        {
    AccessIteratorTT < hface_STI > :: OuterHandle mof (this->containerPll (),l) ;
    AccessIteratorTT < hface_STI > :: InnerHandle mif (this->containerPll (),l) ;
        
    OuterIteratorType wo (mof) ;
    for (wo.first () ; ! wo.done () ; wo.next ()) 
      {
        pair < ElementPllXIF_t *, int > p = wo.item ().accessPllX ().accessOuterPllX () ;
        p.first->readDynamicState (osv [l], p.second) ;
        p.first->readDynamicState (osv [l], gs) ;
      }
        
    InnerIteratorType wi (mif);
    for (wi.first () ; ! wi.done () ; wi.next ()) 
      {
        pair < ElementPllXIF_t *, int > p = wi.item ().accessPllX ().accessOuterPllX () ;
        p.first->readDynamicState (osv [l], p.second) ;
        p.first->readDynamicState (osv [l], gs ) ;
      }
        }
      } 
  }
      } 
    catch (Parallel ::  AccessPllException) 
      {
  cerr << "  FEHLER Parallel :: AccessPllException entstanden in: " << __FILE__ << " " << __LINE__ << endl ;
      }
    assert (debugOption (20) ? (cout << "**INFO GitterDunePll :: duneExchangeData () used " << (float)(clock () - start)/(float)(CLOCKS_PER_SEC) << " sec. " << endl, 1) : 1 ) ;
  }
}

void GitterDunePll :: duneExchangeData (GatherScatterType & gs, bool leaf) 
{
  //if(leaf) 
  //  this->duneExchangeDataLeaf(gs);
  //else 
  this->duneExchangeDataAll(gs);
  return; 
}

pair < IteratorSTI < GitterPll :: vertex_STI > *, IteratorSTI < GitterPll :: vertex_STI > *> 
GitterDunePll :: borderIteratorTT (const vertex_STI * v, int link )
{
  // return default vertex iterator 
  return this->iteratorTT(v, link);
}

pair < IteratorSTI < GitterPll :: hedge_STI > *, IteratorSTI < GitterPll :: hedge_STI > *> 
GitterDunePll :: borderIteratorTT (const hedge_STI * e, int link )
{
  // return edge iterator over all edges 
  is_def_true< hedge_STI > * s = 0;
  return this->createEdgeIteratorTT(s, link);
}

pair < IteratorSTI < GitterPll :: hface_STI > *, IteratorSTI < GitterPll :: hface_STI > *> 
GitterDunePll :: borderIteratorTT (const hface_STI * f, int link )
{
  // return default vertex iterator 
  is_def_true< hface_STI > * s = 0;
  return this->createFaceIteratorTT(s, link);
}

template <class ObjectStreamType, class HItemType> 
void GitterDunePll :: sendSlaves (
    ObjectStreamType & sendBuff, 
    HItemType * fakeItem ,
    GatherScatterType & dataHandle, const int link )
{
  // temporary buffer 
  SmallObjectStream osTmp; 

  pair < IteratorSTI < HItemType > *, IteratorSTI < HItemType > *> 
    a = borderIteratorTT (fakeItem, link ); //ueber alle meine Slave-Knoten 
 
  IteratorSTI < HItemType > & iter = *(a.second);
  for (iter.first (); ! iter.done () ; iter.next ()) 
  {
    HItemType & item = iter.item();
    // gather all data on slaves 
    if ( dataHandle.containsItem(item) ) 
    {
      // write marker that show data is transmitted 
      sendBuff.writeObject( transmittedData );

      osTmp.clear();
      // write data to fake buff to determine size of data package
      dataHandle.sendData(osTmp,item);

      int s = osTmp.size();
      // first write size 
      sendBuff.writeObject(s);
      // then write bytes 
      sendBuff.writeStream(osTmp);
    } 
    else 
    {
      // write noData marker 
      sendBuff.writeObject( noData );
    }
  }
  delete a.first;
  delete a.second;      

  return ;
}

template <class ObjectStreamType, class HItemType> 
void GitterDunePll :: unpackOnMaster (
    ObjectStreamType & recvBuff, 
    HItemType * determType,
    GatherScatterType & dataHandle ,
    const int nl, const int link )
{
  int hasdata;

  typedef SmallObjectStream BufferType;
  typedef vector< BufferType > DataBufferType;

  pair < IteratorSTI < HItemType > *, IteratorSTI < HItemType > *> 
    a = borderIteratorTT (determType, link);
 
  IteratorSTI < HItemType > & iter = *(a.first);

  // for all master items 
  for (iter.first (); ! iter.done () ; iter.next ()) 
  {
    HItemType & item = iter.item();
   
    // read data marker 
    recvBuff.readObject(hasdata);
    
    //DataType & data = masterData[idx];
    //int idx = vx.getIndex(); 
    //const size_t s = vertexData.size(vx);
    //if( data.size() <= nlData ) data.resize(nlData);
    
    item.reserveBuffer( nl + 1 );
    DataBufferType & data = item.commBuffer();

    if ( dataHandle.containsItem( item ) ) 
    {
      // pack master data 
      BufferType & mData = data[nl]; 
      mData.clear();
        
      // write master data to fake buffer 
      dataHandle.sendData(mData,item);
    }

    // if data has been send, read data 
    if (hasdata != noData) 
    {
      // pack slave data to tmnp buffer 
      BufferType & v = data[link]; 
      v.clear();

      int dataSize; 
      recvBuff.readObject(dataSize);
      recvBuff.readStream(v,dataSize);
    }
  }
  delete a.first;
  delete a.second;

  return ;
}

template <class ObjectStreamType, class HItemType> 
void GitterDunePll :: sendMaster (
    ObjectStreamType & sendBuff, 
    HItemType * determType,
    GatherScatterType & dataHandle ,
    const int nl , 
    const int myLink )
{
  typedef SmallObjectStream BufferType;
  typedef vector< BufferType > DataBufferType;

  pair < IteratorSTI < HItemType > *, IteratorSTI < HItemType > *> 
    a = borderIteratorTT (determType , myLink ); //ueber alle meine Slave-Knoten
 
  IteratorSTI < HItemType > & iter = *(a.first);

  // for all master items 
  for (iter.first (); ! iter.done () ; iter.next ()) 
  {
    HItemType & item = iter.item();
    DataBufferType & dataBuff = item.commBuffer();
    
    //int idx = vx.getIndex();
    //DataType & data = masterData[idx];
    
    // scatter on master 
    if ( dataHandle.containsItem( item ) ) 
    {
      for(int link = 0; link<nl; ++link)
      {
        BufferType & localBuff = dataBuff[link];
        if( localBuff.size() > 0 ) 
        {
          localBuff.resetReadPosition();
          dataHandle.recvData(localBuff, item);
        }
      }
    } 
   
    // pack for slaves 
    {
      // write data marker 
      sendBuff.writeObject(transmittedData);

      for(int link = 0; link<nl; ++link)
      {
        // if myLink == link then write master data
        // instead of data of link 
        // we do not send link i its own data
        int l = (link == myLink) ? nl : link;

        BufferType & localBuff = dataBuff[l];
        int s = localBuff.size();
        sendBuff.writeObject(s);
        // if buffer size > 0 write hole buffer to stream 
        if( s > 0 ) sendBuff.writeStream( localBuff );
      }
    } 
  }
  delete a.first;
  delete a.second;     

  return ;
}

template <class ObjectStreamType, class HItemType> 
void GitterDunePll :: unpackOnSlaves (
    ObjectStreamType & recvBuff, 
    HItemType * determType,
    GatherScatterType & dataHandle ,
    const int nOtherlinks, const int myLink)
{
  int hasdata;

  pair < IteratorSTI < HItemType > *, IteratorSTI < HItemType > *> 
    a = borderIteratorTT (determType, myLink );

  // get slave iterator 
  IteratorSTI < HItemType > & iter = *(a.second);
  
  for (iter.first (); ! iter.done () ; iter.next ()) 
  {
    // read data marker 
    recvBuff.readObject(hasdata);

    if (hasdata != noData) 
    {
      HItemType & item = iter.item();
      if( dataHandle.containsItem( item ) )
      {
        // for number of recived data, do scatter 
        for(int link = 0; link<nOtherlinks; ++link)
        {
          int s;
          recvBuff.readObject(s);
          if(s > 0) dataHandle.recvData(recvBuff, item );
        }
      }
      else 
      {
        // for number of recived data, do remove  
        for(int link = 0; link<nOtherlinks; ++link)
        {
          int s;
          recvBuff.readObject(s); 
          // if no data for link exists, s == 0
          // otherwise remove s bytes from stream by increasing 
          // read byte counter 
          if(s > 0) recvBuff.removeObject( s );
        }
      }
    }
  }
  delete a.first;
  delete a.second;
}

template <class ObjectStreamType, class HItemType> 
void GitterDunePll :: sendFaces (
    ObjectStreamType & sendBuff, 
    HItemType * fakeItem ,
    IteratorSTI < HItemType > * iter , 
    GatherScatterType & faceData )
{
  // temporary object buffer  
  SmallObjectStream osTmp; 
  
  for (iter->first () ; ! iter->done () ; iter->next ()) 
  {
    hface_STI & face = iter->item();
    if ( faceData.containsItem( face ) ) 
    {
      sendBuff.writeObject(transmittedData);
      osTmp.clear();
      faceData.sendData(osTmp, face );

      int size = osTmp.size();
      // determin size of data to be able to remove 
      sendBuff.writeObject(size);
      if( size > 0 ) sendBuff.writeStream( osTmp );
    }
    else 
    {
      sendBuff.writeObject(noData);
    }
  }
}

template <class ObjectStreamType, class HItemType> 
void GitterDunePll :: unpackFaces (
    ObjectStreamType & recvBuff, 
    HItemType * fakeItem ,
    IteratorSTI < HItemType > * iter , 
    GatherScatterType & faceData )
{
  int hasdata;
  for (iter->first () ; ! iter->done () ; iter->next ()) 
  {
    recvBuff.readObject(hasdata);
    if (hasdata != noData) 
    {
      hface_STI & face = iter->item();
      int size; 
      recvBuff.readObject(size); 
      if( size > 0 )
      {
        // if entity is not contained just remove data from stream 
        if ( faceData.containsItem( face ) ) 
          faceData.recvData(recvBuff , face );
        else 
          recvBuff.removeObject( size );
      }
    }
  }
}

////////////////////////////////////////////////////////
//
// communication of higher codim data (vertices,edges,faces)
//
////////////////////////////////////////////////////////

void GitterDunePll :: doBorderBorderComm( 
    vector< ObjectStream > & osvec ,
    GatherScatterType & vertexData , 
    GatherScatterType & edgeData,  
    GatherScatterType & faceData )
{
  const int nl = mpAccess ().nlinks ();
  
  const bool containsVertices = vertexData.contains(3,3);
  const bool containsEdges    = edgeData.contains(3,2);
  const bool containsFaces    = faceData.contains(3,1);

  const bool haveVerticesOrEdges = containsVertices || containsEdges;
   
  assert ((debugOption (5) && containsVertices) ? (cout << "**INFO GitterDunePll :: borderBorderComm (): (containsVertices)=true " << endl, 1) : 1) ;
  assert ((debugOption (5) && containsEdges)    ? (cout << "**INFO GitterDunePll :: borderBorderComm (): (containsEdges)=true " << endl, 1) : 1) ;
  assert ((debugOption (5) && containsFaces)    ? (cout << "**INFO GitterDunePll :: borderBorderComm (): (containsFaces)=true " << endl, 1) : 1) ;
   
  // buffer on master 
  //typedef vector< double > BuffType;
  //typedef SmallObjectStream BuffType;
  //typedef vector< BuffType > DataType;
  //map < int , DataType > masterData;

  {
    // gather all data from slaves 
    for (int link = 0; link < nl ; ++link )  
    {
      ObjectStream & sendBuff = osvec[link];
      sendBuff.clear();

      if (containsVertices)
      {
        vertex_STI * determType = 0;
        sendSlaves(sendBuff,determType,vertexData , link);
      }
      
      if (containsEdges) 
      {
        hedge_STI * determType = 0;
        sendSlaves(sendBuff,determType, edgeData , link);
      }
      
      if (containsFaces) 
      {
        hface_STI * determType = 0;
        pair < IteratorSTI < hface_STI > * , IteratorSTI < hface_STI > * >
          iterpair = borderIteratorTT(determType , link );
       
        // pack all faces that we are master on 
        sendFaces( sendBuff, determType, iterpair.first  , faceData ); 
        // pack also all faces that we are not master on 
        sendFaces( sendBuff, determType, iterpair.second , faceData ); 

        delete iterpair.first;
        delete iterpair.second;
      } 
    }
   
    /////////////////////////////////////////////////////
    // den anderen Partitionen die Slave-Daten senden
    /////////////////////////////////////////////////////
    osvec = mpAccess ().exchange (osvec);

    // now get all sended data and store on master item in local buffers
    for (int link = 0; link < nl; ++link) 
    { 
      ObjectStream & recvBuff = osvec[link];
      
      if (containsVertices) 
      {
        vertex_STI * determType = 0;
        unpackOnMaster(recvBuff,determType,vertexData,nl,link);
      }

      if (containsEdges) 
      {
        hedge_STI * determType = 0;
        unpackOnMaster(recvBuff,determType,edgeData,nl,link);
      }

      if (containsFaces) 
      {
        hface_STI * determType = 0;
        pair < IteratorSTI < hface_STI > * , IteratorSTI < hface_STI > * >
          iterpair = borderIteratorTT( determType , link );

        // first unpack slave data 
        unpackFaces(recvBuff,determType,iterpair.second,faceData);
        // then unpack all master data 
        unpackFaces(recvBuff,determType,iterpair.first ,faceData);

        delete iterpair.first;
        delete iterpair.second;
      }
    }
  }

  // now get all data from the local buffer of the master 
  // and send this data to the slaves (only for vertices and edges)
  if( haveVerticesOrEdges )
  {
    for (int link = 0; link < nl; ++link ) 
    {
      ObjectStream & sendBuff = osvec[link];
      sendBuff.clear();
      
      // write Number of my links 
      sendBuff.writeObject(nl); 

      if (containsVertices) 
      {
        vertex_STI * determType = 0;
        sendMaster(sendBuff,determType,vertexData,nl, link );
      }
      
      if (containsEdges) 
      {
        hedge_STI * determType = 0;
        sendMaster(sendBuff,determType,edgeData,nl, link );
      }
    }
   
    ///////////////////////////////////////////////////
    // exchange all gathered data 
    ///////////////////////////////////////////////////
    osvec = mpAccess ().exchange (osvec);
   
    // now unpack all data on slave items 
    for (int link = 0; link < nl; ++link) 
    { 
      ObjectStream & recvBuff = osvec[link];
      
      int nOtherlinks;
      recvBuff.readObject(nOtherlinks); // read number of links 

      if (containsVertices) 
      {
        vertex_STI * determType = 0;
        unpackOnSlaves(recvBuff,determType,vertexData, nOtherlinks, link );
      }
      
      if (containsEdges) 
      {
        hedge_STI * determType = 0;
        unpackOnSlaves(recvBuff,determType, edgeData, nOtherlinks, link );
      }
    }
  } // end second loop over vertices and edges 

  return ;
}

////////////////////////////////////////////////////////
//
// communicate data
// 
////////////////////////////////////////////////////////
void GitterDunePll :: ALUcomm ( 
             GatherScatterType & vertexData , 
             GatherScatterType & edgeData,  
             GatherScatterType & faceData ,
             GatherScatterType & elementData ) 
{
  const int nl = mpAccess ().nlinks ();

  const bool containsVertices = vertexData.contains(3,3);
  const bool containsEdges    = edgeData.contains(3,2);
  const bool containsFaces    = faceData.contains(3,1);
  const bool containsElements = elementData.contains(3,0);

  const bool haveHigherCodimData = containsVertices || 
    containsEdges ||  
    containsFaces ;
   
  assert ((debugOption (5) && containsVertices) ? (cout << "**INFO GitterDunePll :: ALUcomm (): (containsVertices)=true " << endl, 1) : 1) ;
  assert ((debugOption (5) && containsEdges)    ? (cout << "**INFO GitterDunePll :: ALUcomm (): (containsEdges)=true " << endl, 1) : 1) ;
  assert ((debugOption (5) && containsFaces)    ? (cout << "**INFO GitterDunePll :: ALUcomm (): (containsFaces)=true " << endl, 1) : 1) ;
  assert ((debugOption (5) && containsElements) ? (cout << "**INFO GitterDunePll :: ALUcomm (): (containsElements)=true " << endl, 1) : 1) ;
   
  typedef is_def_true < hface_STI > SendRule_t;
  typedef Insert < AccessIteratorTT < hface_STI > :: InnerHandle,
    TreeIterator < hface_STI, SendRule_t > > InnerSendIteratorType;
  typedef Insert < AccessIteratorTT < hface_STI > :: OuterHandle,
    TreeIterator < hface_STI, SendRule_t > > OuterSendIteratorType;
  typedef is_def_true < hface_STI > RecvRule_t;
  typedef Insert < AccessIteratorTT < hface_STI > :: InnerHandle,
    TreeIterator < hface_STI, RecvRule_t > > InnerRecvIteratorType;
  typedef Insert < AccessIteratorTT < hface_STI > :: OuterHandle,
    TreeIterator < hface_STI, RecvRule_t > > OuterRecvIteratorType;
  typedef IteratorSTI < hface_STI > IteratorType;

  // identifier for no data transmitted 
  const int noData = 0;

  // identifier for no data transmitted 
  const int transmittedData = 1;

  // create vector of message buffers 
  // this vector is created here, that the buffer is allocated only once 
  vector < ObjectStream > vec (nl) ;
 
  // if data on entities of higer codim exists
  // then communication if more complicated 
  if (haveHigherCodimData )
  {
    doBorderBorderComm( vec, vertexData, edgeData, faceData );
#if 0
    // buffer on master 
    //typedef vector< double > BuffType;
    //typedef SmallObjectStream BuffType;
    //typedef vector< BuffType > DataType;
    //map < int , DataType > masterData;

    {
      // gather all data from slaves 
      for (int link = 0; link < nl ; ++link )  
      {
        ObjectStream & sendBuff = vec[link];
        sendBuff.clear();

        if (containsVertices)
        {
          vertex_STI * determType = 0;
          sendSlaves(sendBuff,determType,vertexData , link);
        }
        
        if (containsEdges) 
        {
          hedge_STI * determType = 0;
          sendSlaves(sendBuff,determType, edgeData , link);
        }
        
        if (containsFaces) 
        {
          hface_STI * determType = 0;
          pair < IteratorSTI < hface_STI > * , IteratorSTI < hface_STI > * >
            iterpair = borderIteratorTT(determType , link );
         
          // pack all faces that we are master on 
          sendFaces( sendBuff, determType, iterpair.first  , faceData ); 
          // pack also all faces that we are not master on 
          sendFaces( sendBuff, determType, iterpair.second , faceData ); 

          delete iterpair.first;
          delete iterpair.second;
        } 
      }
     
      /////////////////////////////////////////////////////
      // den anderen Partitionen die Slave-Daten senden
      /////////////////////////////////////////////////////
      vec = mpAccess ().exchange (vec);

      // now get all sended data and store on master item in local buffers
      for (int link = 0; link < nl; ++link) 
      { 
        ObjectStream & recvBuff = vec[link];
        
        if (containsVertices) 
        {
          vertex_STI * determType = 0;
          unpackOnMaster(recvBuff,determType,vertexData,nl,link);
        }

        if (containsEdges) 
        {
          hedge_STI * determType = 0;
          unpackOnMaster(recvBuff,determType,edgeData,nl,link);
        }

        if (containsFaces) 
        {
          hface_STI * determType = 0;
          pair < IteratorSTI < hface_STI > * , IteratorSTI < hface_STI > * >
            iterpair = borderIteratorTT( determType , link );

          // first unpack slave data 
          unpackFaces(recvBuff,determType,iterpair.second,faceData);
          // then unpack all master data 
          unpackFaces(recvBuff,determType,iterpair.first ,faceData);

          delete iterpair.first;
          delete iterpair.second;
        }
      }
    }

    // now get all data from the local buffer of the master 
    // and send this data to the slaves (only for vertices and edges)
    if( containsVertices || containsEdges )
    {
      for (int link = 0; link < nl; ++link ) 
      {
        ObjectStream & sendBuff = vec[link];
        sendBuff.clear();
        
        // write Number of my links 
        sendBuff.writeObject(nl); 

        if (containsVertices) 
        {
          vertex_STI * determType = 0;
          sendMaster(sendBuff,determType,vertexData,nl, link );
        }
        
        if (containsEdges) 
        {
          hedge_STI * determType = 0;
          sendMaster(sendBuff,determType,edgeData,nl, link );
        }
      }
     
      ///////////////////////////////////////////////////
      // exchange all gathered data 
      ///////////////////////////////////////////////////
      vec = mpAccess ().exchange (vec);
     
      // now unpack all data on slave items 
      for (int link = 0; link < nl; ++link) 
      { 
        ObjectStream & recvBuff = vec[link];
        
        int nOtherlinks;
        recvBuff.readObject(nOtherlinks); // read number of links 

        if (containsVertices) 
        {
          vertex_STI * determType = 0;
          unpackOnSlaves(recvBuff,determType,vertexData, nOtherlinks, link );
        }
        
        if (containsEdges) 
        {
          hedge_STI * determType = 0;
          unpackOnSlaves(recvBuff,determType, edgeData, nOtherlinks, link );
        }
      }
    } // end second loop over vertices and edges 
#endif  
    
  } // end haveHigherCodimData 

    {
      for (int l = 0 ; l < nl ; l ++) 
      {   
        vec[l].clear();
        {
          AccessIteratorTT < hface_STI > :: InnerHandle mif1_(containerPll (), l);
          InnerSendIteratorType w(mif1_);
          for (w.first () ; ! w.done () ; w.next ()) 
          {
            pair < ElementPllXIF_t *, int > p = 
              w.item ().accessPllX ().accessInnerPllX () ;

            if (w.item().isInteriorLeaf()) 
            {
              vec[l].writeObject(transmittedData);
              if (containsVertices) p.first->VertexData2os(vec[l], vertexData);
              if (containsEdges)    p.first->EdgeData2os(vec[l], edgeData);
              if (containsFaces)    p.first->FaceData2os(vec[l], faceData);
              if (containsElements) 
              {
                p.first->writeDynamicState (vec [l], elementData) ;
              }
            }     
            else vec[l].writeObject(noData);
          }
        }
        
      {
        AccessIteratorTT < hface_STI > :: OuterHandle mif1_(containerPll (), l);
        OuterSendIteratorType w(mif1_);
        for (w.first () ; ! w.done () ; w.next ()) 
        {
          pair < ElementPllXIF_t *, int > p = 
            w.item ().accessPllX ().accessInnerPllX () ;

          if (w.item().isInteriorLeaf()) 
          {
            vec[l].writeObject(transmittedData);
            if (containsVertices) p.first->VertexData2os(vec[l], vertexData);
            if (containsEdges)    p.first->EdgeData2os(vec[l], edgeData);
            if (containsFaces)    p.first->FaceData2os(vec[l], faceData);
            if (containsElements) 
              p.first->writeDynamicState (vec [l], elementData) ;
          }
          else vec[l].writeObject(noData);
        }
      }
    }

    ///////////////////////////////////////////
    // exchange data 
    ///////////////////////////////////////////
    vec = mpAccess ().exchange (vec) ;     
    
    //all ghost cells get new data
    for (int l = 0 ; l < nl ; l ++ ) 
    {  
      {
        AccessIteratorTT < hface_STI > :: OuterHandle mif1_(containerPll (), l);
        OuterRecvIteratorType w(mif1_);

        int hasdata;        
        for (w.first () ; ! w.done () ; w.next ()) 
        {
          pair < ElementPllXIF_t *, int > p = 
            w.item ().accessPllX ().accessOuterPllX () ;
          vec[l].readObject(hasdata);

          if (hasdata != noData) 
          {
            assert(p.first->checkGhostLevel());
            assert(p.first->ghostLeaf());

            Gitter :: helement_STI * ghost = p.first->getGhost().first;
            assert( ghost );
          
            if (containsVertices) 
              ghost->os2VertexData(vec[l], vertexData);
            if (containsEdges)    
              ghost->os2EdgeData(vec[l], edgeData);
            if (containsFaces)    
              ghost->os2FaceData(vec[l], faceData);
            if (containsElements) 
              p.first->readDynamicState (vec [l], elementData);
          }
        }
      }


      { 
        AccessIteratorTT < hface_STI > :: InnerHandle mif1_(containerPll (), l);
        InnerRecvIteratorType w(mif1_);
        int hasdata;        
        for (w.first () ; ! w.done () ; w.next ()) 
        {
          pair < ElementPllXIF_t *, int > p = 
            w.item ().accessPllX ().accessOuterPllX () ;
        
          vec[l].readObject(hasdata);
          if (hasdata != noData) 
          {
            assert(p.first->checkGhostLevel());
            assert(p.first->ghostLeaf());
          
            Gitter :: helement_STI * ghost = p.first->getGhost().first;
            assert( ghost );
          
            if (containsVertices) 
              ghost->os2VertexData(vec[l], vertexData);
            if (containsEdges)    
              ghost->os2EdgeData(vec[l], edgeData);
            if (containsFaces)    
              ghost->os2FaceData(vec[l], faceData);
            if (containsElements) 
              p.first->readDynamicState (vec [l], elementData) ;
          }
        }
      }
    }
  } // end element communication 
}

bool GitterDunePll :: refine () {
  assert (debugOption (5) ? (cout << "**INFO GitterDunePll :: refine () " << endl, 1) : 1) ;
  const int nl = mpAccess ().nlinks () ;
  bool state = false ;
  vector < vector < hedge_STI * > > innerEdges (nl), outerEdges (nl) ;
  vector < vector < hface_STI * > > innerFaces (nl), outerFaces (nl) ;
  {
    // Erst die Zeiger auf alle Fl"achen und Kanten mit paralleler
    // Mehrdeutigkeit sichern, da die LeafIteratorTT < . > nach dem 
    // Verfeinern auf gitter nicht mehr stimmen werden. Die Technik
    // ist zul"assig, da keine mehrfache Verfeinerung entstehen kann.
  
    {for (int l = 0 ; l < nl ; l ++) {
  //cout << "refinepll \n";
  LeafIteratorTT < hface_STI > fw (*this,l) ;
  LeafIteratorTT < hedge_STI > dw (*this,l) ;
  for (fw.outer ().first () ; ! fw.outer().done () ; fw.outer ().next ())
    outerFaces [l].push_back (& fw.outer ().item ()) ;
  for (fw.inner ().first () ; ! fw.inner ().done () ; fw.inner ().next ())
    innerFaces [l].push_back (& fw.inner ().item ()) ;
  for (dw.outer ().first () ; ! dw.outer().done () ; dw.outer ().next ())
    outerEdges [l].push_back (& dw.outer ().item ()) ;
  for (dw.inner ().first () ; ! dw.inner ().done () ; dw.inner ().next ())
    innerEdges [l].push_back (& dw.inner ().item ()) ;
      }}
    // jetzt normal verfeinern und den Status der Verfeinerung
    // [unvollst"andige / vollst"andige Verfeinerung] sichern.
    
    __STATIC_phase = 1 ;
    
    state = Gitter :: refine () ;
       
    // Phase des Fl"achenausgleichs an den Schnittfl"achen des
    // verteilten Gitters. Weil dort im sequentiellen Fall pseudorekursive
    // Methodenaufrufe vorliegen k"onnen, muss solange iteriert werden,
    // bis die Situation global station"ar ist.
  
    __STATIC_phase = 2 ;
  
    bool repeat (false) ;
    _refineLoops = 0 ;
    do {
      repeat = false ;
      {
        vector < ObjectStream > osv (nl) ;
        try {
    //cout << "refinepll 2 \n";
    for (int l = 0 ; l < nl ; l ++) {
            {for (vector < hface_STI * > :: const_iterator i = outerFaces [l].begin () ;
      i != outerFaces [l].end () ; (*i ++)->accessPllX ().accessOuterPllX ().first->getRefinementRequest (osv [l])) ; }
            {for (vector < hface_STI * > :: const_iterator i = innerFaces [l].begin () ;
      i != innerFaces [l].end () ; (*i ++)->accessPllX ().accessOuterPllX ().first->getRefinementRequest (osv [l])) ; }
          }
  } catch (Parallel :: AccessPllException) {
          cerr << "**FEHLER (FATAL) AccessPllException in " << __FILE__ << " " << __LINE__ << endl ; abort () ;
        }
  
        osv = mpAccess ().exchange (osv) ;
  
        try {
    for (int l = 0 ; l < nl ; l ++) {
            {for (vector < hface_STI * > :: const_iterator i = innerFaces [l].begin () ;
      i != innerFaces [l].end () ; repeat |= (*i ++)->accessPllX ().accessOuterPllX ().first->setRefinementRequest (osv [l])) ; }
            {for (vector < hface_STI * > :: const_iterator i = outerFaces [l].begin () ;
      i != outerFaces [l].end () ; repeat |= (*i ++)->accessPllX ().accessOuterPllX ().first->setRefinementRequest (osv [l])) ; }
          }
  } catch (Parallel :: AccessPllException) {
          cerr << "**FEHLER (FATAL) AccessPllException in " << __FILE__ << " " << __LINE__ << endl ; abort () ;
        }
      }

      _refineLoops ++ ;
    } while (mpAccess ().gmax (repeat ? 1 : 0)) ;

    // Jetzt noch die Kantensituation richtigstellen, es gen"ugt ein Durchlauf,
    // weil die Verfeinerung einer Kante keine Fernwirkungen hat. Vorsicht: Die
    // Kanten sind bez"uglich ihrer Identifikation sternf"ormig organisiert, d.h.
    // es muss die Verfeinerungsinformation einmal am Eigent"umer gesammelt und
    // dann wieder zur"ucktransportiert werden, eine einfache L"osung, wie bei
    // den Fl"achen (1/1 Beziehung) scheidet aus.

    __STATIC_phase = 3 ;

    {
      vector < ObjectStream > osv (nl) ;
      {for (int l = 0 ; l < nl ; l ++) 
    for (vector < hedge_STI * > :: const_iterator i = outerEdges [l].begin () ;
         i != outerEdges [l].end () ; (*i ++)->accessPllX ().getRefinementRequest (osv [l])) ;
      }
      osv = mpAccess ().exchange (osv) ;
      {for (int l = 0 ; l < nl ; l ++)
    for (vector < hedge_STI * > :: const_iterator i = innerEdges [l].begin () ;
         i != innerEdges [l].end () ; (*i ++)->accessPllX ().setRefinementRequest (osv [l])) ;
      }
    }   // ~vector < ObjectStream > ... 
    {
      vector < ObjectStream > osv (nl) ;
      {for (int l = 0 ; l < nl ; l ++)
    for (vector < hedge_STI * > :: const_iterator i = innerEdges [l].begin () ;
         i != innerEdges [l].end () ; (*i ++)->accessPllX ().getRefinementRequest (osv [l])) ;
      }
      osv = mpAccess ().exchange (osv) ;
      {for (int l = 0 ; l < nl ; l ++)
    for (vector < hedge_STI * > :: const_iterator i = outerEdges [l].begin () ;
         i != outerEdges [l].end () ; (*i ++)->accessPllX ().setRefinementRequest (osv [l])) ;
      }
    }   // ~vector < ObjectStream > ... 
  }
  
  __STATIC_phase = -1 ;
  
  return state ;
}

void GitterDunePll :: coarse () {
  assert (debugOption (20) ? (cout << "**INFO GitterDunePll :: coarse () " << endl, 1) : 1) ;
  const int nl = mpAccess ().nlinks () ;
  
  {
    vector < vector < hedge_STI * > > innerEdges (nl), outerEdges (nl) ;
    vector < vector < hface_STI * > > innerFaces (nl), outerFaces (nl) ;
  
    for (int l = 0 ; l < nl ; l ++) {
    
      // Zun"achst werden f"ur alle Links die Zeiger auf Gitterojekte mit
      // Mehrdeutigkeit gesichert, die an der Wurzel einer potentiellen
      // Vergr"oberungsoperation sitzen -> es sind die Knoten in der Hierarchie,
      // deren Kinder alle Bl"atter sind. Genau diese Knoten sollen gegen"uber
      // der Vergr"oberung blockiert werden und dann die Vergr"oberung falls
      // sie zul"assig ist, sp"ater durchgef"uhrt werden (pending) ;
    
      AccessIteratorTT < hface_STI > :: InnerHandle mfwi (containerPll (),l) ;
      AccessIteratorTT < hface_STI > :: OuterHandle mfwo (containerPll (),l) ;
      AccessIteratorTT < hedge_STI > :: InnerHandle mdwi (containerPll (),l) ;
      AccessIteratorTT < hedge_STI > :: OuterHandle mdwo (containerPll (),l) ;
      
      // Die inneren und a"usseren Iteratoren der potentiell vergr"oberungsf"ahigen
      // Fl"achen "uber den Grobgitterfl"achen. In den Elementen passiert erstmal
      // nichts, solange nicht mit mehrfachen Grobgitterelementen gearbeitet wird.
      
      Insert < AccessIteratorTT < hface_STI > :: InnerHandle, 
        TreeIterator < hface_STI, childs_are_leafs < hface_STI > > > fwi (mfwi) ;
      Insert < AccessIteratorTT < hface_STI > :: OuterHandle, 
        TreeIterator < hface_STI, childs_are_leafs < hface_STI > > > fwo (mfwo) ;
      
      // Die inneren und a"usseren Iteratoren der potentiell vergr"oberungsf"ahigen
      // Kanten "uber den Grobgitterkanten.
      
      Insert < AccessIteratorTT < hedge_STI > :: InnerHandle, 
        TreeIterator < hedge_STI, childs_are_leafs < hedge_STI > > > dwi (mdwi) ;
      Insert < AccessIteratorTT < hedge_STI > :: OuterHandle, 
        TreeIterator < hedge_STI, childs_are_leafs < hedge_STI > > > dwo (mdwo) ;

      // Die inneren und a"usseren Iteratoren der potentiell vergr"oberungsf"ahigen
      // Kanten "uber den Grobgitterfl"achen. Diese Konstruktion wird beim Tetraeder-
      // gitter notwendig, weil dort keine Aussage der Form:
      //

      Insert < AccessIteratorTT < hface_STI > :: InnerHandle, 
        TreeIterator < hface_STI, has_int_edge < hface_STI > > > efi (mfwi) ;
      Insert < AccessIteratorTT < hface_STI > :: OuterHandle, 
        TreeIterator < hface_STI, has_int_edge < hface_STI > > > efo (mfwo) ;
      Wrapper < Insert < AccessIteratorTT < hface_STI > :: InnerHandle, 
        TreeIterator < hface_STI, has_int_edge < hface_STI > > >, InternalEdge > eifi (efi) ;
      Wrapper < Insert < AccessIteratorTT < hface_STI > :: OuterHandle, 
        TreeIterator < hface_STI, has_int_edge < hface_STI > > >, InternalEdge > eifo (efo) ;
      Insert < Wrapper < Insert < AccessIteratorTT < hface_STI > :: InnerHandle, 
        TreeIterator < hface_STI, has_int_edge < hface_STI > > >, InternalEdge >,
  TreeIterator < hedge_STI, childs_are_leafs < hedge_STI > > > dfi (eifi) ;
      Insert < Wrapper < Insert < AccessIteratorTT < hface_STI > :: OuterHandle, 
        TreeIterator < hface_STI, has_int_edge < hface_STI > > >, InternalEdge >,
  TreeIterator < hedge_STI, childs_are_leafs < hedge_STI > > > dfo (eifo) ;

      // Die 'item ()' Resultatwerte (Zeiger) werden in Vektoren gesichert, weil die
      // Kriterien die zur Erzeugung der Iteratoren angewendet wurden (Filter) nach
      // einer teilweisen Vergr"oberung nicht mehr g"ultig sein werden, d.h. die 
      // Iterationsobjekte "andern w"ahrend der Vergr"oberung ihre Eigenschaften.
      // Deshalb werden sie auch am Ende des Blocks aufgegeben. Der Vektor 'cache'
      // ist zul"assig, weil kein Objekt auf das eine Referenz im 'cache' vorliegt
      // beseitigt werden kann. Sie sind alle ein Niveau darunter.

      for (fwi.first () ; ! fwi.done () ; fwi.next ()) innerFaces [l].push_back (& fwi.item ()) ;
      for (fwo.first () ; ! fwo.done () ; fwo.next ()) outerFaces [l].push_back (& fwo.item ()) ;
      for (dwo.first () ; ! dwo.done () ; dwo.next ()) outerEdges [l].push_back (& dwo.item ()) ;
      for (dfo.first () ; ! dfo.done () ; dfo.next ()) outerEdges [l].push_back (& dfo.item ()) ;
      for (dwi.first () ; ! dwi.done () ; dwi.next ()) innerEdges [l].push_back (& dwi.item ()) ;
      for (dfi.first () ; ! dfi.done () ; dfi.next ()) innerEdges [l].push_back (& dfi.item ()) ;
    }
    try {
      // Erstmal alles was mehrdeutig ist, gegen die drohende Vergr"oberung sichern.
      // Danach werden sukzessive die Fl"achenlocks aufgehoben, getestet und
      // eventuell vergr"obert, dann das gleiche Spiel mit den Kanten.

      for (int l = 0 ; l < nl ; l ++) {
        {for (vector < hedge_STI * > :: iterator i = outerEdges [l].begin () ;
        i != outerEdges [l].end () ; (*i ++)->accessPllX ().lockAndTry ()) ; }
        {for (vector < hedge_STI * > :: iterator i = innerEdges [l].begin () ;
        i != innerEdges [l].end () ; (*i ++)->accessPllX ().lockAndTry ()) ; }
        {for (vector < hface_STI * > :: iterator i = outerFaces [l].begin () ;
        i != outerFaces [l].end () ; (*i ++)->accessPllX ().accessOuterPllX ().first->lockAndTry ()) ; }
        {for (vector < hface_STI * > :: iterator i = innerFaces [l].begin () ;
        i != innerFaces [l].end () ; (*i ++)->accessPllX ().accessOuterPllX ().first->lockAndTry ()) ; }
      }
      
      // Gitter :: coarse () ist elementorientiert, d.h. die Vergr"oberung auf Fl"achen und
      // Kanten wird nur durch Vermittlung eines sich vergr"obernden Knotens in der Element-
      // hierarchie angestossen. In allen gegen Vergr"oberung 'gelockten' Fl"achen und Kanten
      // wird die angeforderte Operation zur"uckgewiesen, um erst sp"ater von aussen nochmals
      // angestossen zu werden.
      
      __STATIC_phase = 4 ;
      
      Gitter :: coarse () ;
      
    } catch (Parallel :: AccessPllException) {
      cerr << "**FEHLER (FATAL) AccessPllException beim Vergr\"obern der Elementhierarchie oder\n" ;
      cerr << "  beim locken der Fl\"achen- bzw. Kantenb\"aume aufgetreten. In " << __FILE__ << " " << __LINE__ << endl ;
      abort () ;
    }
    try {
    
      // Phase des Fl"achenausgleichs des verteilten Vergr"oberungsalgorithmus
      // alle Schnittfl"achenpaare werden daraufhin untersucht, ob eine
      // Vergr"oberung in beiden Teilgittern durchgef"uhrt werden darf,
      // wenn ja, wird in beiden Teilgittern vergr"obert und der Vollzug
      // getestet.
  
      __STATIC_phase = 5 ;
    
      vector < vector < int > > clean (nl) ;
      {
        vector < vector < int > > inout (nl) ;
        {for (int l = 0 ; l < nl ; l ++)
      for (vector < hface_STI * > :: iterator i = outerFaces [l].begin () ; i != outerFaces [l].end () ; i ++)
        inout [l].push_back ((*i)->accessPllX ().accessOuterPllX ().first->lockAndTry ()) ;
        }
        inout = mpAccess ().exchange (inout) ;
        {for (int l = 0 ; l < nl ; l ++) {
      clean [l] = vector < int > (innerFaces [l].size (), long (true)) ;
      vector < int > :: iterator j = clean [l].begin (), k = inout [l].begin () ;
      for (vector < hface_STI * > :: iterator i = innerFaces [l].begin () ; i != innerFaces [l].end () ; i ++, j++, k++) {
        assert (j != clean [l].end ()) ; assert (k != inout [l].end ()) ;
        (*j) &= (*k) && (*i)->accessPllX ().accessOuterPllX ().first->lockAndTry () ;
      }
    }}
      }
      {
        vector < vector < int > > inout (nl) ;
        {for (int l = 0 ; l < nl ; l ++) {
      vector < int > :: iterator j = clean [l].begin () ;
      for (vector < hface_STI * > :: iterator i = innerFaces [l].begin () ; i != innerFaces [l].end () ; i ++, j++) {
        inout [l].push_back (*j) ;
        (*i)->accessPllX ().accessOuterPllX ().first->unlockAndResume (bool (*j)) ;
      }
    }}
      
        inout = mpAccess ().exchange (inout) ;
      
        {for (int l = 0 ; l < nl ; l ++) {
      vector < int > :: iterator j = inout [l].begin () ;
      for (vector < hface_STI * > :: iterator i = outerFaces [l].begin () ; i != outerFaces [l].end () ; i ++, j++) {
        assert (j != inout [l].end ()) ;
        (*i)->accessPllX ().accessOuterPllX ().first->unlockAndResume (bool (*j)) ;
      }
    }}
      }
    } catch (Parallel :: AccessPllException) {
      cerr << "**FEHLER (FATAL) AccessPllException beim Vergr\"obern der Fl\"achenb\"aume\n" ;
      cerr << "  aufgetreten. In " << __FILE__ << " " << __LINE__ << endl ;
      abort () ;
    }
    try {
    
      // Phase des Kantenausgleichs im parallelen Vergr"oberungsalgorithmus:
  
      __STATIC_phase  = 6 ;
    
      // Weil hier jede Kante nur eindeutig auftreten darf, muss sie in einem
      // map als Adresse hinterlegt werden, dann k"onnen die verschiedenen
      // Refcounts aus den verschiedenen Links tats"achlich global miteinander
      // abgemischt werden. Dazu werden zun"achst alle eigenen Kanten auf ihre
      // Vergr"oberbarkeit hin untersucht und dieser Zustand (true = vergr"oberbar
      // false = darf nicht vergr"obert werden) im map 'clean' hinterlegt. Dazu
      // kommt noch ein zweiter 'bool' Wert, der anzeigt ob die Kante schon ab-
      // schliessend vergr"obert wurde oder nicht. 
    
      map < hedge_STI *, pair < bool, bool >, less < hedge_STI * > > clean ;
      
      {for (int l = 0 ; l < nl ; l ++)
    for (vector < hedge_STI * > :: iterator i = innerEdges [l].begin () ; i != innerEdges [l].end () ; i ++)
      if (clean.find (*i) == clean.end ()) clean [*i] = pair < bool, bool > ((*i)->accessPllX ().lockAndTry (), true) ;
      }
      {
        vector < vector < int > > inout (nl) ;
        {for (int l = 0 ; l < nl ; l ++)
      for (vector < hedge_STI * > :: iterator i = outerEdges [l].begin () ; i != outerEdges [l].end () ; i ++)
        inout [l].push_back ((*i)->accessPllX ().lockAndTry ()) ;
  }
        inout = mpAccess ().exchange (inout) ;
        {for (int l = 0 ; l < nl ; l ++) {
      vector < int > :: const_iterator j = inout [l].begin () ;
      for (vector < hedge_STI * > :: iterator i = innerEdges [l].begin () ; i != innerEdges [l].end () ; i ++, j++) {
        assert (j != inout [l].end ()) ;
        assert (clean.find (*i) != clean.end ()) ;
        if (*j == false) clean [*i] = pair < bool, bool > (false, clean[*i].second) ; 
      }
    }}
      }
      {
        vector < vector < int > > inout (nl) ;
        {for (int l = 0 ; l < nl ; l ++) {
      for (vector < hedge_STI * > :: iterator i = innerEdges [l].begin () ; i != innerEdges [l].end () ; i ++) {
        assert (clean.find (*i) != clean.end ()) ;
        pair < bool, bool > & a = clean [*i] ;
        inout [l].push_back (a.first) ;
        if (a.second) {
  
    // Wenn wir hier sind, kann die Kante tats"achlich vergr"obert werden, genauer gesagt,
    // sie wird es auch und der R"uckgabewert testet den Vollzug der Aktion. Weil aber nur
    // einmal vergr"obert werden kann, und die Iteratoren 'innerEdges [l]' aber eventuell
    // mehrfach "uber eine Kante hinweglaufen, muss diese Vergr"oberung im map 'clean'
    // vermerkt werden. Dann wird kein zweiter Versuch unternommen.
  
    a.second = false ;
#ifndef NDEBUG
    bool b = 
#endif
      (*i)->accessPllX ().unlockAndResume (a.first) ;
    assert (b == a.first) ;
        }
      }
    }}
        inout = mpAccess ().exchange (inout) ;
        {for (int l = 0 ; l < nl ; l ++) {
      vector < int > :: iterator j = inout [l].begin () ;
      for (vector < hedge_STI * > :: iterator i = outerEdges [l].begin () ; i != outerEdges [l].end () ; i ++, j++) {
        assert (j != inout [l].end ()) ;
      
        // Selbe Situation wie oben, aber der Eigent"umer der Kante hat mitgeteilt, dass sie
        // vergr"obert werden darf und auch wird auf allen Teilgebieten also auch hier. Der
        // Vollzug der Vergr"oberung wird durch den R"uckgabewert getestet.
      
#ifndef NDEBUG
        bool b = 
#endif
    (*i)->accessPllX ().unlockAndResume (bool (*j)) ;
        assert (b == bool (*j)) ;
      }
    }}
      }
    } catch (Parallel :: AccessPllException) {
      cerr << "**FEHLER (FATAL) AccessPllException beim Vergr\"obern der Kantenb\"aume\n" ;
      cerr << "  aufgetreten. In " << __FILE__ << " " << __LINE__ << endl ;
      abort () ;
    }
  }
  
  __STATIC_phase = -1 ;
  
  return ;
}

#endif
