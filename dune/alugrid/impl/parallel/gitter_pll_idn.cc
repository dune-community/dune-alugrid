// (c) bernhard schupp 1997 - 1998
#include <config.h>

#include <list>
#include <map>
#include <set>
#include <vector>

#include "gitter_pll_sti.h"

namespace ALUGrid
{

  // vorsicht mit optimierenden compilern bei dem untenstehnden template !

  template < class A > void identify (typename AccessIterator < A >::Handle mi, 
                                      std::vector< std::pair< std::list< typename AccessIterator < A >::Handle >, 
                                                      std::list< typename AccessIterator < A >::Handle > > > & tt, 
                                      const MpAccessLocal & c) 
  {
    typedef std::set< std::vector< int > > lp_map_t;

    typedef std::map< typename LinkedObject::Identifier, 
                  std::pair< typename AccessIterator < A >::Handle, 
                  typename lp_map_t::const_iterator > > lmap_t;

    const int me = c.myrank (), nl = c.nlinks ();
    
    lp_map_t linkagePatternMap;
    lmap_t look;
    
    {
      std::vector< std::vector< int > > inout (nl);
      lp_map_t::const_iterator meIt = linkagePatternMap.insert (std::vector< int >  (1L, me)).first;
      {
        std::vector< int > count (nl);
        for(int k=0; k<nl; ++k) count[k] = 0;
        
        typename AccessIterator < A >::Handle micopy ( mi );
        for (micopy.first (); ! micopy.done (); micopy.next ()) 
        {
          std::vector< int > estimate = micopy.item ().accessPllX ().estimateLinkage ();
          if ( ! estimate.empty() ) 
          {
            std::vector< int >::const_iterator iEnd =  estimate.end ();
            for (std::vector< int >::const_iterator i = estimate.begin (); 
                 i != iEnd; ++i )
            {
              count[c.link (*i)] += 4;
            }
          }
        }
        // reserve memory 
        for(int k=0; k<nl; ++k) inout[ k ].reserve( count[k] );
      }
      
      {
        for (mi.first (); ! mi.done (); mi.next ()) 
        {
          std::vector< int > estimate = mi.item ().accessPllX ().estimateLinkage ();
          if (estimate.size ()) 
          {
            LinkedObject::Identifier id = mi.item ().accessPllX ().getIdentifier ();
            look [id].first = mi;
            look [id].second = meIt;
                  {
              std::vector< int >::const_iterator iEnd = estimate.end ();
              for (std::vector< int >::const_iterator i = estimate.begin (); 
                   i != iEnd; ++i )
                id.write (inout [c.link (*i)]);
            }
          }
        }
      }
      
      // exchange data 
      inout = c.exchange (inout);

      std::vector< int > d = c.dest ();
      { 
        for (int l = 0; l < nl; ++l ) 
        {
          std::vector< int >::const_iterator pos = inout [l].begin (), end = inout [l].end ();
          while (pos != end) 
          {
            typename LinkedObject::Identifier id;
            id.read (pos,end);
            typename lmap_t::iterator hit = look.find (id);
            if (hit != look.end ()) 
            {
              std::vector< int > lpn (*(*hit).second.second);
              if (find (lpn.begin (), lpn.end (), d [l]) == lpn.end ()) 
              {
                      lpn.push_back (d [l]);
                      std::sort (lpn.begin(), lpn.end() );
                      (*hit).second.second = linkagePatternMap.insert (lpn).first;
                    }
            }
          }
        } 
      }
    }

    tt = std::vector< std::pair< std::list< typename AccessIterator < A >::Handle >, 
                         std::list< typename AccessIterator < A >::Handle > > > (nl);
    {
      std::vector< std::vector< int > > inout (nl);
      {
        std::vector<int> count(nl);
        for(int k=0; k<nl; k++) count[k] = 0;

        for (typename lmap_t::const_iterator pos = look.begin (); 
          pos != look.end (); pos ++) 
        {
          const std::vector< int > & lk (*(*pos).second.second);
          if (* lk.begin () == me) 
          {
            typename std::vector< int >::const_iterator iEnd = lk.end ();
            for (typename std::vector< int >::const_iterator i = lk.begin (); 
                 i != iEnd; i ++) 
            {
              if (*i != me) 
              {
                int l = c.link (*i);
                count[l] += 4; 
              }
            }
          }
        }
        // reserve memory 
        for(int k=0; k<nl; ++k) inout[k].reserve( count[k] );
      }
      
      {
        const typename lmap_t::const_iterator lookEnd = look.end ();
        for (typename lmap_t::const_iterator pos = look.begin (); 
             pos != lookEnd; ++pos) 
        {
          const std::vector< int > & lk (*(*pos).second.second);
          if (* lk.begin () == me) 
          {
            typename LinkedObject::Identifier id = (*pos).second.first.item ().accessPllX ().getIdentifier ();
            { 
              typename std::vector< int >::const_iterator iEnd = lk.end ();
              for (typename std::vector< int >::const_iterator i = lk.begin (); 
                   i != iEnd; ++i) 
              {
                if (*i != me) 
                {
                  int l = c.link (*i);
                  tt [l].first.push_back ((*pos).second.first);
                  id.write (inout [l]);
                }
              } 
            }
          }
        }
      }

      // exchange data 
      inout = c.exchange (inout);
      
      {
        for (int i = 0; i < nl; ++i ) 
        {
          typename std::vector< int >::const_iterator pos = inout [i].begin (), end = inout [i].end ();
          while (pos != end) 
          {
            typename LinkedObject::Identifier id;
            id.read (pos,end);
            assert (look.find (id) != look.end ());
            tt [i].second.push_back ((*look.find (id)).second.first);
          } 
        }
      }
    }
    return;
  }

  std::set< int > GitterPll::MacroGitterPll::secondScan () 
  {
    std::set< int > s;
    {
      AccessIterator < vertex_STI >::Handle w (*this);
      for ( w.first (); ! w.done (); w.next ()) 
      {
        vertex_STI& vertex = w.item();
        // only border vertices can have linkage 
        if( vertex.isBorder() ) 
        {
          const std::vector< int > l = w.item ().accessPllX ().estimateLinkage ();
          const std::vector< int >::const_iterator iEnd = l.end ();
          for (std::vector< int >::const_iterator i = l.begin (); i != iEnd; ++i ) 
          {
            s.insert ( *i );
          }
        }
      }
    }
    return s;
  }

  void GitterPll::MacroGitterPll::vertexLinkageEstimateGCollect (MpAccessLocal & mpAccess) 
  {
    typedef std::map< int, AccessIterator < vertex_STI >::Handle > map_t;
    map_t vxmap;
    const int np = mpAccess.psize (), me = mpAccess.myrank ();

    ObjectStream os;
    {
      ObjectStream ostmp;
      AccessIterator < vertex_STI >::Handle w (*this);

      const int estimate = 0.25 * w.size() + 1;
      // reserve memory 
      ostmp.reserve( estimate * sizeof(int) );
      int size = 0;
      for (w.first (); ! w.done (); w.next ()) 
      {
        vertex_STI& vertex = w.item();

        // only insert border vertices 
        if( vertex.isBorder() )
        {
          ++size;
          int id = vertex.ident ();
          ostmp.writeObject( id );
          vxmap[ id ] = w;
        }
      }
      os.writeObject( size );
      os.writeStream( ostmp );
    }

    // exchange data 
    std::vector< ObjectStream > osv = mpAccess.gcollect( os );

    // free memory 
    os.reset();

    {
      map_t::const_iterator vxmapEnd = vxmap.end();
      for (int i = 0; i < np; ++i ) 
      {
        if (i != me) 
        {
          ObjectStream& osv_i = osv[ i ]; 

          int size;
          osv_i.readObject ( size );
          for (int j = 0; j < size; ++j ) 
          {
            int id;
            osv_i.readObject( id );
            map_t::const_iterator hit = vxmap.find (id);
            if( hit != vxmapEnd ) 
            {
              std::vector< int > s = (*hit).second.item ().accessPllX ().estimateLinkage ();
                    if (find (s.begin (), s.end (), i) == s.end ()) 
              {
                s.push_back( i );
                      (*hit).second.item ().accessPllX ().setLinkage (s);
                    }
                  }
          }
          // free memory 
          osv_i.reset();
        }
      }
    }
    return;
  }

  void GitterPll::MacroGitterPll::vertexLinkageEstimateBcast (MpAccessLocal & mpAccess) 
  {
    typedef std::map< int, vertex_STI* > map_t;
    map_t vxmap;
    const int np = mpAccess.psize ();
    const int me = mpAccess.myrank ();

    // exchange data 
    {
      AccessIterator < vertex_STI >::Handle w (*this);
      for (w.first (); ! w.done (); w.next ()) 
      {
        vertex_STI& vertex = w.item();

        // only insert border vertices 
        if( vertex.isBorder() )
        {
          vxmap[ vertex.ident() ] = &vertex;
        }
      }
    }

    // get max size for all ranks needed in later bcast cycle 
    const int mysize  = vxmap.size();
    const int maxSize = mpAccess.gmax( mysize ) + 1; // +1 for the size 

    // create buffer 
    std::vector< int > borderIds( maxSize, 0 ); 

    // loop over all ranks 
    for (int rank = 0; rank < np; ++rank ) 
    {
      const map_t::const_iterator vxmapEnd = vxmap.end();
      // fill buffer  
      if( rank == me ) 
      {
        // store the size on first position 
        borderIds[ 0 ] = mysize;
        int count = 1;
        for(map_t::const_iterator vx = vxmap.begin(); vx != vxmapEnd; ++vx, ++count )
        {
          borderIds[ count ] = (*vx).first;
        }
      }

      // send border vertex list to all others 
      mpAccess.bcast( &borderIds[ 0 ], maxSize, rank );

      // check connectivy for receives vertex ids 
      if( rank != me ) 
      {
        // get size which is the first entry 
        const int idSize = borderIds[ 0 ];
        // loop over all sended vertex ids 
        for (int j = 1; j <= idSize; ++j ) 
        {
          const int id = borderIds[ j ];
          map_t::const_iterator hit = vxmap.find (id);
          if( hit != vxmapEnd ) 
          {
            std::vector< int > s = (*hit).second->accessPllX ().estimateLinkage ();
            if (find (s.begin (), s.end (), rank) == s.end ()) 
            {
              s.push_back( rank );
              (*hit).second->accessPllX ().setLinkage (s);
            }
          }
        }
      }
    }
  }

  void GitterPll::MacroGitterPll::vertexLinkageEstimate (MpAccessLocal & mpAccess) 
  {
    // for small processor numbers use gcollect( MPI_Allgather ) version 
    // this method should be faster (log p), 
    // but is more memory consuming O( p ) 
    if( ALUGridExternalParameters::useAllGather( mpAccess ) ) 
    {
      vertexLinkageEstimateGCollect ( mpAccess );
    }
    else 
    {
      // for larger processor numbers use bcast ( MPI_Bcast ) version 
      // this method is more time consuming (p log p)
      // but is the memory consumption is only O( 1 )
      vertexLinkageEstimateBcast ( mpAccess );
    }
  }

  void GitterPll::MacroGitterPll::identification (MpAccessLocal & c) 
  {
    // clear all entries and also clear memory be reassigning 
    vertexTT_t().swap( _vertexTT );
    hedgeTT_t ().swap( _hedgeTT );
    hfaceTT_t ().swap( _hfaceTT );

    // make sure the memory was deallocated 
    assert( _vertexTT.capacity() == 0 );
    assert( _hedgeTT.capacity()  == 0 );
    assert( _hfaceTT.capacity()  == 0 );

    c.removeLinkage ();
    
    int lap1 = clock ();
    vertexLinkageEstimate (c);

    c.insertRequestSymetric (secondScan ());
    if (debugOption (2)) c.printLinkage (std::cout);

    int lap2 = clock ();
    identify < vertex_STI > (AccessIterator < vertex_STI >::Handle (*this), _vertexTT, c);
    
    int lap3 = clock ();
    identify < hedge_STI > (AccessIterator < hedge_STI >::Handle (*this), _hedgeTT, c);
    
    int lap4 = clock ();
    identify < hface_STI > (AccessIterator < hface_STI >::Handle (*this), _hfaceTT, c);

    int lap5 = clock ();
    if (debugOption (2)) {
      float u2 = (float)(lap2 - lap1)/(float)(CLOCKS_PER_SEC);
      float u3 = (float)(lap3 - lap2)/(float)(CLOCKS_PER_SEC);
      float u4 = (float)(lap4 - lap3)/(float)(CLOCKS_PER_SEC);
      float u5 = (float)(lap5 - lap4)/(float)(CLOCKS_PER_SEC);
      float u6 = (float)(lap5 - lap1)/(float)(CLOCKS_PER_SEC);
      std::cout.precision (3);
      std::cout << "**INFO GitterPll::MacroGitterPll::identification () [lnk|vtx|edg|fce|all] ";
      std::cout << u2 << " " << u3 << " " << u4 << " " << u5 << " " << u6 << " sec." << std::endl;
    }
  }

} // namespace ALUGrid