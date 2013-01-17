#ifndef _ALUGRID_PARALLEL_h_INCLUDED_
#define _ALUGRID_PARALLEL_h_INCLUDED_

#include "alugrid_serial.h"

#define ALUGRID_EXPORT_MACROGRID_CHANGES

namespace ALUGridSpace {

// the parallel stuff 
#include "parallel/gitter_pll_sti.h"
#include "parallel/gitter_pll_impl.h"
#include "parallel/gitter_pll_ldb.h"
#include "parallel/gitter_tetra_top_pll.h"
#include "parallel/gitter_hexa_top_pll.h"
#include "parallel/mpAccess.h"
#include "parallel/mpAccess_MPI.h"
#include "parallel/mpAccess_STAR.h"
#include "parallel/gitter_pll_mgb.h"

// the duneinterface stuff 
#include "duneinterface/gitter_dune_pll_impl.h"

}
#endif