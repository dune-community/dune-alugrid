// (c) mario ohlberger 1998
#include <config.h>

#include <dune/alugrid/common/alugrid_assert.hh>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include "mapp_tetra_3d.h"

namespace ALUGrid
{

  const alucoord_t quadraturTetra3Dbasis :: _p1 [4] =  {0.25, 0.25, 0.25, 0.25} ;

  const alucoord_t quadraturTetra3Dbasis :: _w2 [4] = {.0416666666,
                                                   .0416666666,
                                                   .0416666666,
                                                   .0416666666} ;

  const alucoord_t quadraturTetra3Dbasis :: _p2 [4][4] = {{.1381966, .1381966, .1381966, .5854102},
                                                      {.1381966, .1381966, .5854102, .1381966},
                                                      {.1381966, .5854102, .1381966, .1381966},
                                                      {.5854102, .1381966, .1381966, .1381966}};
    
  const alucoord_t quadraturTetra3Dbasis :: _w7 [64] = {0.0026134590, 
                                                    0.0048996145,
                                                    0.0048996145,
                                                    0.0026134590,
                                                    0.0039241268,
                                                    0.0073568050,
                                                    0.0073568050,
                                                    0.0039241268, 
                                                    0.0025043094, 
                                                    0.0046949850,
                                                    0.0046949850,
                                                    0.0025043094, 
                                                    0.0006013729,
                                                    0.0011274313,
                                                    0.0011274313,
                                                    0.0006013729,
                                                    0.0033810896,
                                                    0.0063387393,
                                                    0.0063387393,
                                                    0.0033810896,
                                                    0.0050767294,
                                                    0.0095176610,
                                                    0.0095176610,
                                                    0.0050767294,
                                                    0.0032398804,
                                                    0.0060740056,
                                                    0.0060740056,
                                                    0.0032398804,
                                                    0.0007780094, 
                                                    0.0014585828,
                                                    0.0014585828, 
                                                    0.0007780094,
                                                    0.0016175887,
                                                    0.0030325944,
                                                    0.0030325944,
                                                    0.0016175887,
                                                    0.0024288207,
                                                    0.0045534614,
                                                    0.0045534614,
                                                    0.0024288207,
                                                    0.0015500311,
                                                    0.0029059399,
                                                    0.0029059399,
                                                    0.0015500311,
                                                    0.0003722171,
                                                    0.0006978185,
                                                    0.0006978185,
                                                    0.0003722171,
                                                    0.0002439854,
                                                    0.0004574147,
                                                    0.0004574147,
                                                    0.0002439854,
                                                    0.0003663458,
                                                    0.0006868113,
                                                    0.0006868113,
                                                    0.0003663458,
                                                    0.0002337955,
                                                    0.0004383110,
                                                    0.0004383110,
                                                    0.0002337955,
                                                    0.0000561425,
                                                    0.0001052539,
                                                    0.0001052539,
                                                    0.0000561425} ;

  const alucoord_t quadraturTetra3Dbasis :: _p7 [64][4] = {{0.0485005494,0.0543346112,0.0622918076,0.8348730318},
                                                      {0.0485005494,0.0543346112,0.2960729005,0.6010919389},
                                                      {0.0485005494,0.0543346112,0.6010919389,0.2960729005},
                                                      {0.0485005494,0.0543346112,0.8348730300,0.0622918093},
                                                      {0.0485005494,0.2634159753,0.0477749033,0.6403085720},
                                                      {0.0485005494,0.2634159753,0.2270740686,0.4610094066},
                                                      {0.0485005494,0.2634159753,0.4610094066,0.2270740686},
                                                      {0.0485005494,0.2634159753,0.6403085706,0.0477749047},
                                                      {0.0485005494,0.5552859758,0.0275098315,0.3687036433},
                                                      {0.0485005494,0.5552859758,0.1307542021,0.2654592727},
                                                      {0.0485005494,0.5552859758,0.2654592727,0.1307542021},
                                                      {0.0485005494,0.5552859758,0.3687036425,0.0275098323},
                                                      {0.0485005494,0.8185180165,0.0092331459,0.1237482881},
                                                      {0.0485005494,0.8185180165,0.0438851337,0.0890963004},
                                                      {0.0485005494,0.8185180165,0.0890963004,0.0438851337},
                                                      {0.0485005494,0.8185180165,0.1237482879,0.0092331462},
                                                      {0.2386007376,0.0434790928,0.0498465199,0.6680736497},
                                                      {0.2386007376,0.0434790928,0.2369204606,0.4809997090},
                                                      {0.2386007376,0.0434790928,0.4809997090,0.2369204606},
                                                      {0.2386007376,0.0434790928,0.6680736482,0.0498465214},
                                                      {0.2386007376,0.2107880664,0.0382299497,0.5123812464},
                                                      {0.2386007376,0.2107880664,0.1817069135,0.3689042825},
                                                      {0.2386007376,0.2107880664,0.3689042825,0.1817069135},
                                                      {0.2386007376,0.2107880664,0.5123812453,0.0382299508},
                                                      {0.2386007376,0.4443453248,0.0220136390,0.2950402987},
                                                      {0.2386007376,0.4443453248,0.1046308045,0.2124231331},
                                                      {0.2386007376,0.4443453248,0.2124231331,0.1046308045},
                                                      {0.2386007376,0.4443453248,0.2950402980,0.0220136396},
                                                      {0.2386007376,0.6549862048,0.0073884546,0.0990246030},
                                                      {0.2386007376,0.6549862048,0.0351173176,0.0712957400},
                                                      {0.2386007376,0.6549862048,0.0712957400,0.0351173176},
                                                      {0.2386007376,0.6549862048,0.0990246028,0.0073884548},
                                                      {0.5170472951,0.0275786260,0.0316174612,0.4237566177},
                                                      {0.5170472951,0.0275786260,0.1502777622,0.3050963168},
                                                      {0.5170472951,0.0275786260,0.3050963168,0.1502777622},
                                                      {0.5170472951,0.0275786260,0.4237566168,0.0316174621},
                                                      {0.5170472951,0.1337020823,0.0242491141,0.3250015085},
                                                      {0.5170472951,0.1337020823,0.1152560157,0.2339946069},
                                                      {0.5170472951,0.1337020823,0.2339946069,0.1152560157},
                                                      {0.5170472951,0.1337020823,0.3250015078,0.0242491148},
                                                      {0.5170472951,0.2818465779,0.0139631689,0.1871429581},
                                                      {0.5170472951,0.2818465779,0.0663669280,0.1347391990},
                                                      {0.5170472951,0.2818465779,0.1347391990,0.0663669280},
                                                      {0.5170472951,0.2818465779,0.1871429577,0.0139631693},
                                                      {0.5170472951,0.4154553004,0.0046864691,0.0628109354},
                                                      {0.5170472951,0.4154553004,0.0222747832,0.0452226213},
                                                      {0.5170472951,0.4154553004,0.0452226213,0.0222747832},
                                                      {0.5170472951,0.4154553004,0.0628109352,0.0046864693},
                                                      {0.7958514179,0.0116577407,0.0133649937,0.1791258477},
                                                      {0.7958514179,0.0116577407,0.0635238021,0.1289670393},
                                                      {0.7958514179,0.0116577407,0.1289670393,0.0635238021},
                                                      {0.7958514179,0.0116577407,0.1791258473,0.0133649941},
                                                      {0.7958514179,0.0565171087,0.0102503252,0.1373811482},
                                                      {0.7958514179,0.0565171087,0.0487197855,0.0989116879},
                                                      {0.7958514179,0.0565171087,0.0989116879,0.0487197855},
                                                      {0.7958514179,0.0565171087,0.1373811479,0.0102503255},
                                                      {0.7958514179,0.1191391593,0.0059023608,0.0791070620},
                                                      {0.7958514179,0.1191391593,0.0280539153,0.0569555075},
                                                      {0.7958514179,0.1191391593,0.0569555075,0.0280539153},
                                                      {0.7958514179,0.1191391593,0.0791070618,0.0059023610},
                                                      {0.7958514179,0.1756168040,0.0019810139,0.0265507642},
                                                      {0.7958514179,0.1756168040,0.0094157572,0.0191160209},
                                                      {0.7958514179,0.1756168040,0.0191160209,0.0094157572},
                                                      {0.7958514179,0.1756168040,0.0265507642,0.0019810140}};

  const alucoord_t quadraturTriang2Dbasis :: _p1 [3] = { 1.0/3.0, 1.0/3.0, 1.0/3.0 } ;

  const alucoord_t quadraturTriang2Dbasis :: _w3 [7] = {0.025,
                                                    0.025,
                                                    0.025,
                                                    0.066666666,
                                                    0.066666666,
                                                    0.066666666,
                                                    0.225} ;

  const alucoord_t quadraturTriang2Dbasis :: _p3 [7][3] = {{ .0, .0, 1.},
                                                       { .0, 1., .0},
                                                       { 1., .0, .0},
                                                       {.0, .5, .5 },
                                                       { .5, .0, .5},
                                                       { .5, .5, .0},
                                                       {.333333333, .333333333, .333333333}} ;

  const alucoord_t quadraturTriang2Dbasis :: _w5 [7] = {.11250000,
                                                    .06296959,
                                                    .06296959,
                                                    .06296959,
                                                    .06619708,
                                                    .06619708,
                                                    .06619708} ;

  const alucoord_t quadraturTriang2Dbasis :: _p5 [7][3] = {{.333333333, .333333333, .333333333},
                                                       {.101286507, .101286507, .797426985},
                                                       {.101286507, .797426985, .101286507},
                                                       {.797426985, .101286507, .101286507},
                                                       {.470142064, .470142064, .059715871},
                                                       {.470142064, .059715871, .470142064},
                                                       {.059715871, .470142064, .470142064}} ;

  const alucoord_t quadraturTriang2Dbasis :: _w7 [16] = {.023568368192,
                                                     .044185088508,
                                                     .044185088508,
                                                     .023568368192,
                                                     .035388067903,
                                                     .066344216097,
                                                     .066344216097,
                                                     .035388067903,
                                                     .022584049285,
                                                     .042339724515,
                                                     .042339724515,
                                                     .022584049285,
                                                     .005423225903,
                                                     .010167259547,
                                                     .010167259547,
                                                     .005423225903} ;
                                                       
  const alucoord_t quadraturTriang2Dbasis :: _p7 [16][3] = {{.057104196, .065466992, .877428812},
                                                        {.057104196, .311164552, .631731251},
                                                        {.057104196, .631731250, .311164553},
                                                        {.057104196, .877428808, .065466995},
                                                        {.276843013, .050210121, .672946865},
                                                        {.276843013, .238648659, .484508327},
                                                        {.276843013, .484508326, .238648660},
                                                        {.276843013, .672946863, .050210123},
                                                        {.583590432, .028912083, .387497484},
                                                        {.583590432, .137419104, .278990463},
                                                        {.583590432, .278990463, .137419105},
                                                        {.583590432, .387497483, .028912084},
                                                        {.860240136, .009703785, .130056079},
                                                        {.860240136, .046122079, .093637784},
                                                        {.860240136, .093637784, .046122079},
                                                        {.860240136, .130056078, .009703785}} ;


  void LinearMapping :: inverse() {
          //  Kramer - Regel
    alucoord_t val = 1.0 / det () ;
    Dfi[0][0] = ( Df[1][1] * Df[2][2] - Df[1][2] * Df[2][1] ) * val ;
    Dfi[0][1] = ( Df[0][2] * Df[2][1] - Df[0][1] * Df[2][2] ) * val ;
    Dfi[0][2] = ( Df[0][1] * Df[1][2] - Df[0][2] * Df[1][1] ) * val ;
    Dfi[1][0] = ( Df[1][2] * Df[2][0] - Df[1][0] * Df[2][2] ) * val ;
    Dfi[1][1] = ( Df[0][0] * Df[2][2] - Df[0][2] * Df[2][0] ) * val ;
    Dfi[1][2] = ( Df[0][2] * Df[1][0] - Df[0][0] * Df[1][2] ) * val ;
    Dfi[2][0] = ( Df[1][0] * Df[2][1] - Df[1][1] * Df[2][0] ) * val ;
    Dfi[2][1] = ( Df[0][1] * Df[2][0] - Df[0][0] * Df[2][1] ) * val ;
    Dfi[2][2] = ( Df[0][0] * Df[1][1] - Df[0][1] * Df[1][0] ) * val ;
  }

  void LinearMapping :: world2map (const alucoord_t (&wld)[3], alucoord_t (&map)[4]) {
    map [0] = map [1] = map [2] = map [3] = .0 ;
    alucoord_t upd [3] ;
    map2world (map, upd) ;
    inverse () ;
    alucoord_t u0 = wld [0] - upd [0] ;
    alucoord_t u1 = wld [1] - upd [1] ;
    alucoord_t u2 = wld [2] - upd [2] ;
    alucoord_t c0 = Dfi [0][0] * u0 + Dfi [0][1] * u1 + Dfi [0][2] * u2 ;
    alucoord_t c1 = Dfi [1][0] * u0 + Dfi [1][1] * u1 + Dfi [1][2] * u2 ;
    alucoord_t c2 = Dfi [2][0] * u0 + Dfi [2][1] * u1 + Dfi [2][2] * u2 ;
    map [0] = c0 ;
    map [1] = c1 ;
    map [2] = c2 ;
  }

} // namespace ALUGrid
