from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

from dune.grid.grid_generator import register

register( ALU        = "dune.alugrid.alu",
          ALUConform = "dune.alugrid.aluconform",
          ALUCube    = "dune.alugrid.alucube",
          ALUSimplex = "dune.alugrid.alusimplex" )
