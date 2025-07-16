//
// Example based on the GEOS-Chem carbon gases mechanism
//

#include atoms.kpp      { Periodic table information              }

#DEFVAR
CH4          = IGNORE;  { Active methane species                  }
CO           = IGNORE;  { Active carbon monoxide species          }
PCOfromCH4   = IGNORE;  { Tracks P(CO) from CH4   for diagnostics }
PCOfromNMVOC = IGNORE;  { Tracks P(CO) from NMVOC for diagnostics }
LCH4byOH     = IGNORE;  { Dummy spc to track loss of CH4 by OH    }
LCH4byCl     = IGNORE;  { Dummy spc to track loss of CH4 by Cl    }
LCObyOH      = IGNORE;  { Dummy spc to track loss of CO  by OH    }

#DEFFIX
FixedOH      = IGNORE;  { Externally-supplied OH concentration    }
FixedCl      = IGNORE;  { Externally-supplied Cl concentration    }
DummyCH4     = IGNORE;  { Dummy placeholder for CH4   reactant    }
DummyNMVOC   = IGNORE;  { Dummy placeholder for NMVOC reactant    }
