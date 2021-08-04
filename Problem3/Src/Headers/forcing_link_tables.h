#pragma once
#include "external.h"
#include "service.h"
#include "enum_classes.h"


// topaz -> model
inline std::map<std::string, MomentumSolverType> TopazVariableToModel =
{
    {"mass",             ModelVariableNotation::m    },
    {"hight",            ModelVariableNotation::h    },
    {"concentration",    ModelVariableNotation::a    },
    {"velocity x",       ModelVariableNotation::u    },
    {"velocity y",       ModelVariableNotation::v    },
    {"stress 1",         ModelVariableNotation::sig1 },
    {"stress 2",         ModelVariableNotation::sig2 },
    {"stress 12",        ModelVariableNotation::sig12},
    {"strain rate 1",    ModelVariableNotation::eps1 },
    {"strain rate 2",    ModelVariableNotation::eps2 },
    {"strain rate 12",   ModelVariableNotation::eps12},
    {"delta",            ModelVariableNotation::del  },
    {"pressure",         ModelVariableNotation::P    },
    {"pressure 0",       ModelVariableNotation::P0   },
    {"air velocity x",   ModelVariableNotation::ua   },
    {"air velocity y",   ModelVariableNotation::va   },
    {"water velocity x", ModelVariableNotation::uw   },
    {"water velocity y", ModelVariableNotation::vw   },
    {"water level",      ModelVariableNotation::hw   }
}