#pragma once
#include "external.h"

// advection solver types
enum AdvectionSolverType
{
    TG2,
    CG2,
    TTG2,
    TTG3,
    TTG4
};

// adv sol type <-> adv sol name
inline std::map<AdvectionSolverType, std::string> AdvectionSolverTypeToName = 
{
    {AdvectionSolverType::TG2,  "TG2" },
    {AdvectionSolverType::CG2,  "CG2" },
    {AdvectionSolverType::TTG2, "TTG2"},
    {AdvectionSolverType::TTG3, "TTG3"},
    {AdvectionSolverType::TTG4, "TTG4"}
};

inline std::map<std::string, AdvectionSolverType> AdvectionSolverNameToType = 
{
    {"TG2",  AdvectionSolverType::TG2 },
    {"CG2",  AdvectionSolverType::CG2 },
    {"TTG2", AdvectionSolverType::TTG2},
    {"TTG3", AdvectionSolverType::TTG3},
    {"TTG4", AdvectionSolverType::TTG4}
};

// adv sol type -> is single step
inline std::map<AdvectionSolverType, bool> IsAdvectionSolverSingleStep = 
{
    {AdvectionSolverType::TG2,  true },
    {AdvectionSolverType::CG2,  true },
    {AdvectionSolverType::TTG2, false},
    {AdvectionSolverType::TTG3, false},
    {AdvectionSolverType::TTG4, false}
};

// momentum solver types
enum MomentumSolverType
{
    mEVP,
    aEVP,
    mEVPopt,
    JFNKS,
    Picard
};

// momentum sol type <-> momentum sol name
inline std::map<MomentumSolverType, std::string> MomentumSolverTypeToName = 
{
    {MomentumSolverType::mEVP,    "mEVP"   },
    {MomentumSolverType::aEVP,    "aEVP"   },
    {MomentumSolverType::mEVPopt, "mEVPopt"},
    {MomentumSolverType::JFNKS,   "JFNKS"  },
    {MomentumSolverType::Picard,  "Picard" }
};

inline std::map<std::string, MomentumSolverType> MomentumSolverNameToType = 
{
    {"mEVP",    MomentumSolverType::mEVP    },
    {"aEVP",    MomentumSolverType::aEVP    },
    {"mEVPopt", MomentumSolverType::mEVPopt },
    {"JFNKS",   MomentumSolverType::JFNKS   },
    {"Picard",  MomentumSolverType::Picard  }
};


// all about model variables 
enum ModelVariableNotation
{
    m,
    h,
    a,
    u_ice,
    sig1,
    sig2,
    sig12,
    eps1,
    eps2,
    eps12,
    del,
    P,
    P0,
    u_air,
    u_water,
    hw
};

// all model variables
inline constexpr std::initializer_list<ModelVariableNotation> 
ModelVariableNotationList = 
{
    m,
    h,
    a,
    u_ice,
    sig1,
    sig2,
    sig12,
    eps1,
    eps2,
    eps12,
    del,
    P,
    P0,
    u_air,
    u_water,
    hw
};

// node variable list
inline constexpr std::initializer_list<ModelVariableNotation> 
NodeModelVariableNotationList = 
{
    m,
    h,
    a,
    u_ice,
    u_air,
    u_water,
    hw
};

// triangle variable list
inline constexpr std::initializer_list<ModelVariableNotation> 
TriangleModelVariableNotationList = 
{
    sig1,
    sig2,
    sig12,
    eps1,
    eps2,
    eps12,
    del,
    P,
    P0
};

// vector variables list
inline constexpr std::initializer_list<ModelVariableNotation> 
NodeVectorNotationList = 
{
    u_ice,
    u_air,
    u_water,
};


// variable notation <-> variable name
inline std::map<ModelVariableNotation, std::string> ModelVariableNotationToName =
{
    {ModelVariableNotation::m,       "ice mass"          },
    {ModelVariableNotation::h,       "ice height"        },
    {ModelVariableNotation::a,       "ice concentration" },
    {ModelVariableNotation::u_ice,   "ice velocity"      },
    {ModelVariableNotation::sig1,    "ice stress 1"      },
    {ModelVariableNotation::sig2,    "ice stress 2"      },
    {ModelVariableNotation::sig12,   "ice stress 12"     },
    {ModelVariableNotation::eps1,    "ice strain rate 1" },
    {ModelVariableNotation::eps2,    "ice strain rate 2" },
    {ModelVariableNotation::eps12,   "ice strain rate 12"},
    {ModelVariableNotation::del,     "ice delta"         },
    {ModelVariableNotation::P,       "ice pressure"      },
    {ModelVariableNotation::P0,      "ice pressure 0"    },
    {ModelVariableNotation::u_air,   "air velocity"      },
    {ModelVariableNotation::u_water, "water velocity"    },
    {ModelVariableNotation::hw,      "water level"       }
};

inline std::map<std::string, ModelVariableNotation> ModelVariableNameToNotation =
{
    {"ice mass",             ModelVariableNotation::m      },
    {"ice height",            ModelVariableNotation::h     },
    {"ice concentration",    ModelVariableNotation::a      },
    {"ice velocity",         ModelVariableNotation::u_ice  },
    {"ice stress 1",         ModelVariableNotation::sig1   },
    {"ice stress 2",         ModelVariableNotation::sig2   },
    {"ice stress 12",        ModelVariableNotation::sig12  },
    {"ice strain rate 1",    ModelVariableNotation::eps1   },
    {"ice strain rate 2",    ModelVariableNotation::eps2   },
    {"ice strain rate 12",   ModelVariableNotation::eps12  },
    {"ice delta",            ModelVariableNotation::del    },
    {"ice pressure",         ModelVariableNotation::P      },
    {"ice pressure 0",       ModelVariableNotation::P0     },
    {"air velocity",         ModelVariableNotation::u_air  },
    {"water velocity",       ModelVariableNotation::u_water},
    {"water level",          ModelVariableNotation::hw     }
};

enum NodeCoordsNotation
{
    model,
    geo,
    Cartesian
};

inline std::map<std::string, NodeCoordsNotation> ModelCoordsNameToNotation =
{
    {"model coords",     NodeCoordsNotation::model     },
    {"geo coords",       NodeCoordsNotation::geo       },
    {"Cartesian coords", NodeCoordsNotation::Cartesian }
};

inline std::map<NodeCoordsNotation, std::string> ModelCoordsNotationToName =
{
    {NodeCoordsNotation::model,     "model coords"     },
    {NodeCoordsNotation::geo,       "geo coords"       },
    {NodeCoordsNotation::Cartesian, "Cartesian coords" }
};

enum SurfaceType
{
    plane,
    sphere
};

inline std::map<std::string, SurfaceType> SurfaceNameToType =
{
    {"plane",   SurfaceType::plane },
    {"sphere", SurfaceType::sphere },
};

inline std::map<SurfaceType, std::string> SurfaceTypeToName =
{
    { SurfaceType::plane,  "plane"  },
    { SurfaceType::sphere, "sphere" },
};

enum CoordsType
{
    Cartesian3D,
    Cartesian2D,
    SphericalRadians,
    SphericalDegrees,
    SphericalSpecial
};

inline std::map<std::string, CoordsType> CoordsNameToType =
{
    {"Cartesian3D",      CoordsType::Cartesian3D     },
    {"Cartesian2D",      CoordsType::Cartesian2D     },
    {"SphericalRadians", CoordsType::SphericalRadians},
    {"SphericalDegrees", CoordsType::SphericalDegrees},
    {"SphericalSpecial", CoordsType::SphericalSpecial}
};

inline std::map<CoordsType, std::string> CoordsTypeToName =
{
    {CoordsType::Cartesian3D     , "Cartesian3D"     },
    {CoordsType::Cartesian2D     , "Cartesian2D"     },
    {CoordsType::SphericalRadians, "SphericalRadians"},
    {CoordsType::SphericalDegrees, "SphericalDegrees"},
    {CoordsType::SphericalSpecial, "SphericalSpecial"}
};

enum MomentumBC
{
    no_slip
};

inline std::map<std::string, MomentumBC> MomentumNameToBC =
{
    {"no-slip", MomentumBC::no_slip}
};

inline std::map<MomentumBC, std::string> MomentumBCToName =
{
    {MomentumBC::no_slip, "no-slip"}
};