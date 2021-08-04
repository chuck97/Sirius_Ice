#pragma once
#include "external.h"
#include "mesh.h"
#include "advection_assembling.h"

void FCT_Procedure(IceMesh& n,
                   const std::vector<std::vector<std::vector<double>>>& M_C_minus_M_L,
                   INMOST::Tag m_tag,
                   INMOST::Tag m_low_tag,
                   INMOST::Tag m_high_tag,
                   double fct_cd);

void CalculateAlpha(IceMesh& n,
                    INMOST::Tag& Fe_trian_tag,
                    INMOST::Tag& alpha_trian_tag,
                    INMOST::Tag m_tag,
                    INMOST::Tag m_low_tag);