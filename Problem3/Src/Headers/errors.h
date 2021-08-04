#pragma once
#include "external.h"
#include "mesh.h"
#include "numerical_integration.h"
#include "function.h"

namespace INMOST_ICE
{
    void PrintErrors(IceMesh& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag);
}

double integral_over_subdomain_f(IceMesh& n, INMOST::Tag& m);
double integral_over_subdomain_abs_f1_minus_f2(IceMesh& n, INMOST::Tag& m1, INMOST::Tag& m2);
double integral_over_subdomain_squared_f1_minus_f2(IceMesh& n, INMOST::Tag& m1, INMOST::Tag& m2);

double subdomain_max_f(IceMesh& n, INMOST::Tag& m);
double subdomain_min_f(IceMesh& n, INMOST::Tag& m);
double subdomain_max_abs_f(IceMesh& n, INMOST::Tag& m);
double subdomain_max_abs_f1_minus_f2(IceMesh& n, INMOST::Tag& m1, INMOST::Tag& m2);
double max_f(IceMesh& n, INMOST::Tag& m);
double min_f(IceMesh& n, INMOST::Tag& m);

double l1_error(IceMesh& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag);
double l2_error(IceMesh& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag);
double linf_error(IceMesh& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag);
double phi_max_error(IceMesh& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag);
double phi_min_error(IceMesh& n, INMOST::Tag& m_tag, INMOST::Tag& m_init_tag);

double IntegralMass(IceMesh& n, INMOST::Tag& m_tag);
