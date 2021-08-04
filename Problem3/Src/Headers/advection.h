#pragma once

#include "external.h"
#include "enum_classes.h"
#include "config.h"
#include "mesh.h"
#include "errors.h"
#include "advection_assembling.h"
#include "fct.h"

enum class StepNumber
{
    first,
    second
};

class AdvectionSolver
{
public:
    AdvectionSolver(IceMesh& im,
                    AdvectionParams& ap,
                    ModelParams& mp,
                    MeshParams& mep,
                    ModelVariableNotation transported_scalar_,
                    ModelVariableNotation transporting_velocity_,
                    INMOST::Solver& solver,
                    bool is_verbose_advection_);
    void AssembleLHS();
    void AssembleRHS();
    void Evaluate();

private:
    void AssembleSingleStepRHS();
    void AssembleDoubleStepRHS(StepNumber step_num);
    void LogInit();
    void LogStep();

    IceMesh& ice_mesh;
    AdvectionParams& advection_params;
    ModelParams& model_params;
    MeshParams& mesh_params;
    ModelVariableNotation transported_scalar;
    ModelVariableNotation transporting_velocity;
    INMOST::Tag m_tag;
    INMOST::Tag u_tag;
    INMOST::Tag m_high_tag;
    INMOST::Tag m_low_tag;
    INMOST::Tag m_half_tag;
    INMOST::Sparse::Matrix LHS;
    INMOST::Sparse::Matrix LHS_low;
    INMOST::Sparse::Vector RHS;
    INMOST::Sparse::Vector RHS_low;
    INMOST::Solver& Sol;
    std::vector<std::vector<std::vector<double>>> M_C_minus_M_L;
    bool is_verbose_advection = true;
    double init_mass_integral;
};