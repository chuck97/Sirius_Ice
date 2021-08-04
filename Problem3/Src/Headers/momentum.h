#pragma once

#include "external.h"
#include "enum_classes.h"
#include "config.h"
#include "mesh.h"
#include "advection_assembling.h"
#include "function.h"
#include "momentum_assembling.h"

class MomentumSolver 
{
public:
    MomentumSolver(IceMesh& im,
                   MomentumParams& mp,
                   ModelParams& modp,
                   MeshParams& mep,
                   INMOST::Solver& solver,
                   bool is_verbose_momentum_);

    void Evaluate();

private:
    
    void LogInit();
    void AssembleLHS();
    void mEVP_evaluation();

    void AssignSigma();
    void AssignVelocity();
    void CalculateP_0();
    void CalculateVarepsilonDelta();
    void AssembleLevelVector();

    void UpdateTemporalSigma();
    void UpdateTemporalVelocity();
    void UpdateSigma();
    void AssembleForceVector();
    void UpdateVelocity();
    void UpdateScalars();
    void UpdateShearDeformation();
    void LogMevpError(int stepn);

    IceMesh& ice_mesh;
    MomentumParams& momentum_params;
    ModelParams& model_params;
    MeshParams& mesh_params;
    INMOST::Tag u_old_tag;
    INMOST::Tag u_n_tag;
    INMOST::Tag u_new_tag;
    INMOST::Tag sigma_tag;
    INMOST::Tag delta_tag;
    INMOST::Tag P_0_tag;
    INMOST::Tag varepsilon_tag;
    INMOST::Tag sigma_diff_tag;
    INMOST::Tag u_diff_tag;
    INMOST::Tag lumped_mass_matrix_entry_tag;
    INMOST::Tag Nx_matrix_entries_tag;
    INMOST::Tag Ny_matrix_entries_tag;
    INMOST::Tag Force_vector_tag;
    INMOST::Tag Level_vector_tag;
    INMOST::Tag shear_deformation_tag;
    INMOST::Solver& Sol;
    bool is_verbose_momentum = true;
};