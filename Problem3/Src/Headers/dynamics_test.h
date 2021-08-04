#pragma once
#include "external.h"
#include "enum_classes.h"
#include "config.h"
#include "mesh.h"
#include "config.h"

class DynamicsTest
{
public:
    DynamicsTest(IceMesh& im,
                 DynamicsTestParams& dtp,
                 ModelParams& mp);
    void AssignInitialScalars();
    void AssignInitialVectors();
    void UpdateWaterVelocity();
    void UpdateAirVelocity(double current_time_hours);
private:
    IceMesh& ice_mesh;
    DynamicsTestParams& dynamics_test_params;
    ModelParams& model_params;
};