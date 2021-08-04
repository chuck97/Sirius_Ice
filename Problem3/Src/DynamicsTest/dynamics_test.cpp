#include "dynamics_test.h"

using namespace std;
using namespace INMOST;

DynamicsTest::DynamicsTest(IceMesh& im,
                           DynamicsTestParams& dtp,
                           ModelParams& mp):
                           ice_mesh(im),
                           dynamics_test_params(dtp),
                           model_params(mp)
{
};

void DynamicsTest::AssignInitialScalars()
{
    // initialize a,h
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
            nodeit != ice_mesh.GetMesh()->EndNode();
            ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            // assign a
            nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a]) =
            dynamics_test_params.GetInitialIceConcentration();

            // assign h
            double x = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::Cartesian])[0];
            double y = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::Cartesian])[1];
            double h_background = dynamics_test_params.GetInitialHeightBackground();
            double h_scale = dynamics_test_params.GetInitialHeightScaleFactor();
            double coeff_x = dynamics_test_params.GetInitialHeightXfactor();
            double coeff_y = dynamics_test_params.GetInitialHeightYfactor();

            nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::h]) = 
            h_background + h_scale*(sin(coeff_x*x) + sin(coeff_y*y));

            // calculate m
            nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::m]) =
            nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::a])*
            nodeit->Real(ice_mesh.GetData().NodeData[ModelVariableNotation::h])*
            model_params.GetIceDensity();
        }
    };
    // exchange a, h, m
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::a], NODE, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::h], NODE, 0);
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::m], NODE, 0);
    BARRIER
};

void DynamicsTest::AssignInitialVectors()
{
    // initialize u ice
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
            nodeit != ice_mesh.GetMesh()->EndNode();
            ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_ice])[0] = 0.0;
            nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_ice])[1] = 0.0;
        }
    };
    // exchange u ice
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::u_ice], NODE, 0);
    BARRIER
};

void DynamicsTest::UpdateWaterVelocity()
{
    // update u water
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
            nodeit != ice_mesh.GetMesh()->EndNode();
            ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            double x = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::Cartesian])[0];
            double y = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::Cartesian])[1];
            double v_water_max = dynamics_test_params.GetMaxWaterSpeed();
            double L = dynamics_test_params.GetDomainSize();

            nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_water])[0] = 
            v_water_max*(2.0*y - L)/L;
            nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_water])[1] =
            v_water_max*(L - 2.0*x)/L;
        }
    };
    // exchange u water
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::u_water], NODE, 0);
    BARRIER
}

void DynamicsTest::UpdateAirVelocity(double current_time_hours)
{
    double t = current_time_hours;
    double m_x = 256000.0 + 51200.0*(t/24.0);
    double m_y = m_x;
    double alpha = dynamics_test_params.GetAirConvergenceAngle();
    double v_air_max = dynamics_test_params.GetMaxAirSpeed();
    double air_reduction = dynamics_test_params.GetAirReductionFactor();

    // update u air
    for(Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
            nodeit != ice_mesh.GetMesh()->EndNode();
            ++nodeit)
    {
        if (nodeit->GetStatus() != Element::Ghost)
        {
            double x = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::Cartesian])[0];
            double y = nodeit->RealArray(ice_mesh.GetCoords()[NodeCoordsNotation::Cartesian])[1];
            double r = sqrt((m_x-x)*(m_x-x) + (m_y-y)*(m_y-y));
            double s = air_reduction*exp(-0.00001*r);

            nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_air])[0] = 
            -s*v_air_max*(cos(alpha)*(x-m_x)/1000.0 + sin(alpha)*(y-m_y)/1000.0);
            nodeit->RealArray(ice_mesh.GetData().NodeData[ModelVariableNotation::u_air])[1] =
            -s*v_air_max*(-sin(alpha)*(x-m_x)/1000.0 + cos(alpha)*(y-m_y)/1000.0);
        }
    };
    // exchange u air
    BARRIER
    ice_mesh.GetMesh()->ExchangeData(ice_mesh.GetData().NodeData[ModelVariableNotation::u_air], NODE, 0);
    BARRIER
}