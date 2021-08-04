#include "advection.h"

using namespace std;
using namespace INMOST;

AdvectionSolver::AdvectionSolver(IceMesh& im,
                                 AdvectionParams& ap,
                                 ModelParams& mp,
                                 MeshParams& mep,
                                 ModelVariableNotation transported_scalar_,
                                 ModelVariableNotation transporting_velocity_,
                                 INMOST::Solver& solver,
                                 bool is_verbose_advection_)
    : ice_mesh(im),
      advection_params(ap),
      model_params(mp),
      mesh_params(mep),
      Sol(solver),
      is_verbose_advection(is_verbose_advection_),
      transported_scalar(transported_scalar_),
      transporting_velocity(transporting_velocity_)
{ 
    m_tag = ice_mesh.GetData().NodeData[transported_scalar];
    u_tag = ice_mesh.GetData().NodeData[transporting_velocity];
    BARRIER
    LHS.SetInterval(ice_mesh.GetIdIntervalGlobal().IdMin, ice_mesh.GetIdIntervalGlobal().IdMax);
    LHS_low.SetInterval(ice_mesh.GetIdIntervalGlobal().IdMin, ice_mesh.GetIdIntervalGlobal().IdMax);
    RHS.SetInterval(ice_mesh.GetIdIntervalGlobal().IdMin, ice_mesh.GetIdIntervalGlobal().IdMax);
    RHS_low.SetInterval(ice_mesh.GetIdIntervalGlobal().IdMin, ice_mesh.GetIdIntervalGlobal().IdMax);
    init_mass_integral = IntegralMass(ice_mesh, m_tag);
    LogInit();
}

void AdvectionSolver::AssembleLHS()
{
    INMOST::Tag M_L_entire = ice_mesh.GetMesh()->CreateTag("M_L_entire", DATA_REAL, NODE, NONE, 1);

    for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
            trianit != ice_mesh.GetMesh()->EndCell();
            ++trianit) 
    {
        // get node coords of trian
        INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
        double v0_x = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[0];
        double v0_y = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[1];
        double v1_x = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[2];
        double v1_y = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[3];
        double v2_x = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[4];
        double v2_y = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[5];
        std::vector<std::vector<double>> local_node_coords = {{v0_x, v0_y},
                                                              {v1_x, v1_y},
                                                              {v2_x, v2_y}};
        
        // find global ids of local nodes
        std::vector<int> nums_global = { local_nodes[0].Integer(ice_mesh.GetData().NodeIdGlobal),
                                         local_nodes[1].Integer(ice_mesh.GetData().NodeIdGlobal),
                                         local_nodes[2].Integer(ice_mesh.GetData().NodeIdGlobal) }; 

        // assemble local 3x3 stiffness matrix
        std::vector<std::vector<double>> localLHS = LocalStiffnessMatrixAssembling(local_node_coords);
        if (advection_params.GetIsFct())
        {
            if (trianit->GetStatus()!= Element::Ghost)
            {
                M_C_minus_M_L.push_back(localLHS);
            } 
        }

        // assemble global LHS mass matrix
        for (int i = 0; i < 3; ++i)
        {
            if(local_nodes[i].GetStatus() != Element::Ghost)
            {
                for (int j = 0; j < 3; ++j)
                {
                    LHS[nums_global[i]][nums_global[j]] += localLHS[i][j];
                }
            }
        }
        if (advection_params.GetIsFct())
        {
            for (unsigned int i = 0; i < 3; ++i)
            {
                if(local_nodes[i].GetStatus() != Element::Ghost)
                {
                    for (unsigned int j = 0; j < 3; ++j)
                    {
                        LHS_low[nums_global[i]][nums_global[i]] += localLHS[i][j];
                        local_nodes[i]->Real(M_L_entire) += localLHS[i][j];
                    }
                }
            }
            if (trianit->GetStatus() != Element::Ghost)
            {
                for (unsigned int i = 0; i < 3; ++i)
                {
                    for (unsigned int j = 0; j < 3; ++j)
                    {
                        M_C_minus_M_L.back()[i][i] -= localLHS[i][j];
                    }
                }
            }
        } 
    }   
    BARRIER
    if (advection_params.GetIsFct())
    {
        ice_mesh.GetMesh()->ExchangeData(M_L_entire, NODE, 0);
        BARRIER
    }

    // multiply MC_minus_ML by ML^-1 for FCT
    if (advection_params.GetIsFct())
    {
        int t_num = 0;
        for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
                trianit != ice_mesh.GetMesh()->EndCell();
                ++trianit) 
        {
            if (trianit->GetStatus() != Element::Ghost)
            {
                INMOST::ElementArray<Node> local_nodes = trianit->getNodes(); 
                for (unsigned int i = 0; i < 3; ++i)
                {
                    for (unsigned int j = 0; j < 3; ++j)
                    {
                        M_C_minus_M_L[t_num][i][j] /= -local_nodes[i]->Real(M_L_entire);
                    }
                }
                ++t_num;
            }
        }
    }
    BARRIER
    ice_mesh.GetMesh()->DeleteTag(M_L_entire, NODE);
}

void AdvectionSolver::AssembleSingleStepRHS()
{
    INMOST::Tag id = ice_mesh.GetData().NodeIdGlobal;

    for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
        trianit != ice_mesh.GetMesh()->EndCell();
        ++trianit) 
    {
        // get node coords of trian
        INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
        double v0_x = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[0];
        double v0_y = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[1];
        double v1_x = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[2];
        double v1_y = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[3];
        double v2_x = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[4];
        double v2_y = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[5];
        std::vector<std::vector<double>> local_node_coords = {{v0_x, v0_y},
                                                              {v1_x, v1_y},
                                                              {v2_x, v2_y}};
        // get u values of trian
        std::vector<std::vector<double>> local_u_geo;
        if (mesh_params.GetCoordsType() == CoordsType::Cartesian2D)
        {
            local_u_geo = {{local_nodes[0].RealArray(u_tag)[0], local_nodes[0].RealArray(u_tag)[1]},
                           {local_nodes[1].RealArray(u_tag)[0], local_nodes[1].RealArray(u_tag)[1]},
                           {local_nodes[2].RealArray(u_tag)[0], local_nodes[2].RealArray(u_tag)[1]}};
        }
        else
        {
            INMOST_ICE_ERR("only Cartesian2D coords are avalible for Advection Solver right now");
        }
  
        std::vector<double> u0_local = 
        {
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[0]*local_u_geo[0][0] +
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[1]*local_u_geo[0][1],
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[2]*local_u_geo[0][0] +
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[3]*local_u_geo[0][1]
        };

        std::vector<double> u1_local = 
        {
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[4]*local_u_geo[1][0] +
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[5]*local_u_geo[1][1],
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[6]*local_u_geo[1][0] +
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[7]*local_u_geo[1][1]
        };
        std::vector<double> u2_local = 
        {
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[8]*local_u_geo[2][0] +
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[9]*local_u_geo[2][1],
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[10]*local_u_geo[2][0] +
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[11]*local_u_geo[2][1]
        };
        std::vector<std::vector<double>> local_u = {u0_local, u1_local, u2_local};
        
        // get m values of trian
        std::vector<double> local_m = 
        {local_nodes[0].Real(m_tag),
         local_nodes[1].Real(m_tag),
         local_nodes[2].Real(m_tag)};
        std::vector<double> localRHS;
        
        // assemble local load vector of size 3
        if (advection_params.GetAdvectionSolverType() == AdvectionSolverType::TG2)
        {
            localRHS = LocalTG2RhsAssembling(local_node_coords,
                                             local_u,
                                             local_m,
                                             model_params.GetTimeStepHours()*3600.0);
        }
        else if (advection_params.GetAdvectionSolverType() == AdvectionSolverType::CG2)
        {
            localRHS = LocalCG2RhsAssembling(local_node_coords,
                                             local_u,
                                             local_m,
                                             model_params.GetTimeStepHours()*3600.0);
        }
        else
        {
            INMOST_ICE_ERR("Only avalible single step solvers: TG2, CG2");
        }
        // assemble global load vector
        for (unsigned int i = 0; i < 3; ++i)
        {
            if(local_nodes[i].GetStatus() != Element::Ghost)
            {
                RHS[local_nodes[i].Integer(id)] += localRHS[i];
            }
        }
    }   
}   

void AdvectionSolver::AssembleDoubleStepRHS(StepNumber step_num)
{
    INMOST::Tag id = ice_mesh.GetData().NodeIdGlobal;
    
    for(Mesh::iteratorCell trianit = ice_mesh.GetMesh()->BeginCell();
        trianit != ice_mesh.GetMesh()->EndCell();
        ++trianit) 
    {
        // get node coords of trian
        INMOST::ElementArray<Node> local_nodes = trianit->getNodes();
        double v0_x = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[0];
        double v0_y = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[1];
        double v1_x = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[2];
        double v1_y = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[3];
        double v2_x = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[4];
        double v2_y = trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.LocalNodeCoords)[5];
        std::vector<std::vector<double>> local_node_coords = {{v0_x, v0_y},
                                                              {v1_x, v1_y},
                                                              {v2_x, v2_y}};


        // get u values of trian
        std::vector<std::vector<double>> local_u_geo;
        if (mesh_params.GetCoordsType() == CoordsType::Cartesian2D)
        {
            local_u_geo = {{local_nodes[0].RealArray(u_tag)[0], local_nodes[0].RealArray(u_tag)[1]},
                           {local_nodes[1].RealArray(u_tag)[0], local_nodes[1].RealArray(u_tag)[1]},
                           {local_nodes[2].RealArray(u_tag)[0], local_nodes[2].RealArray(u_tag)[1]}};
        }
        else
        {
            INMOST_ICE_ERR("only Cartesian2D coords are avalible for Advection Solver right now");
        }
  
        std::vector<double> u0_local = 
        {
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[0]*local_u_geo[0][0] +
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[1]*local_u_geo[0][1],
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[2]*local_u_geo[0][0] +
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[3]*local_u_geo[0][1]
        };

        std::vector<double> u1_local = 
        {
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[4]*local_u_geo[1][0] +
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[5]*local_u_geo[1][1],
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[6]*local_u_geo[1][0] +
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[7]*local_u_geo[1][1]
        };
        std::vector<double> u2_local = 
        {
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[8]*local_u_geo[2][0] +
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[9]*local_u_geo[2][1],
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[10]*local_u_geo[2][0] +
            trianit->RealArray(ice_mesh.GetLocalBasisData().triangle_data.VecTransMatriciesFromNodalToTriangle)[11]*local_u_geo[2][1]
        };
        std::vector<std::vector<double>> local_u = {u0_local, u1_local, u2_local};

        // get m values of trian
        std::vector<double> local_m = {local_nodes[0].Real(m_tag),
                                       local_nodes[1].Real(m_tag),
                                       local_nodes[2].Real(m_tag)};

        std::vector<double> local_m_half;
        if (step_num == StepNumber::second)
        {
            local_m_half = 
            {local_nodes[0].Real(m_half_tag),
             local_nodes[1].Real(m_half_tag),
             local_nodes[2].Real(m_half_tag)};
        }

        std::vector<double> localRHS;
        
        // assemble local load vector of size 3
        if (advection_params.GetAdvectionSolverType() == AdvectionSolverType::TTG2)
        {
            if (step_num == StepNumber::first)
            {
                localRHS = LocalTTG2RhsAssembling(local_node_coords,
                                                  local_u,
                                                  local_m,
                                                  local_m,
                                                  model_params.GetTimeStepHours()*3600.0,
                                                  1);
            }
            else
            {
                localRHS = LocalTTG2RhsAssembling(local_node_coords,
                                                  local_u,
                                                  local_m,
                                                  local_m_half,
                                                  model_params.GetTimeStepHours()*3600.0,
                                                  2);
            }
        }
        else if (advection_params.GetAdvectionSolverType() == AdvectionSolverType::TTG3)
        {
            if (step_num == StepNumber::first)
            {
                localRHS = LocalTTG3RhsAssembling(local_node_coords,
                                                  local_u,
                                                  local_m,
                                                  local_m,
                                                  model_params.GetTimeStepHours()*3600.0,
                                                  1);
            }
            else
            {
                localRHS = LocalTTG3RhsAssembling(local_node_coords,
                                                  local_u,
                                                  local_m,
                                                  local_m_half,
                                                  model_params.GetTimeStepHours()*3600.0,
                                                  2);
            }
        }
        else if (advection_params.GetAdvectionSolverType() == AdvectionSolverType::TTG4)
        {
            if (step_num == StepNumber::first)
            {
                localRHS = LocalTTG4RhsAssembling(local_node_coords,
                                                  local_u,
                                                  local_m,
                                                  local_m,
                                                  model_params.GetTimeStepHours()*3600.0,
                                                  1);
            }
            else
            {
                localRHS = LocalTTG4RhsAssembling(local_node_coords,
                                                  local_u,
                                                  local_m,
                                                  local_m_half,
                                                  model_params.GetTimeStepHours()*3600.0,
                                                  2);
            }
        }
        else
        {
            INMOST_ICE_ERR("Only avalible double step solvers: TTG2, TTG3, TTG4");
        }
            
        // assemble global load vector
        for (unsigned int i = 0; i < 3; ++i)
        {
            if(local_nodes[i].GetStatus() != Element::Ghost)
            {
                RHS[local_nodes[i].Integer(id)] += localRHS[i];
            }
        }
    }
    BARRIER
}

void AdvectionSolver::AssembleRHS()
{
    if (IsAdvectionSolverSingleStep[advection_params.GetAdvectionSolverType()])
    {
        AssembleSingleStepRHS();
    }
    else 
    {
        AssembleDoubleStepRHS(StepNumber::first);
    }
}

void AdvectionSolver::Evaluate()
{
    m_high_tag = ice_mesh.GetMesh()->CreateTag("m_high", DATA_REAL, NODE, NONE, 1);
    if (advection_params.GetIsFct())
    {
        m_low_tag = ice_mesh.GetMesh()->CreateTag("m_low", DATA_REAL, NODE, NONE, 1);
    }
    if (!IsAdvectionSolverSingleStep[advection_params.GetAdvectionSolverType()])
    {
        m_half_tag = ice_mesh.GetMesh()->CreateTag("m_half", DATA_REAL, NODE, NONE, 1);
    }


    INMOST::Tag id = ice_mesh.GetData().NodeIdGlobal;
    INMOST::Sparse::Vector SOL;
    unsigned int idmin = ice_mesh.GetIdIntervalGlobal().IdMin;
    unsigned int idmax = ice_mesh.GetIdIntervalGlobal().IdMax;
    SOL.SetInterval(idmin, idmax);
    Sol.SetMatrix(LHS);
    BARRIER
    Sol.Solve(RHS, SOL);
    
    // update m high value
    for( Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
             nodeit != ice_mesh.GetMesh()->EndNode();
             ++nodeit)
    {
        if(nodeit->GetStatus() != Element::Ghost)
        {
            if (!IsAdvectionSolverSingleStep[advection_params.GetAdvectionSolverType()])
            {
                nodeit->Real(m_half_tag) = SOL[nodeit->Integer(id)];
            }
            else
            {
                nodeit->Real(m_high_tag) = SOL[nodeit->Integer(id)];
            }
        }
    }
    BARRIER
    if (!IsAdvectionSolverSingleStep[advection_params.GetAdvectionSolverType()])
    {
        ice_mesh.GetMesh()->ExchangeData(m_half_tag, NODE, 0);       
    }
    else
    {
        ice_mesh.GetMesh()->ExchangeData(m_high_tag, NODE, 0);
    }
    BARRIER

    // make second step if TTG2, TTG3 or TTG4 is used

    if (!IsAdvectionSolverSingleStep[advection_params.GetAdvectionSolverType()])
    {
        // Reset RHS and SOL vector
        SOL.Clear();
        RHS.Clear();
        SOL.SetInterval(idmin, idmax);
        RHS.SetInterval(idmin, idmax);
        BARRIER
        AssembleDoubleStepRHS(StepNumber::second);
        BARRIER
        Sol.Solve(RHS, SOL);

        // finally update m_high value
        for( Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
             nodeit != ice_mesh.GetMesh()->EndNode();
             ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                nodeit->Real(m_high_tag) = SOL[nodeit->Integer(id)];
            }
        }
        BARRIER
        ice_mesh.GetMesh()->ExchangeData(m_high_tag, NODE, 0);
    }

    // Reset SOL vector
    SOL.Clear();
    SOL.SetInterval(idmin, idmax);
    BARRIER
    
    // Calculate low order solution if FCT is used
    if (advection_params.GetIsFct())
    {
        // make old MASS vector and tmp vectors
        Sparse::Vector MASS, tmp_vec_1, tmp_vec_2, tmp_vec_3;
        MASS.SetInterval(idmin, idmax);
        tmp_vec_1.SetInterval(idmin, idmax);
        tmp_vec_2.SetInterval(idmin, idmax);
        tmp_vec_3.SetInterval(idmin, idmax);

        // fill mass vector
        for( Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
             nodeit != ice_mesh.GetMesh()->EndNode();
             ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                MASS[nodeit->Integer(id)] = nodeit->Real(m_tag);
            }
        }
        BARRIER
        mat_mult_vec(LHS, MASS, tmp_vec_1);
        BARRIER
        mat_mult_vec(LHS_low, MASS, tmp_vec_2);
        BARRIER
        vec_mult_num(tmp_vec_1, (advection_params.GetFctCd() - 1.0), idmin, idmax);
        BARRIER
        vec_mult_num(tmp_vec_2, (1.0 - advection_params.GetFctCd()), idmin, idmax);
        BARRIER
        vec_plus_vec(tmp_vec_1, tmp_vec_2, tmp_vec_3, idmin, idmax);
        BARRIER
        vec_plus_vec(tmp_vec_3, RHS, RHS_low, idmin, idmax);
        BARRIER


        // solve low order system 
        Sol.SetMatrix(LHS_low);
        Sol.Solve(RHS_low, SOL);

        // update m_low value
        for( Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
             nodeit != ice_mesh.GetMesh()->EndNode();
             ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                nodeit->Real(m_low_tag) = SOL[nodeit->Integer(id)];
            }
        }
        BARRIER
        ice_mesh.GetMesh()->ExchangeData(m_low_tag, NODE, 0);

        // clear RHS_low vector
        RHS_low.Clear();
        RHS_low.SetInterval(idmin, idmax); 
    }

    // clear RHS and SOL vector
    SOL.Clear();
    RHS.Clear();
    SOL.SetInterval(idmin, idmax);
    RHS.SetInterval(idmin, idmax);

    // apply FCT (if used) and update m finally
    if (advection_params.GetIsFct())
    {
        FCT_Procedure(ice_mesh, M_C_minus_M_L, m_tag, m_low_tag, m_high_tag, advection_params.GetFctCd());
    }
    else
    {
       for( Mesh::iteratorNode nodeit = ice_mesh.GetMesh()->BeginNode();
             nodeit != ice_mesh.GetMesh()->EndNode();
             ++nodeit)
        {
            if(nodeit->GetStatus() != Element::Ghost)
            {
                nodeit->Real(m_tag) = nodeit->Real(m_high_tag);
            }
        }
        ice_mesh.GetMesh()->ExchangeData(m_tag, NODE, 0); 
    }
    BARRIER

    // clear all temporal data
    ice_mesh.GetMesh()->DeleteTag(m_high_tag);
    if (advection_params.GetIsFct())
    {
        ice_mesh.GetMesh()->DeleteTag(m_low_tag);
    }
    if (!IsAdvectionSolverSingleStep[advection_params.GetAdvectionSolverType()])
    {
        ice_mesh.GetMesh()->DeleteTag(m_half_tag);
    }

    BARRIER

    // log if verbose
    if (is_verbose_advection)
    {
        LogStep();   
    }
    BARRIER
}

void AdvectionSolver::LogInit()
{
    if (ice_mesh.GetMesh()->GetProcessorRank() == 0)
    {
        cout << "=============================================" << endl;
        cout << "Advection Solver Info" << endl;
        cout << "Transported scalar:" << ModelVariableNotationToName[transported_scalar] << endl;
        cout << "Transporting velocity:" << ModelVariableNotationToName[transporting_velocity] << endl;
        cout << "Advection solver type:" << AdvectionSolverTypeToName[advection_params.GetAdvectionSolverType()] << endl;
        cout << "Is FCT Applied:" << advection_params.GetIsFct() << endl;
        if (advection_params.GetIsFct())
        {
            cout << "Is FCT Applied:" << advection_params.GetFctCd() << endl;    
        }
        cout << "=============================================" << endl;
    }
    BARRIER
}

void AdvectionSolver:: LogStep()
{
    double relative_mass_integral = IntegralMass(ice_mesh, m_tag)/init_mass_integral;
    double minimal_mass = min_f(ice_mesh, m_tag);
    double maximal_mass = max_f(ice_mesh, m_tag);

    BARRIER

    if (ice_mesh.GetMesh()->GetProcessorRank() == 0)
    {
        cout << "=============================================" << endl;
        cout << "Current " << ModelVariableNotationToName[transported_scalar] 
             << " integral/Initial " 
             << ModelVariableNotationToName[transported_scalar] 
             << " integral = " << relative_mass_integral << ";" << endl;
        
        cout << "Minimal " << ModelVariableNotationToName[transported_scalar] 
             << " : " << minimal_mass << endl;

        cout << "Maximal " << ModelVariableNotationToName[transported_scalar] 
             << " : " << maximal_mass << endl;
        cout << "=============================================" << endl;
    }

    BARRIER
}