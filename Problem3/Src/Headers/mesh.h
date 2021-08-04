#pragma once
#include "external.h"
#include "enum_classes.h"
#include "config.h"
#include "coords_rotation.h"
#include "mat_vec.h"

struct MeshInfo
{
    int NumberOfTriangles;
    int NumberOfNodes;
    int NumberOfBoundaryNodes;
    int NumberOfNotBoundaryNodes;
    double MinimalTriangleSize;
    double MaximalTriangleSize;
};

struct MeshData
{
    std::map<ModelVariableNotation, INMOST::Tag> NodeData;
    std::map<ModelVariableNotation, INMOST::Tag> TriangleData;
    INMOST::Tag NodeIdGlobal;
    INMOST::Tag NodeIdNoBnd;
    INMOST::Tag IsNodeBnd;
};

struct IdInterval
{
    int IdMin;
    int IdMax;
};

struct TriangleBasisData
{
    INMOST::Tag VecTransMatriciesFromNodalToTriangle;
    INMOST::Tag VecTransMatriciesFromTriangleToNodal;
    INMOST::Tag LocalNodeCoords;
};

struct NodalBasisData
{
    INMOST::Tag LocalNodeCoords;
};

struct LocalBasisData
{
    TriangleBasisData triangle_data;
    NodalBasisData nodal_data;
};

class IceMesh
{
public:
    IceMesh(const MeshParams& mp,
            const OutputParams& op,
            const ModelParams& modp);  
    ~IceMesh();                                                            
    INMOST::Mesh* GetMesh();
    MeshData& GetData();
    MeshInfo& GetMeshInfo();
    std::map<NodeCoordsNotation, INMOST::Tag>& GetCoords();
    IdInterval& GetIdIntervalGlobal();
    IdInterval& GetIdIntervalNoBnd();
    LocalBasisData& GetLocalBasisData();
    void PrintPMF(const std::string& filename) const;     
    void PrintPVTU(const std::string& filename) const;    
    void PrintNETCDF(const std::string& filename) const; // to do 

private:
    void Partition();
    void AssignVariables();
    void AssignCoords();
    void AssignIdIntervals();
    void CalculateMeshInfo();
    void SetBoundaryNodes(bool display);
    void AssembleLocalBasisData(bool display);

    INMOST::Mesh* ice_mesh;
    MeshInfo mesh_info;
    MeshParams mesh_params;
    OutputParams output_info;
    ModelParams model_params;
    MeshData data;
    std::map<NodeCoordsNotation, INMOST::Tag> node_coords;
    IdInterval id_interval_global;
    IdInterval id_interval_no_bnd;
    LocalBasisData local_basis_data;
};