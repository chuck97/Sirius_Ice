/**
 * \class CGMParallelComm
 * \brief Parallel communications in CGM
 * \author Hong-Jun Kim, copied from MBParallelComm.hpp
 *
 *  This class implements methods to communicate geometry between processors
 *
 */

#ifndef CGM_PARALLEL_COMM_HPP
#define CGM_PARALLEL_COMM_HPP

#include "CGMProcConfig.hpp"
#include "CGMmpi.h"
#include "DLIList.hpp"
#include <vector>

class GeometryQueryTool;
class RefEntity;

class CGMParallelComm 
{
public:

    //! constructor
  CGMParallelComm(MPI_Comm comm = MPI_COMM_WORLD);
  
  //! constructor taking buffer, for testing
  CGMParallelComm(std::vector<unsigned char> &tmp_buff,
		  MPI_Comm comm = MPI_COMM_WORLD);


  
  //! destructor
  ~CGMParallelComm();
	
  //! return partition ref entity list
  DLIList<RefEntity*> &partition_surf_list() {return partitioningSurfList;}
  const DLIList<RefEntity*> &partition_surf_list() const {return partitioningSurfList;}
  DLIList<RefEntity*> &partition_body_list() {return partitioningBodyList;}
  const DLIList<RefEntity*> &partition_body_list() const {return partitioningBodyList;}
  
  //! Get proc config for this communication object
  const CGMProcConfig &proc_config() const {return procConfig;}
  
  //! Get proc config for this communication object
  CGMProcConfig &proc_config() {return procConfig;}

  CubitStatus broadcast_entities(const unsigned int from_proc,
				 DLIList<RefEntity*> &ref_entity_list);

  CubitStatus scatter_entities(const unsigned int from_proc,
			       DLIList<RefEntity*> &ref_entity_list);

  CubitStatus write_buffer(DLIList<RefEntity*> &ref_entity_list,
			  char* pBuffer,
			  int& n_buffer_size,
			  bool b_export_buffer);

  CubitStatus read_buffer(DLIList<RefEntity*> &ref_entity_list,
			    const char* pBuffer,
			    const int n_buffer_size);
    
  CubitStatus bcast_buffer(const unsigned int from_proc);

  CubitStatus append_to_buffer(DLIList<RefEntity*> &ref_entity_list,
			       int add_size);

private:  


  CubitStatus check_size(int& target_size, const CubitBoolean keep = CUBIT_FALSE);
  
    //! CGM query tool interface associated with this writer
  GeometryQueryTool *gqt;

    //! Proc config object, keeps info on parallel stuff
  CGMProcConfig procConfig;
  
    //! data buffer used to communicate
  std::vector<unsigned char> myBuffer;

  char* m_pBuffer;
  
  int m_nBufferSize;
  
  int m_currentPosition;
  
  DLIList<RefEntity*> partitioningSurfList; // ref entity list containing all parts
  DLIList<RefEntity*> partitioningBodyList; // ref entity list containing all parts

};

#endif
