//-------------------------------------------------------------------------
// Copyright Notice
//
// Copyright (c) 1996 
// by Malcolm J. Panthaki, DBA, and the University of New Mexico.
//-------------------------------------------------------------------------
//
//-------------------------------------------------------------------------
// Filename      : ModelQueryEngine.hpp
//
// Purpose       : This class provides the interface to query the model.
//
// Special Notes : This is a singleton class. 
//
//                 The primary objects that this class deals with, at the
//                 interface, are TopologyEntity's.  These are entities
//                 that are represented as nodes in the DAG.  Each node of
//		   the DAG contains a pointer to a TopologyEntity.
//
//                 All of the query functions rely on the "type" of 
//                 objects. Some of them rely on the relationships
//                 between different "type"s of objects. Currently all
//                 such relationships, i.e. the entity relation diagram 
//                 (ERD) of the different TopologyEntity classes is maintained
//                 by a class called ModelERD.
//
// Creator       : Xuechen Liu 
//
// Creation Date : 06/08/96
//
// Owner         : Malcolm J. Panthaki
//-------------------------------------------------------------------------

#ifndef MODEL_QUERY_ENGINE_HPP
#define MODEL_QUERY_ENGINE_HPP


// ********** BEGIN STANDARD INCLUDES      **********
// ********** END STANDARD INCLUDES        **********

// ********** BEGIN MOTIF INCLUDES         **********
// ********** END MOTIF INCLUDES           **********

// ********** BEGIN OPEN INVENTOR INCLUDES **********
// ********** END OPEN INVENTOR INCLUDES   **********

// ********** BEGIN ACIS INCLUDES          **********
// ********** END ACIS INCLUDES            **********

// ********** BEGIN CUBIT INCLUDES         **********

#include "CubitDefines.h"
#include "TopologyEntity.hpp"
#include "DLIList.hpp"

// ********** END CUBIT INCLUDES           **********

// ********** BEGIN MACROS DEFINITIONS     **********
// ********** END MACROS DEFINITIONS       **********

// ********** BEGIN FORWARD DECLARATIONS   **********

class TopologyEntity;

// ********** END FORWARD DECLARATIONS     **********

// ********** BEGIN ENUM DEFINITIONS       **********

enum CubitSearchDirection 
{ 
   CUBIT_SEARCH_NO_DIRECTION = -1,
   CUBIT_SEARCH_PARENT = 0, 
   CUBIT_SEARCH_CHILD = 1 
};
// ********** END ENUM DEFINITIONS         **********


class CUBIT_GEOM_EXPORT ModelQueryEngine 
{
   public:
      
      static ModelQueryEngine* instance();
      //R ModelQueryEngine* 
      //R- A pointer to the only instance of this class.
      //- This function controls access and creation of the sole instance 
      //- of this class. It ensures that only one instance can ever get 
      //- created. It returns a pointer to the only instance of the class.

      ~ModelQueryEngine() ;
      //- Destructor

      static void delete_instance()
      {
        if(instance_)
          delete instance_;
        instance_ = NULL;
      }

      //HEADER- Query functions on single source objects
      
      CubitStatus query_model( TopologyEntity & source_object,
                               DagType target_type,
                               DLIList<TopologyEntity*>& result_set );

      CubitStatus query_model_and_append( TopologyEntity& source_object,
                                          DagType target_type,
                                          DLIList<TopologyEntity*>& result_set );
      //R CubitStatus
      //R- CUBIT_SUCCESS/FAILURE
      //I sourceObject
      //I- A reference to the ModEnt on which query is to be done.
      //I targetType
      //I- The type of the ModEnts to query.
      //I- It can be any kind of data type: root base type, intermediate
      //I- base type, or leaf type.
      //O resultModEntSet
      //O- Reference to a set of ModEnts where output of the query will be
      //O- put/appended.
      //- This function queries the given object, "sourceObject", for all 
      //- the ModEnts of the given type, "targetType" that it is related to and 
      //- puts/apppends the result of the query in "resultModEntSet".
      //- The return value is CUBIT_SUCCESS if the source object has either
      //- a parent-child or child-parent relationship with objects of the
      //- given type, CUBIT_FAILURE otherwise.

      


      //HEADER- Query functions on sets of source objects

      CubitStatus query_model( DLIList<TopologyEntity*>& source_set,
                               DagType target_type,
                               DLIList<TopologyEntity*>& result_set );
      
      CubitStatus query_model_and_append( DLIList<TopologyEntity*>& source_set,
                                          DagType target_type,
                                          DLIList<TopologyEntity*>& result_set );
      //R CubitStatus
      //R- CUBIT_SUCCESS/FAILURE
      //I sourceObjectSet
      //I- A reference to a set of ModEnts on which query is to be done.
      //I targetType
      //I- The type of the ModEnts to query.
      //I- It can be any kind of data type: root base type, intermediate
      //I- base type, or leaf type.
      //O resultModEntSet
      //O- Reference to a set of ModEnts where output of the query will be
      //O- put/appended.
      //- This function queries the given objects, "sourceObjectSet", for all 
      //- the ModEnts of the given type, "targetType" that they are related to and 
      //- puts/apppends the result of the query in "resultModEntSet".
      //- The return value is CUBIT_SUCCESS if each of the source objects 
      //- has either a parent-child or child-parent relationship with 
      //- objects of the given type, CUBIT_FAILURE otherwise.


      bool encountered( TopologyEntity* );
        //- Mark node as encountered if it was not already encountered.
        //- Return the previous value of the encountered flag.
  
      class BeginQuery {
        public:
          inline BeginQuery()
            { ModelQueryEngine::instance()->inc_query_call_stack(); }
          inline ~BeginQuery()
            { ModelQueryEngine::instance()->dec_query_call_stack(); }
          inline void* operator new(size_t /*size*/) throw()
            { assert(0); return (void*)0; }
      };

   protected:
   
      CubitStatus query_append_children ( TopologyEntity& source_object,
                                          DagType child_type,
                                          DLIList<TopologyEntity*>& result_set );
      CubitStatus query_append_parents ( TopologyEntity& source_object,
                                         DagType parent_type,
                                         DLIList<TopologyEntity*>& result_set );
                                          
   
      
      friend class ModelQueryEngine::BeginQuery;

      void inc_query_call_stack();
      void dec_query_call_stack();
      
      int queryCallStackDepth;
      
      DLIList<TopologyEntity*> encounteredSet;
     // A set of marked ModelEntities.
    
      DLIList<TopologyEntity*> intermediateNodeSets[2];

   private:

      ModelQueryEngine() ;
      //- Constructor. (Not callable by user code. Class is constructed
      //- by the "instance()" member function.

      static ModelQueryEngine* instance_;
      //- The static pointer to unique instance of this class.
};

// ********** BEGIN HELPER CLASSES         **********
// ********** END   HELPER CLASSES         **********

// ********** BEGIN INLINE FUNCTIONS       **********

// ********** END INLINE FUNCTIONS         **********
 
// ********** BEGIN FRIEND FUNCTIONS       **********
// ********** END FRIEND FUNCTIONS         **********
 
// ********** BEGIN EXTERN FUNCTIONS       **********
// ********** END EXTERN FUNCTIONS         **********
 
#endif

