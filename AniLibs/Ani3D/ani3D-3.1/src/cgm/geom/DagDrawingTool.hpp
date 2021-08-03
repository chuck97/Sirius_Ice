
#include "CubitDefines.h"

#include "TopologyEntity.hpp"
//for CubitLink definintion

class RefFace;
class RefEdge;
class RefVertex;
class RefVolume;
template <class X> class DLIList;
class RefEntity;

class RefEntity;
class TopologyEntity;

class TopologyEntity;
class RefVertex;
class RefEdge;
class RefFace;
class RefVolume;
class Body;
class TopologyEntity;

class CUBIT_GEOM_EXPORT DagNodeRow
{
	public:
	
    DagNodeRow( DLIList<TopologyEntity*>& node_list );
		~DagNodeRow();
	
		int length() const { return length_; }
    TopologyEntity*& operator[]( int index );
		
	protected:
    TopologyEntity** array_;
		int length_;
};

class CUBIT_GEOM_EXPORT DagNodeTable
{
	public:
	
		DagNodeTable( int num_rows );
		~DagNodeTable();
		
		int rows() const { return length_; }
		
		DagNodeRow& operator[]( int row );
	
	protected:
		DagNodeRow** array_;
		int length_;
};

class CUBIT_GEOM_EXPORT DagDrawingTool
{
  public:
	
  static DagDrawingTool* instance();
  ~DagDrawingTool();
/*		
		void draw_DAG( DLIList<Body*>&      body_list,           int down );
		void draw_DAG( DLIList<RefVolume*>& volume_list, int up, int down );
		void draw_DAG( DLIList<RefFace*>&   face_list,   int up, int down );
		void draw_DAG( DLIList<RefEdge*>&   edge_list,   int up, int down );
		void draw_DAG( DLIList<RefVertex*>& vertex_list, int up           );
		void draw_DAG( Body*      body,           int down );
		void draw_DAG( RefVolume* volume, int up, int down );
		void draw_DAG( RefFace*   face,   int up, int down );
		void draw_DAG( RefEdge*   edge,   int up, int down );
		void draw_DAG( RefVertex* vertex, int up           );
*/		
    
//    void draw_host_parasite_links( CubitBoolean yes_no );
//    CubitBoolean draw_host_parasite_links() const;
/*		
		//Ways to label RefEntities:
		void label_ref_entity_name() { ref_entity_label_type_ = 0; }
		//Label with RefEntity::entity_name()
		void label_ref_entity_type() { ref_entity_label_type_ = 1; }
		//Label with type (e.g. curve) and id
		void label_ref_entity_short() { ref_entity_label_type_ = 2; }
		//Label with abbreviation and id
		//Body->B, RefVolume->Vm, RefFace->F, RefEdge->E, RefVertex->Vx
		void label_ref_entity_id() { ref_entity_label_type_ = 3; }
		//Label with id only
		
		//Ways to label entities other than RefEntities:
		void label_other_type() { other_entity_label_type_ = 0; }
		//Label with type (e.g. CoEdge)
		void label_other_short() { other_entity_label_type_ = 1; }
		//Label with abbreviations
		//CoVolume->CVm, CoFace->CF, CoEdge->CE, CoVertex->CVx
		//Shell->Sh, Loop->L, Chain->Ch 
		void label_other_none() { other_entity_label_type_ = 2; }
		//Do not show any labels for other entities.
		
		void primary_link_color( int color );
		int primary_link_color() const;
		//Color to draw primary links in.

		void secondary_link_color( int color );
		int secondary_link_color() const;
		//Color to draw secondary (directional) links in.
		
		void bad_link_color( int color );
		int bad_link_color() const;
		//Color to draw broken links in 
		//(hopefully you won't see any of these)
    
    void host_para_link_color( int color );
    int host_para_link_color() const;
		
		void node_color( int color );
		int node_color() const;
		//Color to draw nodes in,
		//(This currently has no effect.)
		
		void label_color( int color );
		int label_color() const;
		//Color to draw node labels in
		
		void background_color( int color );
		int background_color() const;
		//Color of the background
		
		int window_id() const;
		//The last window drawn in my DagDrawingTool
*/	
  
  //! Print the DAG of the specified surfaces.
  void printDag(DLIList<RefFace*> &face_list, int depth);

  //! Print the DAG of the specified body.
  void printDag(DLIList<Body*> &body_list, int depth);

  //! Print the DAG of the entities. 
  void printDag(DLIList<TopologyEntity*> &entity_list, int depth);
  
  //! Print the DAG of the entity. 
  void printDag( TopologyEntity* ME_ptr, int depth );

	protected:
	
		DagDrawingTool();
		static DagDrawingTool* instance_;
	
		void draw_DAG( DLIList<RefEntity*>& enity_list, int up, int down );
    void draw_DAG( DLIList<TopologyEntity*>& node_list, int up, int down );

//		void draw_table( );
		
    //const char* name( TopologyEntity* node_ptr );
		//const char* name( int row, int index );
		//int name_width( int row );
    //int name_width( TopologyEntity* node_ptr );
		
		//direction-indenpendent query stuff:
    void get_relatives( TopologyEntity* node_ptr,
                        DLIList<TopologyEntity*>& result_set,
		                    int direction );
    void get_relatives( DLIList<TopologyEntity*>& source_set,
                        DLIList<TopologyEntity*>& result_set,
												int direction );
    int link_type( TopologyEntity* source_node,
                   TopologyEntity* target_node,
		               int direction );
		int link_type( int source_row, int source_index,
		               int target_row, int target_index );
		               
		void position( int row, int index, float& x, float& y );
                /*
		void draw_node( int row, int index );
		void draw_link( int source_row, int source_index,
		                int target_row, int target_index );
   void draw_host_link( int row, int source_index, int target_index );


		//The only methods that interact with drawing tool
		void draw_link( float x1, float y1, float x2, float y2, int type );
    void draw_host_link( float x1, float x2, float y, int top );
    void draw_node( float x, float y, TopologyEntity* node );
		void draw_ellipse( float cx, float cy, float rx, float ry, int color );
    void draw_arc_link( float x1, float y1, float x2, float y2 );
		void initialize_graphics();
		void finalize_graphics();
		float text_height() const;
		float arrow_size() const;
                */
	
	private:
	
  void print_dag( TopologyEntity* node_ptr, int depth );
  
  void print_dag( TopologyEntity* node_ptr, int indent, int depth );
  
  void print_node( TopologyEntity* node_ptr, int indent, const char* prefix);
  
  int get_id( TopologyEntity* node_ptr );
  
		void one_sort_pass( int start_row, int end_row, float* temp_array );
    
    void sort_row( int row, int rel_row /*1 or -1*/ );
	
    /*
		int window_id_;
		int active_window_;

		int prim_link_color_;
		int sec_link_color_;
		int bad_link_color_;
		int node_color_;
		int label_color_;
		int background_color_;
    int host_para_color_;

		CubitLink link_type_;
    CubitBoolean draw_host_para_links_;
		
		int ref_entity_label_type_;
		// 0  use RefEntity:entity_name()
		// 1  use type and id
		// 2  use abbreviated type and id
		// 3  use only id
		
		int other_entity_label_type_;
		// 0  use type name
		// 1  use abbreviated type name
		// 2  no label
              */
		
		//Working data
		DagNodeTable* node_table;
};

class DagListTool
{
public:
	
};
  
