
#ifndef GSTREE_HPP
#define GSTREE_HPP

#include "GridSearchTreeNode.hpp"

#include <map>

typedef std::map< GridSearchTreeNode* , int, GridSearchTreeNode::GSTNodeComparator > gmap;

class GridSearchTree {

private:
  double epsilon;
  gmap nodemap;
  gmap::iterator pos;
public:

  GridSearchTree(double tolerance) {
    epsilon = tolerance;
  }

  ~GridSearchTree() {}
  

  CubitPoint * fix (CubitPoint * data);

};

#endif


