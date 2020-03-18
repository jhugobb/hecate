#ifndef TWODGRID_H
#define TWODGRID_H

#include <vector>
#include <unordered_map>
#include <list>
#include <glm/glm.hpp>

#include "Triangle.h"
#include "../TriangleMesh.h"
#include "Geometry.h"

typedef std::pair<unsigned int, unsigned int> p_yz; 


class BinTreeNode {
  public:
    BinTreeNode() {};
    ~BinTreeNode() {}; 
    glm::vec3 min_point;
    glm::vec3 max_point;
    int representative;
    bool is_gray = false;
};

class BinTree {
  public:
    BinTree() {};
    std::list<BinTreeNode*> nodes;
    TriangleMesh* mesh_;
    Geo::BBox bbox;
    void subdivide_recursive();
};


class node {
  public:
    node() {};
    // bool is_leaf = true; 
    // node* children[4];
    glm::vec2 min_point;
    glm::vec2 max_point;
    std::vector<int> members;

    std::list<BinTreeNode*> build_bin_tree(TriangleMesh* mesh, Geo::BBox space, double min_node_size);
};

// Quadtree will manage node*
class TwoDGrid {

  public:
    TwoDGrid();
    TwoDGrid(TriangleMesh* mesh, int size, Geo::BBox space);
    ~TwoDGrid();
    void insert(int t);
    node* query(glm::vec2 coords);

  private:
    std::unordered_map<uint, node*> grid;
    TriangleMesh* m;
    int num_nodes;
    double node_size;
    Geo::BBox space;
    // node* insert_recursive(int t, node* cell, TriangleMesh* mesh, int depth);
    // node* query_recursive(node* cell, glm::vec2 coords);
};

#endif