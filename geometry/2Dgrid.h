#ifndef TWODGRID_H
#define TWODGRID_H

#include <vector>
#include <unordered_map>
#include <glm/glm.hpp>

#include "Triangle.h"
#include "../TriangleMesh.h"
#include "Geometry.h"

typedef std::pair<unsigned int, unsigned int> p_yz; 

class node {
  public:
    node() {};
    // bool is_leaf = true; 
    // node* children[4];
    glm::vec2 min_point;
    glm::vec2 max_point;
    std::vector<int> members;
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