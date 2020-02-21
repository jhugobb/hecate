#ifndef QUAD_H
#define QUAD_H

#include <vector>
#include <glm/glm.hpp>

#include "Triangle.h"
#include "../TriangleMesh.h"

class node {
  public:
    node() {};
    bool is_leaf = true; 
    node* children[4];
    glm::vec2 min_point;
    glm::vec2 max_point;
    std::vector<int> members;
};

// Quadtree will manage node*
class Quadtree {

  public:
    Quadtree();
    Quadtree(TriangleMesh* mesh);
    ~Quadtree();
    void insert(int t);
    node* query(glm::vec2 coords);

  private:
    node* root = NULL;
    TriangleMesh* m;
    node* insert_recursive(int t, node* cell, TriangleMesh* mesh, int depth);
    node* query_recursive(node* cell, glm::vec2 coords);
};

#endif