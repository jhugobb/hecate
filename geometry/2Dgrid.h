#ifndef TWODGRID_H
#define TWODGRID_H

#include <vector>
#include <unordered_map>
#include <list>
#include <set>
#include <glm/glm.hpp>

#include "Triangle.h"
#include "../TriangleMesh.h"
#include "Geometry.h"

typedef std::pair<unsigned int, unsigned int> p_yz; 


class BinTreeNode {
  public:
    BinTreeNode() {};
    ~BinTreeNode() {
      triangles.clear();
    }; 
    glm::vec3 min_point;
    glm::vec3 max_point;
    std::multiset<Triangle> triangles;
    // Index of the triangle that represents this cell
    int representative;
    bool is_gray = false;
};

class node {
  public:
    node() {};
    ~node() {
      members.clear();
    }
    glm::vec2 min_point;
    glm::vec2 max_point;
    std::vector<int> members;
    std::list<BinTreeNode*> bin_tree;

    /**
     * Creates the X-Bin-tree of the node
     *
     * @param mesh The triangular mesh of the model
     * @param space The bounding box of the model
     * @param min_node_size the size of a voxel in any direction
     */  
    void build_bin_tree(TriangleMesh* mesh, Geo::BBox space, double min_node_size);
};

// Quadtree will manage node*
class TwoDGrid {

  public:
    TwoDGrid();
    TwoDGrid(TriangleMesh* mesh, int size, Geo::BBox space);
    ~TwoDGrid();

    /**
     * Inserts a triangle in the 2D grid
     *
     * @param t index of the triangle in the TriangleMesh triangles array
     */ 
    void insert(int t);

    /**
     * Finds the 2Dgrid node that corresponds to a voxel
     *
     * @param coords The coordinates of the voxel that corresponds to an X-row of the 2D grid
     * @return the node of the 2D grid that encompases the voxel
     */ 
    node* query(glm::vec2 coords);

    /**
     * Creates the X-Bin-trees of all the nodes
     */ 
    void buildBinTrees();

  private:
    // Array of nodes of the 2Dgrid
    node** grid;

    // Triangle mesh pointer of the model
    TriangleMesh* m;

    // Triangles of the mesh
    std::vector<Triangle*> triangles;

    // Number of nodes in a dimension
    int num_nodes;

    // Size of a 2Dgrid node in any direction
    double node_size;

    // Bounding Box of the model
    Geo::BBox space;
};

#endif