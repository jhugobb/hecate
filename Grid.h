#ifndef GRID_H
#define GRID_H

#include <map>
#include <set>
#include "geometry/Geometry.h"
#include "geometry/2Dgrid.h"

class Grid {

  enum VoxelColor {
    WHITE,
    BLACK,
    GRAY
  };

  struct Voxel {
    VoxelColor color = VoxelColor::WHITE;
    double x;
    glm::vec3 normal;
    
    bool operator <(const Voxel other ) const {
      return x < other.x;
    }
  };

  public:
    Grid() {};
    Grid(unsigned int size, Geo::BBox space_);
    ~Grid() {}; // TODO
    
    /**
     * Voxelizes the model
     *
     * @param mesh Triangle Mesh of the model to be voxelized
     * @param qt 2Dgrid of the model that subdivided the space
     * @param useNaive boolean that dictates if the naive approach should be used
     * @param calculate_b_w boolean that dictates if black and white coloring should be calculated
     * @param filename name of the file where to save the voxelization
     */ 
    void colorGrid(TriangleMesh* mesh, TwoDGrid* qt, bool useNaive, bool useBox, bool calculate_b_w, std::string filename);


    /**
     * Creates a PLY model with the voxelization
     *
     * @param filename The name of the file to be created
     */ 
    void writePLY(std::string filename);

  private:
    // Number of nodes in a dimension
    uint size_;
    
    // Size of a grid node in any direction
    double node_size;

    // Bounding Box of the model
    Geo::BBox space;

    // Vector of voxels
    // TODO change to only have a single row of voxels at any given point in time
    std::vector<Voxel> elements;

    TriangleMesh* mesh_;

    bool testVoxelGray_Naive(int &representative, glm::vec3 voxel_min, glm::vec3 voxel_max, std::vector<int> candidates);
    
    void testVoxelGray_Box(std::vector<Voxel> &voxels, std::vector<int> candidates, glm::vec2 row_coords);

};
#endif