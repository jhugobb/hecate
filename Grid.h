#ifndef GRID_H
#define GRID_H

#include <map>
#include <set>
#include "geometry/Geometry.h"
#include "geometry/Quadtree.h"

class Grid {

  enum VoxelColor {
    WHITE,
    BLACK,
    GRAY
  };

  struct Voxel {
    VoxelColor color = VoxelColor::WHITE;
    float z;
    glm::vec3 normal;
    
    bool operator <(const Voxel *other ) const {
      return z < other->z;
    }
  };

  public:
    Grid() {};
    Grid(unsigned int size, Geo::BBox space_);
    ~Grid() {}; // TODO
    void colorGrid(TriangleMesh* mesh, Quadtree qt);
    void writePLY(std::string filename);

  private:
    unsigned int size_;
    float node_size;
    Geo::BBox space;
    std::map<std::pair<unsigned int, unsigned int>, std::set<Voxel*>> elements;

};
#endif