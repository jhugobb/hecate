#ifndef GRID_H
#define GRID_H

#include <map>
#include <set>

#include "geometry/Quadtree.h"

class Grid {

  enum VoxelColor {
    WHITE,
    BLACK,
    GRAY
  };

  struct Voxel {
    VoxelColor color = VoxelColor::WHITE;
    unsigned int z;
    
    bool operator <(const Voxel& other ) const {
      return z < other.z;
    }
  };

  public:
    Grid() {};
    Grid(unsigned int size);
    ~Grid() {}; // TODO
    void colorGrid(TriangleMesh* mesh, Quadtree qt);

  private:
    unsigned int size_;
    std::map<std::pair<unsigned int, unsigned int>, std::set<Voxel>> elements;

};
#endif