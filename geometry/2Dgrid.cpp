#include "2Dgrid.h"

#include <iostream>

TwoDGrid::TwoDGrid() {}

TwoDGrid::TwoDGrid(TriangleMesh* m_, int size, Geo::BBox space_) {
  glm::vec3 center_box = (space_.maxPoint + space_.minPoint)/2.0f;
  glm::vec3 size_box = space_.maxPoint - space_.minPoint;
  float largest = max(size_box.x, max(size_box.y, size_box.z));
  space.minPoint = center_box - glm::vec3(largest, largest, largest)/2.0f;
  space.maxPoint = center_box + glm::vec3(largest, largest, largest)/2.0f;
  node_size = largest/size;
  num_nodes = size;
  m = m_;
  grid = std::unordered_map<uint, node*>();
}

TwoDGrid::~TwoDGrid() {}

void TwoDGrid::insert(int t) {
  std::vector<Triangle> tris = m->getTriangles();
  std::vector<glm::vec3> verts = m->getVertices();

  Triangle tri = tris[t];

  glm::vec3 v1 = verts[tri.getV1()];
  glm::vec3 v2 = verts[tri.getV2()];
  glm::vec3 v3 = verts[tri.getV3()];

  Geo::BBox tri_box;

  tri_box.addPoint(v1);
  tri_box.addPoint(v2);
  tri_box.addPoint(v3);

  glm::vec3 min_box_tri = tri_box.minPoint;
  glm::vec3 max_box_tri = tri_box.maxPoint;

  // translate to origin
  min_box_tri -= space.minPoint;
  max_box_tri -= space.minPoint;

  int initial_y_grid = std::max(min_box_tri.y / node_size - 1, 0.0);
  int final_y_grid = max_box_tri.y / node_size;

  int initial_z_grid = std::max(min_box_tri.z / node_size - 1, 0.0);
  int final_z_grid = max_box_tri.z / node_size;

  bool inserted = false;
    // x * height * depth + y * depth + z
  for (int y = initial_y_grid; y < final_y_grid + 1 && y < num_nodes; y++) {
    for (int z = initial_z_grid; z < final_z_grid + 1 && z < num_nodes; z++) {

      glm::vec2 min_box = glm::vec2(space.minPoint.y, space.minPoint.z) + glm::vec2(y*node_size,z*node_size);
      glm::vec2 max_box = glm::vec2(space.minPoint.y, space.minPoint.z) + glm::vec2((y+1)*node_size,(z+1)*node_size);

      if (Geo::testQuadTriangle(m, tri, min_box, max_box)) {
        inserted = true;
        uint key = y * num_nodes + z;
        if (grid.find(key) != grid.end()) {
          grid[key]->members.push_back(t);
        } else {
          node* n = new node();
          n->members.push_back(t);
          grid.emplace(key, n);
        }
      }
    }
  }

  // std::cout << inserted << std::endl;
}

node* TwoDGrid::query(glm::vec2 coords) {
  int y = floor(coords.x);
  int z = floor(coords.y);

  uint key = y * num_nodes + z;

  if (grid.find(key) != grid.end()) {
    return grid[key];
  } 

  return NULL;
}

