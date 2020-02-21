#include <algorithm>

#include "Grid.h"
#include "geometry/Geometry.h"

Grid::Grid(unsigned int size) {
  size_ = size;
  for (unsigned int x = 0; x < size; x++) {
    for (unsigned int y = 0; y < size; y++) {
      std::set<Voxel> s;
      for (unsigned int z = 0; z < size; z++) {
        Voxel voxel;
        voxel.z = z;
        s.insert(voxel);
      }
      elements.emplace(std::make_pair(x,y), s);
    }
  }
}

void Grid::colorGrid(TriangleMesh* mesh, Quadtree qt) {
  std::vector<Triangle> triangles = mesh->getTriangles();
  std::vector<glm::vec3> vertices = mesh->getVertices();

  for (unsigned int x = 0; x < size_; x++) {
    for (unsigned int y = 0; y < size_; y++) {
      node* quad_node = qt.query(glm::vec2(x+0.5, y+0.5));
      std::vector<double> intersect_zs;

      for (unsigned int tri_idx : quad_node->members) {
        Triangle t = triangles[tri_idx];
        glm::vec3 intersection_point;
        glm::vec3 v1 = vertices[t.getV1()];
        glm::vec3 v2 = vertices[t.getV2()];
        glm::vec3 v3 = vertices[t.getV3()];

        // Ray in the middle of the node pointing upward
        bool intersects = Geo::rayIntersectsTriangle(glm::vec3(x+0.5, y+0.5,-0.5), glm::vec3(0,0,1), v1, v2, v3, intersection_point);

        if (intersects) {
          intersect_zs.push_back(intersection_point.z);
        }
      }
      std::sort(intersect_zs.begin(), intersect_zs.end());
      if (intersect_zs.size() != 0) {
        unsigned int z_idx = 0;
        std::set<Voxel>::iterator it;
        for (it = elements[std::make_pair(x,y)].begin(); it != elements[std::make_pair(x,y)].end(); ++it) {
          Voxel v = *it;
          while (z_idx <= intersect_zs.size() || v.z > intersect_zs[z_idx]) {
            z_idx++;
          }
          if (z_idx != 0) { 
            bool intersects_node = false;
            for (unsigned int tri_idx : quad_node->members) {
              Triangle t = triangles[tri_idx];
              intersects_node = Geo::testBoxTriangle(mesh, t, glm::vec3(x,y,v.z), glm::vec3(x+1, y+1, v.z+1));
              if (intersects_node) break;
            }
            if (intersects_node) v.color = VoxelColor::GRAY;
            else if (z_idx % 2 == 1) v.color = VoxelColor::BLACK;
          }
        }
      }
    }
  }
}