#include <algorithm>

#include <memory>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>

#include "Grid.h"
#include "geometry/Geometry.h"


// Necessary to Write PLY
template<typename ... Args>
std::string string_format( const std::string& format, Args ... args )
{
    using namespace std;
    size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    unique_ptr<char[]> buf( new char[ size ] ); 
    snprintf( buf.get(), size, format.c_str(), args ... );
    return string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

Grid::Grid(unsigned int size, Geo::BBox space_) {
  size_ = size;
  space = space_;

  // Make voxels cubic
  glm::vec3 center_box = (space.maxPoint + space.minPoint)/2.0f;
  glm::vec3 size_box = space.maxPoint - space.minPoint;
  float largest = max(size_box.x, max(size_box.y, size_box.z));
  space.minPoint = center_box - glm::vec3(largest, largest, largest)/2.0f;
  space.maxPoint = center_box + glm::vec3(largest, largest, largest)/2.0f;
  node_size = largest/size;
  // elements = std::vector<Voxel>();
}

bool Grid::testVoxelGray_Naive(int &representative, glm::vec3 voxel_min, glm::vec3 voxel_max, std::vector<int> candidates) {
  // Naive approach: for each triangle in the node, check if intersects the voxel
  std::vector<Triangle> triangles = mesh_->getTriangles();
  for (int tri_idx : candidates) {
    if (Geo::testBoxTriangle(mesh_, triangles[tri_idx], voxel_min, voxel_max)) {
      representative = tri_idx;
      return true;
    }
  }
  return false;
}

void Grid::testVoxelGray_Box(std::vector<Voxel> &voxels, std::vector<int> candidates, glm::vec2 row_coords) {
  // Box approach: for each triangle in the row, check which voxels it intersects
  std::vector<Triangle> triangles = mesh_->getTriangles();
  std::vector<glm::vec3> verts = mesh_->getVertices();

  for (int tri_idx : candidates) {

    Triangle tri = triangles[tri_idx];

    glm::vec3 v1 = verts[tri.getV1()];
    glm::vec3 v2 = verts[tri.getV2()];
    glm::vec3 v3 = verts[tri.getV3()];

    Geo::BBox tri_box;

    tri_box.addPoint(v1);
    tri_box.addPoint(v2);
    tri_box.addPoint(v3);

    // Save Bbox of triangle
    glm::vec3 min_box_tri = tri_box.minPoint;
    glm::vec3 max_box_tri = tri_box.maxPoint;

    // translate to origin
    min_box_tri -= space.minPoint;
    max_box_tri -= space.minPoint;

    int initial_x_grid = std::max(min_box_tri.x / node_size - 1, 0.0);
    int final_x_grid = max_box_tri.x / node_size;

    glm::vec3 min_box, max_box;
    for (int x = initial_x_grid; x < final_x_grid + 1 && x < size_; x++) {
      
      if (voxels[x].color == VoxelColor::GRAY) continue;

      min_box = glm::vec3(voxels[x].x, row_coords.x, row_coords.y);
      max_box = glm::vec3(voxels[x].x + node_size, row_coords.x + node_size, row_coords.y + node_size);

      if (Geo::testBoxTriangle(mesh_, triangles[tri_idx], min_box, max_box)) {
        voxels[x].color = VoxelColor::GRAY;
        voxels[x].normal = glm::normalize(glm::cross(v3-v2, v3-v1));
      }
    }
  }
}


void Grid::colorGrid(TriangleMesh* mesh, TwoDGrid* qt, bool useNaive, bool useBox, bool calculate_b_w, std::string filename) {
  std::vector<Triangle> triangles = mesh->getTriangles();
  std::vector<glm::vec3> vertices = mesh->getVertices();
  mesh_ = mesh;
  std::ofstream out_fobj(filename);
  out_fobj<< "ply\r\n"
          << "format ascii 1.0\r\n";
  auto file_verts_line = out_fobj.tellp();

  out_fobj<< string_format("element vertex %20d\r\n", 0)
          << "property float x\r\n"
          << "property float y\r\n"
          << "property float z\r\n"
          << "property uchar red\r\n"
          << "property uchar green\r\n"
          << "property uchar blue\r\n"
          << "end_header\r\n";

  uint num_points = 0;

  #pragma omp parallel for
  for (unsigned int y = 0; y < size_; y++) {
    #pragma omp parallel for
    for (unsigned int z = 0; z < size_; z++) {
      glm::vec2 coords;
      std::vector<Voxel> voxels;

      // Coords of voxel in the 2D grid dimensions
      coords.x = space.minPoint.y + y * node_size;
      coords.y = space.minPoint.z + z * node_size;

      // 2Dgrid node that has the triangles for this y z 
      node* quad_node = qt->query(glm::vec2(y,z));
      
      std::list<BinTreeNode*>::iterator list_it;
      std::list<BinTreeNode*> nods;

      // If the node has triangles
      if (quad_node->members.size() != 0) {
        if (!useNaive && !useBox) {
          nods = quad_node->bin_tree;
          list_it = nods.begin();
        }
        if (useBox) {
          for (uint x = 0; x < size_; x++) {
            Voxel v = Voxel();
            v.x = space.minPoint.x + x * node_size;
            voxels.push_back(v);
          }
          testVoxelGray_Box(voxels, quad_node->members, coords);
          for (uint x = 0; x < size_; x++) {
            if (voxels[x].color != VoxelColor::GRAY) continue;
            glm::vec3 n;
            double x_coord = space.minPoint.x + x * node_size;
            unsigned char r, g, b;
            n.x = static_cast<uint8_t>((voxels[x].normal.x * 0.5 + 0.5) * 255);
            n.y = static_cast<uint8_t>((voxels[x].normal.y * 0.5 + 0.5) * 255);
            n.z = static_cast<uint8_t>((voxels[x].normal.z * 0.5 + 0.5) * 255);
            r = glm::clamp(2*abs(int(n.x) - 128), 0, 255);
            g = glm::clamp(2*abs(int(n.y) - 128), 0, 255);
            b = glm::clamp(2*abs(int(n.z) - 128), 0, 255);
            #pragma omp critical 
            {
              out_fobj 
                      << x_coord << " "
                      << space.minPoint.y + y*node_size + node_size/2.f << " "
                      << space.minPoint.z + z*node_size + node_size/2.f << " "
                      << static_cast<int>(r) << " "
                      << static_cast<int>(g) << " "
                      << static_cast<int>(b)
                      <<'\n';
              num_points++;
            }
          }
        } else {
          for (uint x = 0; x < size_; x++) {
            Voxel v;
            // coordinate in world space
            double x_coord = space.minPoint.x + x * node_size;
            bool intersects_node = false;
            int representative = -1;
            if (useNaive) {
              // Naive approach
              intersects_node = testVoxelGray_Naive(representative, 
                                                    glm::vec3(x_coord, coords.x, coords.y),
                                                    glm::vec3(x_coord + node_size, coords.x + node_size, coords.y + node_size),
                                                    quad_node->members);
            } else {
              // 2D Grid Bin tree approach: Follow the bin tree until finding the node that represents this voxel
              while ((*list_it)->max_point.x <= x_coord+node_size){
                ++list_it;
              }
              intersects_node = (*list_it)->is_gray;
              representative = (*list_it)->representative;
            }
            if (intersects_node){
              v.color = VoxelColor::GRAY;
              glm::vec3 normal;
              Triangle t = triangles[representative];
              glm::vec3 v1 = vertices[t.getV1()];
              glm::vec3 v2 = vertices[t.getV2()];
              glm::vec3 v3 = vertices[t.getV3()];

              normal = glm::normalize(glm::cross(v3-v2, v3-v1));

              glm::vec3 n;
              unsigned char r, g, b;
              n.x = static_cast<uint8_t>((normal.x * 0.5 + 0.5) * 255);
              n.y = static_cast<uint8_t>((normal.y * 0.5 + 0.5) * 255);
              n.z = static_cast<uint8_t>((normal.z * 0.5 + 0.5) * 255);
              r = glm::clamp(2*abs(int(n.x) - 128), 0, 255);
              g = glm::clamp(2*abs(int(n.y) - 128), 0, 255);
              b = glm::clamp(2*abs(int(n.z) - 128), 0, 255);
              #pragma omp critical 
              {
                out_fobj 
                        << x_coord << " "
                        << space.minPoint.y + y*node_size + node_size/2.f << " "
                        << space.minPoint.z + z*node_size + node_size/2.f << " "
                        << static_cast<int>(r) << " "
                        << static_cast<int>(g) << " "
                        << static_cast<int>(b)
                        <<'\n';
                num_points++;
              }
            }
            if (calculate_b_w) {
              #pragma omp critical
              {
                voxels.push_back(v);
              }
            }
          }
        }

        if (calculate_b_w) {
          // Extract this to function
          std::vector<double> intersect_xs;
          for (int tri_idx : quad_node->members) {
            glm::vec3 intersection_point;
            glm::vec3 v1, v2, v3;
            Triangle t = triangles[tri_idx];
            v1 = vertices[t.getV1()];
            v2 = vertices[t.getV2()];
            v3 = vertices[t.getV3()];
            // Ray in the middle of the node pointing in the X direction
            bool intersects = Geo::rayIntersectsTriangle(glm::vec3(space.minPoint.x - 0.5f, coords.x + node_size/2.f, coords.y + node_size/2.f), glm::vec3(1,0,0), v1, v2, v3, intersection_point);
            // Keep all the intersections of this row
            if (intersects) {          
              intersect_xs.push_back(intersection_point.x);
            }
          }
          // Number of intersection so far
          uint x_idx = 0;
          for (uint x = 0; x < size_; x++) {
            if (voxels[x].color == VoxelColor::GRAY) continue;
            double x_coord = space.minPoint.x + x * node_size;
            while ( x_idx < intersect_xs.size() && x_coord > intersect_xs[x_idx]) x_idx++;

            // If number of intersections so far is odd, we are inside the model
            if (x_idx % 2 == 1) voxels[x].color = VoxelColor::BLACK;
          }
          #pragma omp critical
          {
            std::cout << "Row Coloring: ";
            for (Voxel v : voxels) {
              if (v.color == VoxelColor::BLACK) std::cout << "B";
              if (v.color == VoxelColor::WHITE) std::cout << "W";
              if (v.color == VoxelColor::GRAY) std::cout << "G";
            }
            std::cout << endl;
          }
        }
      }
    }
  }

  out_fobj.seekp(file_verts_line, std::ios::beg);
  out_fobj << string_format("element vertex %20d\r\n", num_points);
  out_fobj.close();
}

// void Grid::writePLY(std::string filename) {
//   std::ofstream out_fobj(filename);

//   out_fobj<< "ply\r\n"
//           << "format ascii 1.0\r\n";
//   auto file_verts_line = out_fobj.tellp();

//   out_fobj<< string_format("element vertex %20d\r\n", 0)
//           << "property float x\r\n"
//           << "property float y\r\n"
//           << "property float z\r\n"
//           << "property uchar red\r\n"
//           << "property uchar green\r\n"
//           << "property uchar blue\r\n"
//           << "end_header\r\n";

//   uint num_points = 0;

//   unsigned char r, g, b;
//   r = g = b = 255u;
  
//   for (uint y = 0; y < size_; y++) {
//     for (uint z = 0; z < size_; z++) {
//       for (uint x = 0; x < size_; x++) {
//         Voxel* v = elements[y*size_*size_+z*size_+x];
//         if (v->color == VoxelColor::GRAY) {

//           glm::vec3 n;
//           n.x = static_cast<uint8_t>((v->normal.x * 0.5 + 0.5) * 255);
//           n.y = static_cast<uint8_t>((v->normal.y * 0.5 + 0.5) * 255);
//           n.z = static_cast<uint8_t>((v->normal.z * 0.5 + 0.5) * 255);
//           r = glm::clamp(2*abs(int(n.x) - 128), 0, 255);
//           g = glm::clamp(2*abs(int(n.y) - 128), 0, 255);
//           b = glm::clamp(2*abs(int(n.z) - 128), 0, 255);

//           out_fobj 
//                   << v->x << " "
//                   << space.minPoint.y + y*node_size + node_size/2.f << " "
//                   << space.minPoint.z + z*node_size + node_size/2.f << " "
//                   << static_cast<int>(r) << " "
//                   << static_cast<int>(g) << " "
//                   << static_cast<int>(b)
//                   <<'\n';
//           num_points++;
//         }
//       }
//     }
//   }

//   out_fobj.seekp(file_verts_line, std::ios::beg);
//   out_fobj << string_format("element vertex %20d\r\n", num_points);
//   out_fobj.close();
// }