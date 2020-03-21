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
  uint key;

  // Create Voxel Grid
  elements = std::vector<Voxel*>(size*size*size, NULL);
  for (uint y = 0; y < size; y++) {
    for (uint z = 0; z < size; z++) {
      for (uint x = 0; x < size; x++) {
        key = y*size*size + z*size + x;
        Voxel* voxel = new Voxel();
        voxel->x = space.minPoint.x + x * node_size;
        elements[key] = voxel;
      }

    }
  }
  cout << "Coords of min voxel: " << space.minPoint.x << " " << space.minPoint.y << " " << space.minPoint.z << endl
       <<"                      " << space.minPoint.x + node_size << " " << space.minPoint.y + node_size << " " << space.minPoint.z + node_size << endl;
  cout << "Coords of max voxel: " << space.minPoint.x + (size-2)*node_size << " " << space.minPoint.y + (size-2)*node_size << " " << space.minPoint.z + (size-2)*node_size << endl
       <<"                      " << space.minPoint.x + (size-1)*node_size << " " << space.minPoint.y + (size-1)*node_size << " " << space.minPoint.z + (size-1)*node_size << endl;  
}

void Grid::colorGrid(TriangleMesh* mesh, TwoDGrid* qt) {
  std::vector<Triangle> triangles = mesh->getTriangles();
  std::vector<glm::vec3> vertices = mesh->getVertices();
  #pragma omp parallel for
  for (unsigned int y = 0; y < size_; y++) {
    #pragma omp parallel for
    for (unsigned int z = 0; z < size_; z++) {
      glm::vec2 coords;
      // Coords of voxel in the 2D grid dimensions
      coords.x = space.minPoint.y + y * node_size;
      coords.y = space.minPoint.z + z * node_size;
      node* quad_node = qt->query(glm::vec2(y,z));
      
      std::list<BinTreeNode*>::iterator list_it;
      std::list<BinTreeNode*> nods;
      if (quad_node->members.size() != 0) {
        nods = quad_node->bin_tree;
        list_it = nods.begin();
      
        for (uint x = 0; x < size_; x++) {
          Voxel* v = elements[y*size_*size_+z*size_+x];
          bool intersects_node = false;
          if (quad_node->members.size() != 0) {
            while ((*list_it)->max_point.x <= v->x){
              ++list_it;
            }
            intersects_node = (*list_it)->is_gray;
          }
          if (intersects_node){
            v->color = VoxelColor::GRAY;
            Triangle t = triangles[(*list_it)->representative];
            glm::vec3 v1 = vertices[t.getV1()];
            glm::vec3 v2 = vertices[t.getV2()];
            glm::vec3 v3 = vertices[t.getV3()];

            v->normal = glm::normalize(glm::cross(v3-v2, v3-v1));
          }
        }
      }
      // if (quad_node->members.size() != 0) {
      //   for (BinTreeNode* n : nods) {
      //     delete n;
      //   }
      //   nods.clear();
      // }
      // for (Triangle t : triangles) {
      //   // Triangle t = triangles[tri_idx];
      //   glm::vec3 intersection_point;
      //   glm::vec3 v1 = vertices[t.getV1()];
      //   glm::vec3 v2 = vertices[t.getV2()];
      //   glm::vec3 v3 = vertices[t.getV3()];

      //   // Ray in the middle of the node pointing upward
      //   bool intersects = Geo::rayIntersectsTriangle(glm::vec3(coords.x + node_size/2.f, coords.y + node_size/2.f, space.minPoint.z - 0.5f), glm::vec3(0,0,1), v1, v2, v3, intersection_point);

      //   if (intersects) {          
      //     // std::cout<<"int"<<std::endl;
      //     intersect_zs.push_back(intersection_point.z);
      //   }
      // }
      // std::sort(intersect_zs.begin(), intersect_zs.end());
      // if (intersect_zs.size() != 0) {
      //   unsigned int z_idx = 0;
      //   std::set<Voxel*>::iterator it;
      //   for (it = elements[std::make_pair(x,y)].begin(); it != elements[std::make_pair(x,y)].end(); ++it) {
      //     Voxel* v = *it;
      //     while (!(z_idx > intersect_zs.size()) && !(v->z <= intersect_zs[z_idx])) {
      //       z_idx++;
      //     }
      //     if (z_idx != 0) { 
      //       bool intersects_node = false;
      //       for (Triangle t : triangles) {
      //         // Triangle t = triangles[tri_idx];
      //         intersects_node = Geo::testBoxTriangle(mesh, t, glm::vec3(coords,v->z), glm::vec3(coords.x+node_size, coords.y+node_size, v->z+node_size));
      //         if (intersects_node) {
      //           glm::vec3 v1 = vertices[t.getV1()];
      //           glm::vec3 v2 = vertices[t.getV2()];
      //           glm::vec3 v3 = vertices[t.getV3()];

      //           v->normal = glm::normalize(glm::cross(v3-v2, v3-v1));
      //           break;
      //         }
      //       }
      //       if (intersects_node){
      //         v->color = VoxelColor::GRAY;
      //       } else if (z_idx % 2 == 1) {
      //         v->color = VoxelColor::BLACK;
      //       }
          // }
        // }
      // }
    }
  }
}

void Grid::writePLY(std::string filename) {
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

  unsigned char r, g, b;
  r = g = b = 255u;
  
  for (uint y = 0; y < size_; y++) {
    for (uint z = 0; z < size_; z++) {
      for (uint x = 0; x < size_; x++) {
        Voxel* v = elements[y*size_*size_+z*size_+x];
        if (v->color == VoxelColor::GRAY) {

          glm::vec3 n;
          n.x = static_cast<uint8_t>((v->normal.x * 0.5 + 0.5) * 255);
          n.y = static_cast<uint8_t>((v->normal.y * 0.5 + 0.5) * 255);
          n.z = static_cast<uint8_t>((v->normal.z * 0.5 + 0.5) * 255);
          r = glm::clamp(2*abs(int(n.x) - 128), 0, 255);
          g = glm::clamp(2*abs(int(n.y) - 128), 0, 255);
          b = glm::clamp(2*abs(int(n.z) - 128), 0, 255);

          out_fobj 
                  << v->x << " "
                  << space.minPoint.y + y*node_size + node_size/2.f << " "
                  << space.minPoint.z + z*node_size + node_size/2.f << " "
                  << static_cast<int>(r) << " "
                  << static_cast<int>(g) << " "
                  << static_cast<int>(b)
                  <<'\n';
          num_points++;
        }
      }
    }
  }

  out_fobj.seekp(file_verts_line, std::ios::beg);
  out_fobj << string_format("element vertex %20d\r\n", num_points);
  out_fobj.close();
}