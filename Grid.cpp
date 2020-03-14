#include <algorithm>

#include <memory>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>

#include "Grid.h"
#include "geometry/Geometry.h"


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
  glm::vec3 center_box = (space.maxPoint + space.minPoint)/2.0f;
  glm::vec3 size_box = space.maxPoint - space.minPoint;
  float largest = max(size_box.x, max(size_box.y, size_box.z));
  space.minPoint = center_box - glm::vec3(largest, largest, largest)/2.0f;
  space.maxPoint = center_box + glm::vec3(largest, largest, largest)/2.0f;
  node_size = largest/size;
  for (unsigned int y = 0; y < size; y++) {
    for (unsigned int z = 0; z < size; z++) {
      std::set<Voxel*> s;
      for (unsigned int x = 0; x < size; x++) {
        Voxel* voxel = new Voxel();
        voxel->x = space.minPoint.x + x * node_size;
        s.insert(voxel);
      }
      elements.emplace(std::make_pair(y,z), s);
    }
  }
  cout << "Coords of min voxel: " << space.minPoint.x << " " << space.minPoint.y << " " << space.minPoint.z << endl
       <<"                      " << space.minPoint.x + node_size << " " << space.minPoint.y + node_size << " " << space.minPoint.z + node_size << endl;
  cout << "Coords of max voxel: " << space.minPoint.x + (size-2)*node_size << " " << space.minPoint.y + (size-2)*node_size << " " << space.minPoint.z + (size-2)*node_size << endl
       <<"                      " << space.minPoint.x + (size-1)*node_size << " " << space.minPoint.y + (size-1)*node_size << " " << space.minPoint.z + (size-1)*node_size << endl;  
}

void Grid::colorGrid(TriangleMesh* mesh, TwoDGrid qt) {
  std::vector<Triangle> triangles = mesh->getTriangles();
  std::vector<glm::vec3> vertices = mesh->getVertices();
  for (unsigned int y = 0; y < size_; y++) {
    for (unsigned int z = 0; z < size_; z++) {
      glm::vec2 coords;
      coords.x = space.minPoint.y + y * node_size;
      coords.y = space.minPoint.z + z * node_size;
      node* quad_node = qt.query(glm::vec2(y,z));
      // std::vector<double> intersect_zs;
      
      std::set<Voxel*>::iterator it;
      for (it = elements[std::make_pair(y,z)].begin(); it != elements[std::make_pair(y,z)].end(); ++it) {
        Voxel* v = *it;
        bool intersects_node = false;
        if (quad_node != NULL) {
          for (unsigned int tri_idx : quad_node->members) {
            Triangle t = triangles[tri_idx];
            intersects_node = Geo::testBoxTriangle(mesh, t, glm::vec3(v->x,coords.x, coords.y), glm::vec3(v->x+node_size, coords.x+node_size, coords.y+node_size));
            if (intersects_node) {
              glm::vec3 v1 = vertices[t.getV1()];
              glm::vec3 v2 = vertices[t.getV2()];
              glm::vec3 v3 = vertices[t.getV3()];

              v->normal = glm::normalize(glm::cross(v3-v2, v3-v1));
              break;
            }
          }
        }
        if (intersects_node){
          v->color = VoxelColor::GRAY;
        }
      }
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

  unsigned int num_points = 0;
  const unsigned int LINE_SIZE = (20u-1u)*3u + 17u; // (CHARS_PER_COORD-1)*3 position, 3*3 color, 5 spaces, 2 endl, 1 null

  unsigned char r, g, b;
  r = g = b = 255u;
  
  for (unsigned int y = 0; y < size_; y++) {
    for (unsigned int z = 0; z < size_; z++) {
      for (Voxel* v : elements[std::make_pair(y,z)]) {
        if (v->color == VoxelColor::GRAY) {

          glm::vec3 n;
          n.x = static_cast<uint8_t>((v->normal.x * 0.5 + 0.5) * 255);
          n.y = static_cast<uint8_t>((v->normal.y * 0.5 + 0.5) * 255);
          n.z = static_cast<uint8_t>((v->normal.z * 0.5 + 0.5) * 255);
          r = glm::clamp(2*abs(int(n.x) - 128), 0, 255);
          g = glm::clamp(2*abs(int(n.y) - 128), 0, 255);
          b = glm::clamp(2*abs(int(n.z) - 128), 0, 255);

          char line[LINE_SIZE];
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