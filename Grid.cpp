#include <algorithm>

#include <memory>
#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>
#include <random>
#include <bitset>

#include "Grid.h"
#include "geometry/Geometry.h"

// PNG library
#include "lodepng/lodepng.h"


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

void Grid::testVoxelGray_Box(int z, 
                             Voxel* voxels, 
                             std::vector<int> &candidates, 
                             glm::vec2 &row_coords) {
  std::vector<Triangle> triangles = mesh_->getTriangles();
  std::vector<glm::vec3> verts = mesh_->getVertices();
  // Box approach: for each triangle in the row, check which voxels it intersects
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
    for (int x = initial_x_grid; x < final_x_grid + 1 && x < (int) size_; x++) {
      
      if (voxels[z * size_ + x].color == VoxelColor::GRAY) continue;

      min_box = glm::vec3(voxels[z * size_ + x].x, row_coords.x, row_coords.y);
      max_box = glm::vec3(voxels[z * size_ + x].x + node_size, row_coords.x + node_size, row_coords.y + node_size);

      if (Geo::testBoxTriangle(mesh_, triangles[tri_idx], min_box, max_box)) {
        voxels[z * size_ + x].color = VoxelColor::GRAY;
        voxels[z * size_ + x].normal = glm::normalize(glm::cross(v3-v2, v3-v1));
      }
    }
  }
}

void Grid::writePLY(int x, int y, int z, Voxel &voxel, std::ofstream &out_fobj) const {
  glm::vec3 n;
  double x_coord = space.minPoint.x + x * node_size;
  unsigned char r, g, b;
  n.x = static_cast<uint8_t>((voxel.normal.x * 0.5 + 0.5) * 255);
  n.y = static_cast<uint8_t>((voxel.normal.y * 0.5 + 0.5) * 255);
  n.z = static_cast<uint8_t>((voxel.normal.z * 0.5 + 0.5) * 255);
  r = glm::clamp(2*abs(int(n.x) - 128), 0, 255);
  g = glm::clamp(2*abs(int(n.y) - 128), 0, 255);
  b = glm::clamp(2*abs(int(n.z) - 128), 0, 255);

  out_fobj 
          << x_coord << " "
          << space.minPoint.y + y*node_size + node_size/2.f << " "
          << space.minPoint.z + z*node_size + node_size/2.f << " "
          << static_cast<int>(r) << " "
          << static_cast<int>(g) << " "
          << static_cast<int>(b)
          <<'\n';
}

void Grid::calculateBlackWhite(int z, 
                               glm::vec2 coords, 
                               node* quad_node, 
                               Voxel* voxels, 
                               double threshold) {
  std::vector<Triangle> triangles = mesh_->getTriangles();
  std::vector<glm::vec3> vertices = mesh_->getVertices();
  // One row of intersections for each ray
  std::vector<double> intersect_xs;
  glm::vec3 origin;
  // random double engine
  double lower_bound = 0+0.00001;
  double upper_bound = node_size-0.0000001;
  std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
  std::default_random_engine re;

  origin = glm::vec3(space.minPoint.x - 0.5f, coords.x + unif(re), coords.y + unif(re));

  for (int tri_idx : quad_node->members) {
    glm::vec3 v1, v2, v3;
    Triangle t = triangles[tri_idx];
    v1 = vertices[t.getV1()];
    v2 = vertices[t.getV2()];
    v3 = vertices[t.getV3()];
    glm::vec3 rayDirection = glm::vec3(1,0,0);

    // Keep changing rays until all are valid
    bool intersects = false;
    glm::vec3 intersection_point;
    Geo::IntersectionResult res = Geo::rayIntersectsTriangle(origin, rayDirection, v1, v2, v3, threshold, intersection_point);
    while (res == Geo::IntersectionResult::INVALID) {
      origin = glm::vec3(space.minPoint.x - node_size, coords.x + unif(re), coords.y + unif(re));
      res = Geo::rayIntersectsTriangle(origin, rayDirection, v1, v2, v3, threshold, intersection_point);
    }
    intersects = res == Geo::IntersectionResult::INTERSECTS;
    if (intersects) intersect_xs.push_back(intersection_point.x);
  }

  // sort them in ascending X order
  std::sort(intersect_xs.begin(), intersect_xs.end());
  // Number of intersection so far
  uint x_idx = 0;
  for (uint x = 0; x < size_; x++) {
    if (voxels[z * size_ + x].color == VoxelColor::GRAY) continue;
    double x_coord = space.minPoint.x + x * node_size;
    while (x_idx < intersect_xs.size() && x_coord > intersect_xs[x_idx]) x_idx++;
    // If number of intersections so far is odd , we are inside the model
    if (x_idx % 2 == 1){
        voxels[z * size_ + x].color = VoxelColor::BLACK;
    }
  }
}

void Grid::colorGrid(TriangleMesh* mesh, TwoDGrid* qt, ColoringConfiguration config, std::string filename) {
  std::vector<Triangle> triangles = mesh->getTriangles();
  std::vector<glm::vec3> vertices = mesh->getVertices();
  mesh_ = mesh;

  // For PLY writing
  std::ofstream out_fobj(filename);
  std::streampos file_verts_line;
  uint num_points = 0;

  // For Binary file writing
  std::ofstream bin_file("test.hec", std::ios::binary | std::ios::out);

  if (config.writePLY) {
    out_fobj<< "ply\r\n"
            << "format ascii 1.0\r\n";

    file_verts_line = out_fobj.tellp();

    out_fobj<< string_format("element vertex %20d\r\n", 0)
            << "property float x\r\n"
            << "property float y\r\n"
            << "property float z\r\n"
            << "property uchar red\r\n"
            << "property uchar green\r\n"
            << "property uchar blue\r\n"
            << "end_header\r\n";
  }

  for (unsigned int y = 0; y < size_; y++) {
    Voxel voxels[size_*size_];
    #pragma omp parallel for
    for (unsigned int z = 0; z < size_; z++) {
      glm::vec2 coords;

      // Coords of voxel in the 2D grid dimensions
      coords.x = space.minPoint.y + y * node_size;
      coords.y = space.minPoint.z + z * node_size;

      // 2Dgrid node that has the triangles for this y z 
      node* quad_node = qt->query(glm::vec2(y,z));
      
      std::list<BinTreeNode*>::iterator list_it;
      std::list<BinTreeNode*> nods;

      // If the node has triangles
      if (quad_node->members.size() != 0) {
        if (!config.useNaive && !config.useBox) {
          nods = quad_node->bin_tree;
          list_it = nods.begin();
        }
        if (config.useBox) {
          for (uint x = 0; x < size_; x++) {
            Voxel v = Voxel();
            v.x = space.minPoint.x + x * node_size;
            voxels[z*size_ + x] = v;
          }
          testVoxelGray_Box(z, voxels, quad_node->members, coords);
          for (uint x = 0; x < size_; x++) {
            if (voxels[z*size_ + x].color != VoxelColor::GRAY) continue;
            if (config.writePLY) {
              #pragma omp critical 
              {
                writePLY(x, y, z, voxels[z*size_ + x], out_fobj);
                num_points++;
              }
            }
          }
        } else {
          for (uint x = 0; x < size_; x++) {
            Voxel v;
            // coordinate in world space
            double x_coord = space.minPoint.x + x * node_size;
            bool intersects_node = false;
            int representative = -1;
            if (config.useNaive) {
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
              if (config.writePLY) {
                glm::vec3 normal;
                Triangle t = triangles[representative];
                glm::vec3 v1 = vertices[t.getV1()];
                glm::vec3 v2 = vertices[t.getV2()];
                glm::vec3 v3 = vertices[t.getV3()];

                normal = glm::normalize(glm::cross(v3-v2, v3-v1));
                v.normal = normal;
                #pragma omp critical 
                {
                  writePLY(x, y, z, v, out_fobj);
                  num_points++;
                }
              }
            }

            if (config.calculate_black_white) {
              voxels[z*size_ + x] = v;
            }
          }
        }

        if (config.calculate_black_white) {
          calculateBlackWhite(z, coords, quad_node, voxels, config.threshold_raycasting);
        }
      }
    }
    if (config.writePNG) {
      #pragma omp critical 
      {
        saveSliceAsPNG(voxels, y);
      }
    }

    if (config.writeHEC) {
      #pragma omp critical
      {
        // Write Hecate (binary file)
        saveSliceAsHEC(voxels, bin_file);
      }
    }
  }

  if (config.writePLY) {
    out_fobj.seekp(file_verts_line, std::ios::beg);
    out_fobj << string_format("element vertex %20d\r\n", num_points);
  }
  out_fobj.close();
  bin_file.close();

}

void Grid::saveSliceAsPNG(Voxel* voxels, uint y) {

  std::vector<unsigned char> image(size_*size_*4);
  for (uint x = 0; x < size_; x++) {
    for (uint z = 0; z < size_; z++) {
      // int index = z * size_+ x;
      unsigned int value = 0;
      if (voxels[x*size_+z].color == VoxelColor::GRAY) {
        value = 128;
      } else if (voxels[x*size_+z].color == VoxelColor::WHITE) {
        value = 255;
      }
      image[4 * size_ * x + 4 * z + 0] = value;
      image[4 * size_ * x + 4 * z + 1] = value;
      image[4 * size_ * x + 4 * z + 2] = value;
      image[4 * size_ * x + 4 * z + 3] = 255;
    }
  }

  // we're going to encode with a state rather than a convenient function, because enforcing a color type requires setting options
  // lodepng::State state;
  // // input color type
  // state.info_raw.colortype = LCT_GREY;
  // state.info_raw.bitdepth = 1;
  // // output color type
  // state.info_png.color.colortype = LCT_GREY;
  // state.info_png.color.bitdepth = 1;
  // state.encoder.auto_convert = 0; 
  // std::vector<unsigned char> buffer;
  // unsigned error = lodepng::encode(buffer, &image[0], size_, size_ , state);

  std::string filename = "slices/slice_";
  filename.append(std::to_string(y));
  filename.append(".png");
  unsigned error = lodepng::encode(filename, &image[0], size_, size_);


  //if there's an error, display it
  if(error) std::cout << "encoder error " << error << ": "<< lodepng_error_text(error) << std::endl;
  // else lodepng::save_file(buffer, filename);
}

void Grid::saveSliceAsHEC(Voxel* voxels, std::ofstream &bin_file) {

  char* slice_memblock;

  if (bin_file.is_open()) {
    slice_memblock = new char[(size_*size_)/4];
    uint slice_idx = 0;
    for (uint i = 0; i < size_*size_; i+=4) {
      std::bitset<8> bits(0);

      if (voxels[i].color == VoxelColor::BLACK) {
        bits.set(1);
      } else if (voxels[i].color == VoxelColor::GRAY) {
        bits.set(0);
      }

      if (voxels[i+1].color == VoxelColor::BLACK) {
        bits.set(3);
      } else if (voxels[i+1].color == VoxelColor::GRAY) {
        bits.set(2);
      }

      if (voxels[i+2].color == VoxelColor::BLACK) {
        bits.set(5);
      } else if (voxels[i+2].color == VoxelColor::GRAY) {
        bits.set(4);
      }

      if (voxels[i+3].color == VoxelColor::BLACK) {
        bits.set(7);
      } else if (voxels[i+3].color == VoxelColor::GRAY) {
        bits.set(6);
      }

      slice_memblock[slice_idx++] =  static_cast<unsigned char>(bits.to_ulong());
    }

    bin_file.write(slice_memblock, size_*size_/4);

    delete[] slice_memblock;
  } else {
    std::cout << "Unable to open file." << std::endl;
    assert(false);
  }

}