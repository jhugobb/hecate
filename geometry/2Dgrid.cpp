#include "2Dgrid.h"

#include <iostream>
#include <map>

TwoDGrid::TwoDGrid(TriangleMesh* m_, int size, Geo::BBox space_) : triangles(m_->getTriangles()){
  glm::vec3 center_box = (space_.maxPoint + space_.minPoint)/2.0f;
  glm::vec3 size_box = space_.maxPoint - space_.minPoint;
  // Make the cells cubic
  float largest = max(size_box.x, max(size_box.y, size_box.z));
  space.minPoint = center_box - glm::vec3(largest, largest, largest)/2.0f;
  space.maxPoint = center_box + glm::vec3(largest, largest, largest)/2.0f;
  node_size = largest/size;
  num_nodes = size;
  m = m_;

  // Allocate space for the 2Dgrid
  grid = new node*[num_nodes*num_nodes];

  for (int y = 0; y < num_nodes; y++) {
    for (int z = 0; z < num_nodes; z++) {

      glm::vec2 min_box = glm::vec2(space.minPoint.y, space.minPoint.z) + glm::vec2(y*node_size,z*node_size);
      glm::vec2 max_box = glm::vec2(space.minPoint.y, space.minPoint.z) + glm::vec2((y+1)*node_size,(z+1)*node_size);
      grid[y * num_nodes + z] = new node();
      grid[y * num_nodes + z]->min_point = min_box;
      grid[y * num_nodes + z]->max_point = max_box;
      grid[y * num_nodes + z]->members = std::vector<int>();
    }
  }

}

TwoDGrid::~TwoDGrid() {
  for (int i = 0; i < num_nodes*num_nodes; i++) {
    delete grid[i];
  }
  delete[] grid;
}

void TwoDGrid::insert(int t) {

  Triangle* tri = triangles[t];

  // Save Bbox of triangle
  glm::vec3 min_box_tri = tri->tri_bbox.minPoint;
  glm::vec3 max_box_tri = tri->tri_bbox.maxPoint;

  // translate to origin
  min_box_tri -= space.minPoint;
  max_box_tri -= space.minPoint;

  int initial_y_grid = std::max(min_box_tri.y / node_size - 1, 0.0);
  int final_y_grid = max_box_tri.y / node_size;

  int initial_z_grid = std::max(min_box_tri.z / node_size - 1, 0.0);
  int final_z_grid = max_box_tri.z / node_size;

  uint key;

  for (int y = initial_y_grid; y < final_y_grid + 1 && y < num_nodes; y++) {
    for (int z = initial_z_grid; z < final_z_grid + 1 && z < num_nodes; z++) {
      key = y * num_nodes + z;
      glm::vec2 &min_box = grid[key]->min_point;
      glm::vec2 &max_box = grid[key]->max_point;

      // If node intersects the projection in X of the triangle
      if (Geo::testQuadTriangle(m, tri, min_box, max_box)) {
        #pragma omp critical 
        {
          grid[key]->members.push_back(t);
        }
      }
    }
  }

}

void TwoDGrid::saveMinXs() {

  #pragma omp parallel for
  for (uint i = 0; i < triangles.size(); i++) {
    triangles[i]->saveMinX(m);
  }
}

void TwoDGrid::buildBinTree(int y, int z) {
  int key = y*num_nodes + z;
  // If node has triangles, calculate its bin tree
  if (grid[key]->members.size() != 0)
    grid[key]->build_bin_tree(m, space, node_size);
}

node* TwoDGrid::query(glm::vec2 coords) {
  int y = floor(coords.x);
  int z = floor(coords.y);

  uint key = y * num_nodes + z;

  return grid[key];
}

void node::build_bin_tree(TriangleMesh* mesh, Geo::BBox space, double min_node_size) {
  std::list<BinTreeNode*> result;
  std::list<BinTreeNode*>::iterator it;

  std::vector<Triangle*> &triangles = mesh->getTriangles();
  std::map<Triangle, int> reference;

  BinTreeNode* btn = new BinTreeNode();
  btn->min_point = glm::vec3(space.minPoint.x, min_point.x, min_point.y);
  btn->max_point = glm::vec3(space.maxPoint.x+min_node_size, max_point.x, max_point.y);

  for (int t : members) {
    btn->triangles.insert(*triangles[t]);
    reference.emplace(*triangles[t], t);
  }

  result.push_back(btn);



  it = result.begin();
  bool needs_subdivision = true;
  BinTreeNode* current;
  std::multiset<Triangle>::iterator it_tri, current_it_tri;
  
  while(it!=result.end()) {
    needs_subdivision = false;
    current = *it;
    it_tri = current->triangles.begin();
    while(it_tri != current->triangles.end()) {
      current_it_tri = it_tri++;
      Triangle t = *current_it_tri;
      // If triangle intersects node
      if (t.min_x < current->max_point.x && Geo::testBoxTriangle(mesh, &t, current->min_point, current->max_point)) {
        // We need to subdivide only if half the size of the cell is larger than the minimum size
        if (abs(current->max_point.x - current->min_point.x)/2.0 < min_node_size) {
          current->is_gray = true;
          current->representative = reference[t];
          break;
        } else {
          needs_subdivision = true;
          BinTreeNode* nbtn = new BinTreeNode();
          double center = (current->min_point.x + current->max_point.x)/2.0;
          nbtn->min_point = glm::vec3(current->min_point.x, min_point.x, min_point.y);
          nbtn->max_point = glm::vec3(center, max_point.x, max_point.y);
          nbtn->triangles = std::multiset<Triangle>(current->triangles);
          current->min_point.x = center;
          result.insert(it, nbtn);

          // Make sure to go back and process the new cell
          --it;

          break;
        }
      } else {
        current->triangles.erase(current_it_tri);
      }
    }
    // Only advance if we didnt create new cells to be processed
    if (!needs_subdivision) ++it;
  }
  bin_tree = result;
}