#include "Quadtree.h"
#include "Geometry.h"

Quadtree::Quadtree() {}

Quadtree::Quadtree(TriangleMesh* mesh) {
  m = mesh;
}

Quadtree::~Quadtree() {}

node* Quadtree::insert_recursive(int t, node* cell, TriangleMesh* mesh, int depth) {
  if (cell != NULL) {
    if (cell->is_leaf) {
      if (depth >= 8 || cell->members.size() < 17) {
        cell->members.push_back(t);
      } else {
        cell->children[0] = new node();
        cell->children[1] = new node();
        cell->children[2] = new node();
        cell->children[3] = new node();

        cell->is_leaf = false;
        glm::vec2 center = (cell->min_point + cell->max_point) / 2.0f; 
        
        //bottom left quadrant
        cell->children[0]->min_point = cell->min_point;
        cell->children[0]->max_point = center;
        
        //top left quadrant
        cell->children[1]->min_point = glm::vec2(cell->min_point.x, center.y); 
        cell->children[1]->max_point = glm::vec2(center.x, cell->max_point.y);

        //top right quadrant
        cell->children[2]->min_point = center; 
        cell->children[2]->max_point = cell->max_point;

        //bottom right quadrant
        cell->children[3]->min_point = glm::vec2(center.x, cell->min_point.y); 
        cell->children[3]->max_point = glm::vec2(cell->max_point.x, center.y);

        cell->members.push_back(t);
        for (int tri : cell->members) {

          Triangle triangle = mesh->getTriangles()[tri];
          // can a triangle be in different cells?
          if (Geo::testQuadTriangle(mesh, triangle, glm::vec3(cell->children[0]->min_point,0), glm::vec3(cell->children[0]->max_point,0))){
            cell->children[0] = insert_recursive(tri, cell->children[0], mesh, depth+1);
          } 
          if (Geo::testQuadTriangle(mesh, triangle, glm::vec3(cell->children[1]->min_point,0), glm::vec3(cell->children[1]->max_point,0))) {
            cell->children[1] = insert_recursive(tri, cell->children[1], mesh, depth+1);
          } 
          if (Geo::testQuadTriangle(mesh, triangle, glm::vec3(cell->children[2]->min_point,0), glm::vec3(cell->children[2]->max_point,0))) {
            cell->children[2] = insert_recursive(tri, cell->children[2], mesh, depth+1);
          }
          if (Geo::testQuadTriangle(mesh, triangle, glm::vec3(cell->children[3]->min_point,0), glm::vec3(cell->children[3]->max_point,0))) {
            cell->children[3] = insert_recursive(tri, cell->children[3], mesh, depth+1);
          }

        }

        cell->members.clear();
      }
    } else {
      Triangle triangle = mesh->getTriangles()[t];
      if (Geo::testQuadTriangle(mesh, triangle, glm::vec3(cell->children[0]->min_point,0), glm::vec3(cell->children[0]->max_point,0))){
        cell->children[0] = insert_recursive(t, cell->children[0], mesh, depth+1);
      } 
      if (Geo::testQuadTriangle(mesh, triangle, glm::vec3(cell->children[1]->min_point,0), glm::vec3(cell->children[1]->max_point,0))) {
        cell->children[1] = insert_recursive(t, cell->children[1], mesh, depth+1);
      } 
      if (Geo::testQuadTriangle(mesh, triangle, glm::vec3(cell->children[2]->min_point,0), glm::vec3(cell->children[2]->max_point,0))) {
        cell->children[2] = insert_recursive(t, cell->children[2], mesh, depth+1);
      }
      if (Geo::testQuadTriangle(mesh, triangle, glm::vec3(cell->children[3]->min_point,0), glm::vec3(cell->children[3]->max_point,0))) {
        cell->children[3] = insert_recursive(t, cell->children[3], mesh, depth+1);
      }
    }
  } else {
    cell = new node();
    glm::vec3 box[2];
    mesh->getBBox(box);
    cell->min_point = box[0];
    cell->max_point = box[1];
    cell->members = std::vector<int>();
    cell->members.push_back(t);
    cell->is_leaf = true;
  }
  return cell;
}

void Quadtree::insert(int t) {
  root = insert_recursive(t, root, m, 0);
}

node* Quadtree::query_recursive(node* cell, glm::vec2 coords) {
  if (cell != NULL) {
    if (cell->is_leaf)  return cell;
    else {
      if (Geo::isPointInQuad(cell->children[0]->min_point, cell->children[0]->max_point, coords)) {
        return query_recursive(cell->children[0], coords);
      } else if (Geo::isPointInQuad(cell->children[1]->min_point, cell->children[1]->max_point, coords)) {
        return query_recursive(cell->children[1], coords);
      } else if (Geo::isPointInQuad(cell->children[2]->min_point, cell->children[2]->max_point, coords)) {
        return query_recursive(cell->children[2], coords);
      } else query_recursive(cell->children[3], coords);
    }
  } else return cell;
}

node* Quadtree::query(glm::vec2 coords) {
  return query_recursive(root, coords);
}