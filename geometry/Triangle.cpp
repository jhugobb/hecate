#include <glm/common.hpp>
#include <vector>

#include "Triangle.h"
#include "../TriangleMesh.h"

Triangle::Triangle(unsigned int v1_, unsigned int v2_, unsigned int v3_) {
  v1 = v1_;
  v2 = v2_;
  v3 = v3_;
}

unsigned int Triangle::getV1() {
  return v1;
}

unsigned int Triangle::getV2() {
  return v2;
}

unsigned int Triangle::getV3() {
  return v3;
}

Triangle::~Triangle() {}

void Triangle::saveMinX(TriangleMesh* m) {
  std::vector<glm::vec3> verts = m->getVertices();

  glm::vec3 v1_v = verts[v1];
  glm::vec3 v2_v = verts[v2];
  glm::vec3 v3_v = verts[v3];

  min_x = std::min(v1_v.x, std::min(v2_v.x, v3_v.x));
}

void Triangle::calculateTriBbox(std::vector<glm::vec3> &vertices) {
  tri_bbox.addPoint(vertices[v1]);
  tri_bbox.addPoint(vertices[v2]);
  tri_bbox.addPoint(vertices[v3]);
}