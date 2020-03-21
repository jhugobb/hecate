#ifndef GEO_H
#define GEO_H

#include <glm/glm.hpp>

#include "../TriangleMesh.h"
#include "Triangle.h"

namespace Geo {

const float EPSILON = 0.0000001;

class BBox {
  public:
    glm::vec3 minPoint;
    glm::vec3 maxPoint;


    BBox() {
      init();
    }

    void init() {
      minPoint = glm::vec3(std::numeric_limits<double>::infinity());
      maxPoint = glm::vec3(-std::numeric_limits<double>::infinity());
    }

    void addPoint(const glm::vec3 &point) {
      minPoint = glm::min(minPoint, point);
      maxPoint = glm::max(maxPoint, point);
    }

    void addBox(const BBox &bbox) {
      addPoint(glm::vec3(bbox.minPoint.x, bbox.minPoint.y, bbox.minPoint.z));
      addPoint(glm::vec3(bbox.minPoint.x, bbox.minPoint.y, bbox.maxPoint.z));
      addPoint(glm::vec3(bbox.minPoint.x, bbox.maxPoint.y, bbox.minPoint.z));
      addPoint(glm::vec3(bbox.minPoint.x, bbox.maxPoint.y, bbox.maxPoint.z));
      addPoint(glm::vec3(bbox.maxPoint.x, bbox.minPoint.y, bbox.minPoint.z));
      addPoint(glm::vec3(bbox.maxPoint.x, bbox.minPoint.y, bbox.maxPoint.z));
      addPoint(glm::vec3(bbox.maxPoint.x, bbox.maxPoint.y, bbox.minPoint.z));
      addPoint(glm::vec3(bbox.maxPoint.x, bbox.maxPoint.y, bbox.maxPoint.z));
    }

    bool overlaps(const BBox &bbox) const {
      return ((minPoint.x <= bbox.maxPoint.x) && (maxPoint.x >= minPoint.x) &&
              (minPoint.y <= bbox.maxPoint.y) && (maxPoint.y >= minPoint.y) &&
              (minPoint.z <= bbox.maxPoint.z) && (maxPoint.z >= minPoint.z));
    }

    bool contains(const glm::vec3 &point) const {
      return (point.x >= minPoint.x) && (point.x <= maxPoint.x) &&
            (point.y >= minPoint.y) && (point.y <= maxPoint.y) &&
            (point.z >= minPoint.z) && (point.z <= maxPoint.z);
    }

    glm::vec3 size() const {
      return maxPoint - minPoint;
    }

    glm::vec3 center() const {
      return 0.5f * (maxPoint + minPoint);
    }

    BBox &operator=(const BBox &bbox) {
      minPoint = bbox.minPoint;
      maxPoint = bbox.maxPoint;
      return *this;
    }
};

inline BBox operator*(const glm::mat4 &matrix, const BBox &bbox) {
  glm::vec3 points[8];
  points[0] = glm::vec3(bbox.minPoint.x, bbox.minPoint.y, bbox.minPoint.z);
  points[1] = glm::vec3(bbox.minPoint.x, bbox.minPoint.y, bbox.maxPoint.z);
  points[2] = glm::vec3(bbox.minPoint.x, bbox.maxPoint.y, bbox.minPoint.z);
  points[3] = glm::vec3(bbox.minPoint.x, bbox.maxPoint.y, bbox.maxPoint.z);
  points[4] = glm::vec3(bbox.maxPoint.x, bbox.minPoint.y, bbox.minPoint.z);
  points[5] = glm::vec3(bbox.maxPoint.x, bbox.minPoint.y, bbox.maxPoint.z);
  points[6] = glm::vec3(bbox.maxPoint.x, bbox.maxPoint.y, bbox.minPoint.z);
  points[7] = glm::vec3(bbox.maxPoint.x, bbox.maxPoint.y, bbox.maxPoint.z);

  BBox result;
  for (int i = 0; i < 8; ++i) {
    const glm::vec3 &point = points[i];
    const glm::vec4 trPoint4 = matrix * glm::vec4(point.x, point.y, point.z, 1);
    const glm::vec3 trPoint = glm::vec3(trPoint4.x, trPoint4.y, trPoint4.z);
    result.addPoint(trPoint);
  }

  return result;
}

/**
 * Tests if a 2D point is inside a Triangle
 *
 * @param p The point to be tested
 * @param a Vertex A of the triangle
 * @param b Vertex B of the triangle
 * @param c Vertex C of the triangle
 * @return true if the point is inside the triangle, false otherwise
 */ 
bool point2DInsideTriangle(glm::vec2 p, glm::vec2 a, glm::vec2 b, glm::vec2 c);

/**
 * Tests if a 2D Quad intersects a 2D triangle
 *
 * @param mesh Triangle mesh of the model
 * @param t Triangle to test
 * @param quad_min Minimum point of the quad
 * @param quad_max Maximum point of the quad
 * @return true if the triangle and the quad intersect, false otherwise
 */ 
bool testQuadTriangle(TriangleMesh* mesh, Triangle t, glm::vec2 min_point, glm::vec2 max_point);

/**
 * Tests if a 3D Box intersects a 3D triangle
 *
 * @param mesh Triangle mesh of the model
 * @param t Triangle to test
 * @param min_point Minimum point of the Box
 * @param max_point Maximum point of the Box
 * @return true if the triangle and the Box intersect, false otherwise
 */ 
bool testBoxTriangle(TriangleMesh* mesh, Triangle t, glm::vec3 min_point, glm::vec3 max_point, bool test2D = false);

bool rayIntersectsTriangle(glm::vec3 rayOrigin, glm::vec3 rayVector, glm::vec3 v1_tri, glm::vec3 v2_tri, glm::vec3 v3_tri, glm::vec3& outIntersectionPoint);

} // Geo

#endif // GEO_H
