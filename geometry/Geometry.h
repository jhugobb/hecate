#ifndef GEO_H
#define GEO_H

#include <glm/glm.hpp>

#include "../TriangleMesh.h"
#include "Triangle.h"

namespace Geo {

  enum IntersectionResult {
    NOT_INTERSECTS,
    INTERSECTS,
    INVALID
  };

const float EPSILON = 0.0000001;



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
bool testQuadTriangle(TriangleMesh* mesh, Triangle *t, glm::vec2 &min_point, glm::vec2 &max_point);

/**
 * Tests if a 3D Box intersects a 3D triangle
 *
 * @param mesh Triangle mesh of the model
 * @param t Triangle to test
 * @param min_point Minimum point of the Box
 * @param max_point Maximum point of the Box
 * @return true if the triangle and the Box intersect, false otherwise
 */ 
bool testBoxTriangle(TriangleMesh* mesh, Triangle *t, glm::vec3 &min_point, glm::vec3 &max_point);

IntersectionResult rayIntersectsTriangle(glm::vec3 rayOrigin, glm::vec3 rayVector, glm::vec3 v1_tri, glm::vec3 v2_tri, glm::vec3 v3_tri, double treshold, glm::vec3& outIntersectionPoint);

double distPointLine(glm::vec3 point, glm::vec3 lineDir, glm::vec3 pointInLine);

bool isRayInvalid(glm::vec3 rayOrigin, glm::vec3 rayVector, glm::vec3 v1_tri, glm::vec3 v2_tri, glm::vec3 v3_tri, double threshold);

} // Geo

#endif // GEO_H
