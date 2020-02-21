#include "Geometry.h"

namespace Geo {
bool point2DInsideTriangle(glm::vec2 p, glm::vec2 a, glm::vec2 b, glm::vec2 c) {
    float as_x = p.x-a.x;
    float as_y = p.y-a.y;

    bool p_ab = (b.x-a.x)*as_y-(b.y-a.y)*as_x > 0;

    if((c.x-a.x)*as_y-(c.y-a.y)*as_x > 0 == p_ab) return false;

    if((c.x-b.x)*(p.y-b.y)-(c.y-b.y)*(p.x-b.x) > 0 != p_ab) return false;

    return true;
}

bool triangleQuadintersectionTest(TriangleMesh* mesh, Triangle* t, glm::vec2 quad_min, glm::vec2 quad_max) {


  // Create 2d bbox for the triangle, and then the circumscribed circle of that box and the quad. 
  // If the distance of their centers is less than the sum of their radiuses, then they do not intersect.

  Geo::BBox tri_box;

  // project triangle into z = 0 plane
  std::vector<glm::vec3> vertices = mesh->getVertices();
  glm::vec3 v1 = vertices[t->getV1()];
  glm::vec3 v2 = vertices[t->getV2()];
  glm::vec3 v3 = vertices[t->getV3()];

  tri_box.addPoint(v1);
  tri_box.addPoint(v2);
  tri_box.addPoint(v3);

  glm::vec2 tri_box_min_2d = glm::vec2(tri_box.minPoint.x, tri_box.minPoint.y);
  glm::vec2 tri_box_max_2d = glm::vec2(tri_box.maxPoint.x, tri_box.maxPoint.y);
  glm::vec2 center_tri_box =  (tri_box_min_2d + tri_box_max_2d) / 2.0f;
  float radius_tri_box = glm::distance(tri_box_min_2d, center_tri_box);

  glm::vec2 center_quad = (quad_min + quad_max) / 2.0f;
  float radius_quad = glm::distance(quad_min, center_quad);

  if (glm::distance(center_tri_box, center_quad) > radius_quad + radius_tri_box) {
    return false;
  }

  // Check if any of the points of the triangle lie inside the quad
  if (v1.x >= quad_min.x && v1.x <= quad_max.x && v1.y >= quad_min.y && v1.y <= quad_max.y ||
      v2.x >= quad_min.x && v2.x <= quad_max.x && v2.y >= quad_min.y && v2.y <= quad_max.y ||
      v3.x >= quad_min.x && v3.x <= quad_max.x && v3.y >= quad_min.y && v3.y <= quad_max.y)
  {
      return true;
  }
  
  // Check if any of the points of the quad lie inside the triangle

  if (point2DInsideTriangle(quad_min, v1, v2, v3)                          ||
      point2DInsideTriangle(quad_max, v1, v2, v3)                          ||
      point2DInsideTriangle(glm::vec2(quad_min.x, quad_max.y), v1, v2, v3) ||
      point2DInsideTriangle(glm::vec2(quad_max.x, quad_min.y), v1, v2, v3)) 
  {
      return true;
  }

  return false;

}

bool testAABoxAABox(Geo::BBox box1, Geo::BBox box2) {
  bool x1 = box1.maxPoint.x >= box2.minPoint.x;
  bool x2 = box2.maxPoint.x >= box1.minPoint.x;
  bool y1 = box1.maxPoint.y >= box2.minPoint.y;
  bool y2 = box2.maxPoint.y >= box1.minPoint.y;
  bool z1 = box1.maxPoint.z >= box2.minPoint.z;
  bool z2 = box2.maxPoint.z >= box1.minPoint.z;

  return (x1 && x2 && y1 && y2 && z1 && z2);
}

bool testPointsVsAxis(glm::vec3 v1, glm::vec3 v2, glm::vec3 v3, glm::vec3 axis, glm::vec3 h) {
  float p0 = glm::dot(axis, v1);
  float p1 = glm::dot(axis, v2);
  float p2 = glm::dot(axis, v3);

  float min_p = glm::min(p0, glm::min(p1,p2));
  float max_p = glm::max(p0, glm::max(p1,p2));

  const glm::vec3 e0 = glm::vec3(1,0,0);
  const glm::vec3 e1 = glm::vec3(0,1,0);
  const glm::vec3 e2 = glm::vec3(0,0,0);

  float r = h.x * abs(axis.x) + h.y * abs(axis.y) + h.z * abs(axis.z);

  if (min_p > r || max_p < -r ) return false;
  return true;
}


bool testBoxTriangle(TriangleMesh* mesh, Triangle t, glm::vec3 min_point, glm::vec3 max_point, bool test2D) {

  Geo::BBox box;

  box.addPoint(min_point);
  box.addPoint(max_point);

  glm::vec3 v1o, v2o, v3o;
  std::vector<Triangle> tris = mesh->getTriangles();
  std::vector<glm::vec3> vertices = mesh->getVertices();
  
  v1o = vertices[t.getV1()];
  v2o = vertices[t.getV2()];
  v3o = vertices[t.getV3()];


  const glm::vec3 center_box = (box.minPoint + box.maxPoint) / 2.0f;

  glm::vec3 v1 = v1o - center_box;
  glm::vec3 v2 = v2o - center_box;
  glm::vec3 v3 = v3o - center_box;
  
  if (test2D) {
    v1.z = 0;    
    v2.z = 0;    
    v3.z = 0;    
  }

  Geo::BBox tri_bbox;

  tri_bbox.addPoint(v1o);
  tri_bbox.addPoint(v2o);
  tri_bbox.addPoint(v3o);

  if (testAABoxAABox(tri_bbox, box)) return true;

  const glm::vec3 f1 = v2-v1;
  const glm::vec3 f2 = v3-v2;
  const glm::vec3 f3 = v1-v3;

  glm::vec3 tri_normal = glm::normalize(glm::cross(f1, f2));


  // Convert AABB to center-extents representation
  glm::vec3 e = box.maxPoint - center_box; // Compute positive extents

  // Compute the projection interval radius of b onto L(t) = b.c + t * p.n
  float r = e[0]*glm::abs(tri_normal.x) + e[1]*glm::abs(tri_normal.y) + e[2]*glm::abs(tri_normal.z);

  float d = tri_normal.x * v1.x + tri_normal.y * v1.y + tri_normal.z * v1.z;

  // Compute distance of box center from plane
  float s = glm::dot(tri_normal, center_box) - d;
  
  // Intersection occurs when distance s falls within [-r,+r] interval
  if (abs(s) <= r) return true;

  const glm::vec3 e0 = glm::vec3(1,0,0);
  const glm::vec3 e1 = glm::vec3(0,1,0);
  const glm::vec3 e2 = glm::vec3(0,0,0);

  if (!testPointsVsAxis(v1, v2, v3, glm::cross(e0, f1), e) || 
      !testPointsVsAxis(v1, v2, v3, glm::cross(e1, f1), e) ||
      !testPointsVsAxis(v1, v2, v3, glm::cross(e2, f1), e) ||
      !testPointsVsAxis(v1, v2, v3, glm::cross(e0, f2), e) ||
      !testPointsVsAxis(v1, v2, v3, glm::cross(e1, f2), e) ||
      !testPointsVsAxis(v1, v2, v3, glm::cross(e2, f2), e) ||
      !testPointsVsAxis(v1, v2, v3, glm::cross(e0, f3), e) ||
      !testPointsVsAxis(v1, v2, v3, glm::cross(e1, f3), e) ||
      !testPointsVsAxis(v1, v2, v3, glm::cross(e2, f3), e))
      {
        return false;
      }

  return true;
}

bool isPointInQuad(glm::vec2 min_point, glm::vec2 max_point, glm::vec2 query_point) {
  return min_point.x <= query_point.x && min_point.y <= query_point.y 
      && max_point.x >= query_point.x && max_point.y >= query_point.y;
}

bool rayIntersectsTriangle(glm::vec3 rayOrigin, glm::vec3 rayVector, glm::vec3 v1_tri, glm::vec3 v2_tri, glm::vec3 v3_tri, glm::vec3& outIntersectionPoint) {
  glm::vec3 edge1, edge2, h, s, q;
  float a,f,u,v;
  edge1 = v2_tri - v1_tri;
  edge2 = v3_tri - v1_tri;
  h = glm::cross(rayVector, edge2);
  a = glm::dot(edge1, h);
  if (a > -Geo::EPSILON && a < Geo::EPSILON)
    return false;    // This ray is parallel to this triangle.
  f = 1.0/a;
  s = rayOrigin - v1_tri;
  u = f * glm::dot(s, h);
  if (u < 0.0 || u > 1.0)
    return false;
  q = glm::cross(s, edge2);
  v = f * glm::dot(rayVector, q);
  if (v < 0.0 || u + v > 1.0)
      return false;
  // At this stage we can compute t to find out where the intersection point is on the line.
  float t = f * glm::dot(edge2, q);
  if (t > EPSILON) { // ray intersection
    outIntersectionPoint = rayOrigin + rayVector * t;
    return true;
  }
  else // This means that there is a line intersection but not a ray intersection.
    return false;
}



} //GEO