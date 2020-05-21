#include "Geometry.h"
#include <iostream>

namespace Geo {

bool point2DInsideTriangle(glm::vec2 p, glm::vec2 a, glm::vec2 b, glm::vec2 c) {
    float as_x = p.x-a.x;
    float as_y = p.y-a.y;

    bool p_ab = (b.x-a.x)*as_y-(b.y-a.y)*as_x > 0;

    if(((c.x-a.x)*as_y-(c.y-a.y)*as_x > 0) == p_ab) return false;

    if(((c.x-b.x)*(p.y-b.y)-(c.y-b.y)*(p.x-b.x) > 0) != p_ab) return false;

    return true;
}


bool testAABoxAABox_2D(glm::vec2 min_box_1, glm::vec2 max_box_1, glm::vec2 min_box_2, glm::vec2 max_box_2) {
  bool x1 = max_box_1.x >= min_box_2.x;
  bool x2 = max_box_2.x >= min_box_1.x;
  bool y1 = max_box_1.y >= min_box_2.y;
  bool y2 = max_box_2.y >= min_box_1.y;

  return (x1 && x2 && y1 && y2);
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

bool testPointsVsAxis2D(glm::vec2 v1, glm::vec2 v2, glm::vec2 v3, glm::vec2 axis, glm::vec2 h) {
  float p0 = glm::dot(axis, v1);
  float p1 = glm::dot(axis, v2);
  float p2 = glm::dot(axis, v3);

  float min_p = glm::min(p0, glm::min(p1,p2));
  float max_p = glm::max(p0, glm::max(p1,p2));

  float r = h.x * abs(axis.x) + h.y * abs(axis.y);

  if (min_p > r || max_p < -r ) return false;
  return true;
}

bool testPointsVsAxis(glm::vec3 v1, glm::vec3 v2, glm::vec3 v3, glm::vec3 axis, glm::vec3 h) {
  float p0 = glm::dot(axis, v1);
  float p1 = glm::dot(axis, v2);
  float p2 = glm::dot(axis, v3);

  float min_p = glm::min(p0, glm::min(p1,p2));
  float max_p = glm::max(p0, glm::max(p1,p2));

  float r = h.x * abs(axis.x) + h.y * abs(axis.y) + h.z * abs(axis.z);

  if (min_p > r || max_p < -r ) return false;
  return true;
}

bool testQuadTriangle(TriangleMesh* mesh, Triangle *t, glm::vec2 &min_point, glm::vec2 &max_point) {
  Geo::BBox box;

  glm::vec3 min_point_3d(0, min_point.x, min_point.y);
  glm::vec3 max_point_3d(0, max_point.x, max_point.y);

  box.addPoint(min_point_3d);
  box.addPoint(max_point_3d);

  glm::vec3 v1o, v2o, v3o;
  std::vector<glm::vec3> &vertices = mesh->getVertices();
  
  v1o = vertices[t->v1];
  v2o = vertices[t->v2];
  v3o = vertices[t->v3];

  const glm::vec3 center_box = (box.minPoint + box.maxPoint) / 2.0f;

  glm::vec3 v1 = v1o - center_box;
  glm::vec3 v2 = v2o - center_box;
  glm::vec3 v3 = v3o - center_box;
  
  // We ignore the X direction
  v1.x = 0;    
  v2.x = 0;    
  v3.x = 0;    

  glm::vec2 v1_2d = glm::vec2(v1.y, v1.z);
  glm::vec2 v2_2d = glm::vec2(v2.y, v2.z);
  glm::vec2 v3_2d = glm::vec2(v3.y, v3.z);

  if (!testAABoxAABox_2D(glm::vec2(t->tri_bbox.minPoint.y, t->tri_bbox.minPoint.z), 
                         glm::vec2(t->tri_bbox.maxPoint.y, t->tri_bbox.maxPoint.z),
                         glm::vec2(box.minPoint.y, box.minPoint.z),
                         glm::vec2(box.maxPoint.y, box.maxPoint.z)))
    return false;

  const glm::vec2 f1 = v2_2d-v1_2d;
  const glm::vec2 f2 = v3_2d-v2_2d;
  const glm::vec2 f3 = v1_2d-v3_2d;

  glm::vec3 e = box.maxPoint - center_box; // Compute positive extents
  glm::vec2 e_2d = glm::vec2(e.y, e.z);

  glm::vec2 n1 = glm::vec2(f1.y, -f1.x);
  glm::vec2 n2 = glm::vec2(f2.y, -f2.x);
  glm::vec2 n3 = glm::vec2(f3.y, -f3.x);

    if (!testPointsVsAxis2D(v1_2d, v2_2d, v3_2d, n1, e_2d) || 
        !testPointsVsAxis2D(v1_2d, v2_2d, v3_2d, n2, e_2d) || 
        !testPointsVsAxis2D(v1_2d, v2_2d, v3_2d, n3, e_2d))
      {
        return false;
      }
  return true;
}


bool testBoxTriangle(TriangleMesh* mesh, Triangle *t, glm::vec3 &min_point, glm::vec3 &max_point) {

  Geo::BBox box;

  box.addPoint(min_point);
  box.addPoint(max_point);

  glm::vec3 v1o, v2o, v3o;
  std::vector<glm::vec3> &vertices = mesh->getVertices();
  
  v1o = vertices[t->v1];
  v2o = vertices[t->v2];
  v3o = vertices[t->v3];


  const glm::vec3 center_box = (box.minPoint + box.maxPoint) / 2.0f;

  glm::vec3 v1 = v1o - center_box;
  glm::vec3 v2 = v2o - center_box;
  glm::vec3 v3 = v3o - center_box;

  // If their bounding boxes don't intersect, they definitely don't intersect
  if (!testAABoxAABox(t->tri_bbox, box)) return false;
  
  const glm::vec3 f1 = v2-v1;
  const glm::vec3 f2 = v3-v2;
  const glm::vec3 f3 = v1-v3;

  glm::vec3 tri_normal = glm::normalize(glm::cross(f1, f2));

  glm::vec3 vector1, vector2;
  
  float planeDistance = - (tri_normal.x * v1.x + tri_normal.y * v1.y + tri_normal.z + v1.z);

  if(tri_normal.x >= 0) {
    vector1.x = min_point.x;
    vector2.x = max_point.x;
  } else {
    vector1.x = max_point.x;
    vector2.x = min_point.x;
  }
  if(tri_normal.y >= 0) {
    vector1.y = min_point.y;
    vector2.y = max_point.y;
  } else {
    vector1.y = max_point.y;
    vector2.y = min_point.y;
  }
  if(tri_normal.z >= 0) {
    vector1.z = min_point.z;
    vector2.z = max_point.z;
  } else {
    vector1.z = min_point.z;
    vector2.z = max_point.z;
  }
  float posSide = (tri_normal.x * vector2.x)+(tri_normal.y * vector2.y)+(tri_normal.y * vector2.y)+planeDistance;
  float negSide = (tri_normal.x * vector1.x)+(tri_normal.y * vector1.y)+(tri_normal.y * vector1.y)+planeDistance;
  if(posSide <=  0 && negSide >= 0) {
    //box not intersects
    return false;
  }

  glm::vec3 e = box.maxPoint - center_box; // Compute positive extents
  const glm::vec3 e0 = glm::vec3(1,0,0);
  const glm::vec3 e1 = glm::vec3(0,1,0);
  const glm::vec3 e2 = glm::vec3(0,0,1);

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

IntersectionResult rayIntersectsTriangle(glm::vec3 &rayOrigin, glm::vec3 &rayVector, const glm::vec3 &v1_tri, const glm::vec3 &v2_tri, const glm::vec3 &v3_tri, double &threshold, glm::vec3& outIntersectionPoint) {
  glm::vec3 edge1, edge2;
  float f,u,v,a;
  glm::vec3 q,h,s;
  edge1 = v2_tri - v1_tri;
  edge2 = v3_tri - v1_tri;

  glm::vec3 n = glm::cross(edge2, edge1);

  // Check if ray and triangle are parallel
  float nDotRay = glm::dot(n, rayVector);
  if (abs(nDotRay) <= threshold) {
    // They are parallel, but are they invalid?
    glm::vec3 barycenter = (v1_tri + v2_tri + v3_tri) / 3.0f;
    if (distPointLine(barycenter, rayVector, rayOrigin) < threshold) {
      // Ray is invalid
      return IntersectionResult::INVALID;
    } else {
      // It's not invalid, just doesn't intersect
      return IntersectionResult::NOT_INTERSECTS;
    }
  }

  h = glm::cross(rayVector, edge2);
  a = glm::dot(edge1, h);

  f = 1.0/a;
  s = rayOrigin - v1_tri;
  u = f * glm::dot(s, h);
  // barycentric coordinate U
  if (u < 0.0 || u > 1.0) {
    // The point is outside the triangle
    return IntersectionResult::NOT_INTERSECTS;
  }

  q = glm::cross(s, edge1);
  v = f * glm::dot(rayVector, q);
  // barycentric coordinate V
  if (v < 0.0 || u + v > 1.0){
    // The point is outside the triangle
    return IntersectionResult::NOT_INTERSECTS;

  }

  // At this stage we can compute t to find out where the intersection point is on the line.
  float t = f * glm::dot(edge2, q);
  if (t > threshold) { // ray intersection
    // glm::vec3 test = v*v1_tri + u*v2_tri + (1-u-v)*v3_tri;
    // std::cout << test.x << " " << test.y << " " << test.z << std::endl;
    outIntersectionPoint = rayOrigin + rayVector * t;
    // outIntersectionPoint = v*v1_tri + u*v2_tri + (1-u-v)*v3_tri;
    // std::cout << "orig: " << outIntersectionPoint.x << " " << outIntersectionPoint.y << " " << outIntersectionPoint.z << std::endl;
    // Now we check if the intersection is too close to an edge or vertex of the triangle
    float w = 1.0f - u - v;
    if ((u < threshold && v < threshold && w >= threshold) ||
        (v < threshold && w < threshold && u >= threshold) ||
        (w < threshold && u < threshold && v >= threshold)) {
      // Too close to a vertex
      return IntersectionResult::INVALID;
    }

    if ((u < threshold && v >= threshold && w >= threshold) ||
        (v < threshold && w >= threshold && u >= threshold) ||
        (w < threshold && u >= threshold && v >= threshold)) {
      // Too close to an edge
      return IntersectionResult::INVALID;
    }

    return IntersectionResult::INTERSECTS;
  }
  else // This means that there is a line intersection but not a ray intersection. 
    return IntersectionResult::NOT_INTERSECTS;
}

double distPointLine(glm::vec3 point, glm::vec3 lineDir, glm::vec3 pointInLine) {
  return glm::length(glm::cross(pointInLine-point, lineDir)) / glm::length(lineDir);
}

bool isRayInvalid(glm::vec3 rayOrigin, glm::vec3 rayVector, glm::vec3 v1_tri, glm::vec3 v2_tri, glm::vec3 v3_tri, double threshold) {
  glm::vec3 normal = glm::normalize(glm::cross(v3_tri - v2_tri, v3_tri - v1_tri));

  if (abs(glm::dot(normal, rayVector)) <= 0.0001 ) {
    // Ray is invalid if it's too close
    // It's too close to the triangle if it's too close to one of its points or one of its edges
    if (distPointLine(v1_tri, rayVector, rayOrigin) < threshold) return true;
    if (distPointLine(v2_tri, rayVector, rayOrigin) < threshold) return true;
    if (distPointLine(v3_tri, rayVector, rayOrigin) < threshold) return true;
    // TODO distance from ray to each edge
  }

  return false; 
}



} //GEO