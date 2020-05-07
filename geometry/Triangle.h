#ifndef TRI_H
#define TRI_H


class TriangleMesh;

#include "BBox.h"

class Triangle {

  public:
    Triangle(unsigned int v1_, unsigned int v2_, unsigned int v3_);
    ~Triangle();
    // unsigned int getV1();
    // unsigned int getV2();
    // unsigned int getV3();
    unsigned int v1, v2, v3;
    double min_x;
    Geo::BBox tri_bbox;
    void saveMinX(TriangleMesh* m);
    void calculateTriBbox(std::vector<glm::vec3> &vertices);

    
    bool operator < (const Triangle& tr2) const {
      return (min_x < tr2.min_x);
    }
  private:
    
    
};

#endif