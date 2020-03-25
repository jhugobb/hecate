#ifndef TRI_H
#define TRI_H


class TriangleMesh;

class Triangle {

  public:
    Triangle(unsigned int v1_, unsigned int v2_, unsigned int v3_);
    ~Triangle();
    unsigned int getV1();
    unsigned int getV2();
    unsigned int getV3();
    double min_x;
    void saveMinX(TriangleMesh* m);
    
    bool operator < (const Triangle& tr2) const {
      return (min_x < tr2.min_x);
    }
  private:
    unsigned int v1, v2, v3;
    
};

#endif