#ifndef TRI_H
#define TRI_H

class Triangle {

  public:
    Triangle(unsigned int v1_, unsigned int v2_, unsigned int v3_);
    ~Triangle();
    unsigned int getV1();
    unsigned int getV2();
    unsigned int getV3();

  private:
    unsigned int v1, v2, v3;
};

#endif