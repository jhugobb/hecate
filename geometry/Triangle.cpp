#include "Triangle.h"

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