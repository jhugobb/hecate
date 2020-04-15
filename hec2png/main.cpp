#include <iostream>
#include <fstream>
#include <math.h>
#include <bitset>
#include <cassert>

#include "lodepng.h"

using namespace std;

enum Color {
  W,
  B,
  G
};

int main(int argc, char **argv) {
  
  if (argc > 2) {
    cout << "Usage: ./hec2png /path/to/hec/file" << endl; 
    return 0;
  }


  streampos size;
  char * memblock;

  ifstream file (argv[1], ios::in|ios::binary|ios::ate);
  if (file.is_open())
  {
    size = file.tellg();
    memblock = new char [size];
    file.seekg (0, ios::beg);
    file.read (memblock, size);
    file.close();
    
    int resolution = cbrt(size * 4);

    cout << "Resolution: " << resolution  << endl;
    
    std::vector<Color> voxels;

    for (int i = 0; i < size; i++) {
      char c = memblock[i];
      bitset<8> bits(c);
      for (int idx = 0; idx < 8; idx+=2) {
        if (bits[idx] == 0 && bits[idx+1] == 0) voxels.push_back(Color::W);
        else if (bits[idx] == 0 && bits[idx+1] == 1) voxels.push_back(Color::B);
        else if (bits[idx] == 1 && bits[idx+1] == 0) voxels.push_back(Color::G);
        else assert(false);
      }
    }

    cout << "voxels size:" << voxels.size() << endl;
    assert(voxels.size() == resolution*resolution*resolution);

    for (int y = 0; y < resolution; y++) {
      std::vector<unsigned char> image(resolution*resolution*4);
      for (int z = 0; z < resolution; z++)
      for (int x = 0; x < resolution; x++) {
        unsigned int value = 0;
        if (voxels[y*resolution*resolution + z*resolution+x] == Color::G) {
          value = 128;
        } else if (voxels[y*resolution*resolution + z*resolution+x] == Color::W) {
          value = 255;
        }
        image[4 * resolution * z + 4 * x + 0] = value;
        image[4 * resolution * z + 4 * x + 1] = value;
        image[4 * resolution * z + 4 * x + 2] = value;
        image[4 * resolution * z + 4 * x + 3] = 255;
      }
      std::string filename = "slices/slice_";
      filename.append(std::to_string(y));
      filename.append(".png");
      unsigned error = lodepng::encode(filename, &image[0], resolution, resolution);

      if (error) {
        cout << "Error Writing Image" << endl;
        assert(false);
      }
    }
    
    delete[] memblock;
  } else cout << "Unable to open file";
  
  
  
  
  return 0;
}