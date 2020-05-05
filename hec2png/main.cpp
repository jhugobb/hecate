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

void readHEC_normal(char* filename) {
  streampos size;
  char * memblock;

  ifstream file (filename, ios::in|ios::binary|ios::ate);
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
    assert((int) voxels.size() == resolution*resolution*resolution);

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
}

void readHEC_RLE_N_8(char* filename) {
  streampos size;
  char * memblock;

  ifstream file (filename, ios::in|ios::binary|ios::ate);
  if (file.is_open())
  {
    size = file.tellg();
    memblock = new char [size];
    file.seekg (0, ios::beg);
    file.read (memblock, size);
    file.close();
    char c_1, c_2;
    c_1 = memblock[0];
    c_2 = memblock[1];
    bitset<8> res_1(c_1);
    bitset<8> res_2(c_2);

    bitset<16> resolution_b;

    for (uint i = 0; i < 8; i++) {
      resolution_b[i] = res_1[i];
      resolution_b[i+8] = res_2[i];
    }
    int resolution = resolution_b.to_ulong();

    cout << "Resolution: " << resolution  << endl;
    
    std::vector<Color> voxels;
    std::vector<int> count;

    for (int i = 2; i < size; i++) {
      char c = memblock[i];
      bitset<8> bits(c);
      if (bits[7] == 0 && bits[6] == 0) voxels.push_back(Color::W);
      else if (bits[7] == 0 && bits[6] == 1) voxels.push_back(Color::B);
      else if (bits[7] == 1 && bits[6] == 0) voxels.push_back(Color::G);
      else assert(false);

      bits.set(7,0);
      bits.set(6,0);

      count.push_back(bits.to_ulong()+1);
    }
    
    int sum = 0;
    for (int cou : count) {
      sum += cou;
    }
    cout << "voxels size:" << sum << endl;
    assert(sum == resolution*resolution*resolution);


    int runs_idx = 0;
    int n_voxels_written = 0;
    for (int y = 0; y < resolution; y++) {
      std::vector<unsigned char> image(resolution*resolution*4);
      for (int z = 0; z < resolution; z++)
      for (int x = 0; x < resolution; x++) {
        unsigned int value = 0;
        if (voxels[runs_idx] == Color::G) {
          value = 128;
        } else if (voxels[runs_idx] == Color::W) {
          value = 255;
        }
        image[4 * resolution * z + 4 * x + 0] = value;
        image[4 * resolution * z + 4 * x + 1] = value;
        image[4 * resolution * z + 4 * x + 2] = value;
        image[4 * resolution * z + 4 * x + 3] = 255;
        n_voxels_written++;
        if (n_voxels_written+1 > count[runs_idx]) {
          runs_idx++;
          n_voxels_written = 0;
        }
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
}

void readHEC_RLE_N_16(char* filename) {
  streampos size;
  char * memblock;

  ifstream file (filename, ios::in|ios::binary|ios::ate);
  if (file.is_open())
  {
    size = file.tellg();
    memblock = new char [size];
    file.seekg (0, ios::beg);
    file.read (memblock, size);
    file.close();
    const char c_1 = memblock[0];
    const char c_2 = memblock[1];
    bitset<8> res_1(c_1);
    bitset<8> res_2(c_2);

    bitset<16> resolution_b;

    for (uint i = 0; i < 8; i++) {
      resolution_b[i] = res_1[i];
      resolution_b[i+8] = res_2[i];
    }
    int resolution = resolution_b.to_ulong();

    cout << "Resolution: " << resolution  << endl;
    
    std::vector<Color> voxels;
    std::vector<int> count;

    for (int i = 2; i < size; i+=2) {
      char c_first = memblock[i];
      char c_second = memblock[i+1];
      bitset<8> bits_first(c_first);
      bitset<8> bits_second(c_second);
      bitset<16> bits;

      for (uint idx = 0; idx < 8; idx++) {
        bits[idx] = bits_first[idx];
        bits[idx+8] = bits_second[idx];
      }
      if (bits[15] == 0 && bits[14] == 0) voxels.push_back(Color::W);
      else if (bits[15] == 0 && bits[14] == 1) voxels.push_back(Color::B);
      else if (bits[15] == 1 && bits[14] == 0) voxels.push_back(Color::G);
      else assert(false);

      bits.set(15,0);
      bits.set(14,0);

      count.push_back(bits.to_ulong()+1);
    }
    
    int sum = 0;
    for (int cou : count) {
      sum += cou;
    }
    cout << "voxels size:" << sum << endl;
    assert(sum == resolution*resolution*resolution);


    int runs_idx = 0;
    int n_voxels_written = 0;
    for (int y = 0; y < resolution; y++) {
      std::vector<unsigned char> image(resolution*resolution*4);
      for (int z = 0; z < resolution; z++)
      for (int x = 0; x < resolution; x++) {
        unsigned int value = 0;
        if (voxels[runs_idx] == Color::G) {
          value = 128;
        } else if (voxels[runs_idx] == Color::W) {
          value = 255;
        }
        image[4 * resolution * z + 4 * x + 0] = value;
        image[4 * resolution * z + 4 * x + 1] = value;
        image[4 * resolution * z + 4 * x + 2] = value;
        image[4 * resolution * z + 4 * x + 3] = 255;
        n_voxels_written++;
        if (n_voxels_written+1 > count[runs_idx]) {
          runs_idx++;
          n_voxels_written = 0;
        }
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
}

void readHEC_RLE_A_8(char* filename) {
  streampos size;
  char * memblock;

  ifstream file (filename, ios::in|ios::binary|ios::ate);
  if (file.is_open())
  {
    size = file.tellg();
    memblock = new char [size];
    file.seekg (0, ios::beg);
    file.read (memblock, size);
    file.close();
    char c_1, c_2;
    c_1 = memblock[0];
    c_2 = memblock[1];
    bitset<8> res_1(c_1);
    bitset<8> res_2(c_2);

    bitset<16> resolution_b;

    for (uint i = 0; i < 8; i++) {
      resolution_b[i] = res_1[i];
      resolution_b[i+8] = res_2[i];
    }
    int resolution = resolution_b.to_ulong();

    cout << "Resolution: " << resolution  << endl;
    
    std::vector<Color> voxels;
    std::vector<int> count;

    Color order[2] = {W, B};
    Color color = W;
    int n_times_changed = 0;
    bool writeGray = false;
    bool wrote1s = false;
    int last = -1;
    int n_gray = 0;
    for (int i = 2; i < size; i++) {
      char c = memblock[i];
      bitset<8> bits(c);
      // cout << bits.to_ulong() << endl;
      if (bits.all()) {
        // if (last != 0)
        n_times_changed++;
        color = order[n_times_changed % 2];
        // writeGray = false;
        // last = 255;
        continue;
      } else if (bits.to_ulong() == 254) {
        count.push_back(254);
        voxels.push_back(color);
      } else {
        if (bits.any()) {
          count.push_back(bits.to_ulong());
          voxels.push_back(color);
        }
        if (i != size-1) {
          // don't insert gray at the end
          count.push_back(1);
          voxels.push_back(Color::G);
          n_gray++;
        }
        n_times_changed++;
        color = order[n_times_changed % 2];
      }

      // if (last == 0)
      //   writeGray = true;

      // if (writeGray) {
      //   voxels.push_back(Color::G);
      //   count.push_back(1);
      //   n_gray++;
      // }

      // if (bits.any()) {
      //   // if (last == 0) n_times_changed++;
      //   count.push_back(bits.to_ulong());
      //   voxels.push_back(order[n_times_changed % 2]);
      //   // n_times_changed = 0;
      // }
      // n_times_changed++;

      // writeGray = true;
      // last = bits.to_ulong();
    }
    int sum = 0;
    for (int cou : count) {
      sum += cou;
    }
    cout << "voxels size:" << sum << endl;
    cout << "difference: " << sum - resolution*resolution*resolution << endl;
    cout << "total grays: " << n_gray << endl;

    // for (int idx_pls = 0; idx_pls < voxels.size(); idx_pls++) {
    //   cout << "Color: " << voxels[idx_pls] << " Num: " << count[idx_pls] << endl;;
    // }
    // assert(sum == resolution*resolution*resolution);

    int runs_idx = 0;
    int n_voxels_written = 0;
    for (int y = 0; y < resolution; y++) {
      std::vector<unsigned char> image(resolution*resolution*4);
      for (int z = 0; z < resolution; z++)
      for (int x = 0; x < resolution; x++) {
        unsigned int value = 0;
        if (voxels[runs_idx] == Color::G) {
          value = 128;
        } else if (voxels[runs_idx] == Color::W) {
          value = 255;
        }
        image[4 * resolution * z + 4 * x + 0] = value;
        image[4 * resolution * z + 4 * x + 1] = value;
        image[4 * resolution * z + 4 * x + 2] = value;
        image[4 * resolution * z + 4 * x + 3] = 255;
        n_voxels_written++;
        if (n_voxels_written+1 > count[runs_idx]) {
          runs_idx++;
          n_voxels_written = 0;
        }
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
}


int main(int argc, char **argv) {
  
  if (argc != 3) {
    cout << "Usage: ./hec2png /path/to/hec/file type_of_hec" << endl; 
    return 0;
  }
  // Types of HEC
  // 0 -> Normal
  // 1 -> RLE Naive 8bit
  // 2 -> RLE Naive 16bit
  // 3 -> RLE Alternated 8bit
  // 4 -> RLE Alternated 8bit
  // 5 -> Modulus Coding

  int type_of_hec = std::stoi(argv[2]);
  switch (type_of_hec) {
    case 0:
      readHEC_normal(argv[1]);
      break;
    case 1:
      readHEC_RLE_N_8(argv[1]);
      break;
    case 2:  
      readHEC_RLE_N_16(argv[1]);
      break;
    case 3:
      readHEC_RLE_A_8(argv[1]);
    default:
      break;
  }
  
  return 0;
}