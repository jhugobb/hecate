#include <iostream>
#include <fstream>
#include <math.h>
#include <bitset>
#include <cassert>
#include <boost/filesystem.hpp>

#include "lodepng.h"

using namespace std;
using namespace boost::filesystem;

enum Color {
  W,
  B,
  G
};

int readResolution(char mem1, char mem2) {

  bitset<8> res_1(mem1);
  bitset<8> res_2(mem2);

  bitset<16> resolution_b;

  for (uint i = 0; i < 8; i++) {
    resolution_b[i] = res_1[i];
    resolution_b[i+8] = res_2[i];
  }
  
  return resolution_b.to_ulong();
}

void write_slice_PNG(int resolution, std::vector<Color> &voxels, path& entry) {
  std::vector<unsigned char> image(resolution*resolution*4);
  for (int z = 0; z < resolution; z++)
  for (int x = 0; x < resolution; x++) {
    unsigned int value = 0;
    if (voxels[z*resolution+x] == Color::G) {
      value = 128;
    } else if (voxels[z*resolution+x] == Color::W) {
      value = 255;
    }
    image[4 * resolution * z + 4 * x + 0] = value;
    image[4 * resolution * z + 4 * x + 1] = value;
    image[4 * resolution * z + 4 * x + 2] = value;
    image[4 * resolution * z + 4 * x + 3] = 255;
  }
  std::string filename = "slices/slice_";
  filename.append(entry.stem().string());
  filename.append(".png");
  unsigned error = lodepng::encode(filename, &image[0], resolution, resolution);

  if (error) {
    std::cout << "Error Writing Image" << endl;
    assert(false);
  }
}

void write_slice_PNG_RLE(int resolution, std::vector<Color> &voxels, std::vector<int> &count, path& entry) {
  int runs_idx = 0;
  int n_voxels_written = 0;
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
  filename.append(entry.stem().string());
  filename.append(".png");
  unsigned error = lodepng::encode(filename, &image[0], resolution, resolution);

  if (error) {
    cout << "Error Writing Image" << endl;
    assert(false);
  }
}

void readHEC_normal(char* filename, bool write_slices) {
  streampos size;
  char * memblock;
  std::vector<path> paths = std::vector<path>(directory_iterator(filename),boost::filesystem::directory_iterator());
  std::sort(paths.begin(), paths.end());
  bool already_have_resolution = false;
  int resolution = 0;
  double max_time = std::numeric_limits<double>::lowest();
  double sum_time = 0.0;

  for (uint path_idx = 0; path_idx < paths.size(); path_idx++) {
    path& entry = paths[path_idx];
    // std::cout << entry.path().string() << '\n';

    std::ifstream file (entry.string(), ios::in|ios::binary|ios::ate);
    
    timespec begin, end;
    clock_gettime(CLOCK_REALTIME, &begin);

    if (file.is_open())
    {
      size = file.tellg();
      memblock = new char [size];
      file.seekg (0, ios::beg);
      file.read (memblock, size);
      file.close();
      if (!already_have_resolution) {
        resolution = readResolution(memblock[0], memblock[1]);
        std::cout << "Resolution: " << resolution  << endl;
        already_have_resolution = true;
      }

      std::vector<Color> voxels;
      bitset<8> bits;
      for (int i = 2; i < size; i++) {
        bits = bitset<8> (memblock[i]);
        for (int idx = 0; idx < 8; idx+=2) {
          if (bits[idx] == 0 && bits[idx+1] == 0) voxels.push_back(Color::W);
          else if (bits[idx] == 0 && bits[idx+1] == 1) voxels.push_back(Color::B);
          else if (bits[idx] == 1 && bits[idx+1] == 0) voxels.push_back(Color::G);
          else assert(false);
        }
      }

      clock_gettime(CLOCK_REALTIME, &end);
      double time = end.tv_sec - begin.tv_sec + ((end.tv_nsec - begin.tv_nsec) / 1E9);
      sum_time += time;
      if (max_time < time) max_time = time; 
      std::cout << "Decoded " << entry.stem().string() << " -> "
                << time << " s." << std::endl;

      // std::cout << "voxels size:" << voxels.size() << endl;
      assert((int) voxels.size() == resolution*resolution);

      if (write_slices) {
        write_slice_PNG(resolution, voxels, entry);
      }
      delete[] memblock;
    } else std::cout << "Unable to open file";
  }
  std::cout << "Max time of slicing: " << max_time << " s." << endl;
  std::cout << "Total time of slicing: " << sum_time << " s." << endl;
}

void readHEC_RLE_N_8(char* filename, bool write_slices) {
  streampos size;
  char * memblock;
  std::vector<path> paths = std::vector<path>(directory_iterator(filename),boost::filesystem::directory_iterator());
  std::sort(paths.begin(), paths.end());
  bool already_have_resolution = false;
  int resolution = 0;
  double max_time = std::numeric_limits<double>::lowest();
  double sum_time = 0.0;

  for (uint path_idx = 0; path_idx < paths.size(); path_idx++) {
    path& entry = paths[path_idx];
    std::ifstream file (entry.string(), ios::in|ios::binary|ios::ate);

    timespec begin, end;
    clock_gettime(CLOCK_REALTIME, &begin);

    if (file.is_open())
    {
      size = file.tellg();
      memblock = new char [size];
      file.seekg (0, ios::beg);
      file.read (memblock, size);
      file.close();
      if (!already_have_resolution) {
        resolution = readResolution(memblock[0], memblock[1]);
        std::cout << "Resolution: " << resolution  << endl;
        already_have_resolution = true;
      }

      std::vector<Color> voxels;
      std::vector<int> count;

      bitset<8> bits;

      for (int i = 2; i < size; i++) {
        bits = bitset<8>(memblock[i]);
        if (bits[7] == 0 && bits[6] == 0) voxels.push_back(Color::W);
        else if (bits[7] == 0 && bits[6] == 1) voxels.push_back(Color::B);
        else if (bits[7] == 1 && bits[6] == 0) voxels.push_back(Color::G);
        else assert(false);

        bits.set(7,0);
        bits.set(6,0);

        count.push_back(bits.to_ulong()+1);
      }
      clock_gettime(CLOCK_REALTIME, &end);
      double time = end.tv_sec - begin.tv_sec + ((end.tv_nsec - begin.tv_nsec) / 1E9);
      sum_time += time;
      if (max_time < time) max_time = time; 
      std::cout << "Decoded " << entry.stem().string() << " -> "
                << time << " s." << std::endl;

      int sum = 0;
      for (int cou : count) {
        sum += cou;
      }
      // cout << "voxels size:" << sum << endl;
      assert(sum == resolution*resolution);

      if (write_slices) {
        write_slice_PNG_RLE(resolution, voxels, count, entry);
      }
      delete[] memblock;
    } else cout << "Unable to open file";
  }
  std::cout << "Max time of slicing: " << max_time << " s." << endl;
  std::cout << "Total time of slicing: " << sum_time << " s." << endl;
}

void readHEC_RLE_N_16(char* filename, bool write_slices) {
  streampos size;
  char * memblock;
  std::vector<path> paths = std::vector<path>(directory_iterator(filename),boost::filesystem::directory_iterator());
  std::sort(paths.begin(), paths.end());
  bool already_have_resolution = false;
  int resolution = 0;
  double max_time = std::numeric_limits<double>::lowest();
  double sum_time = 0.0;

  for (uint path_idx = 0; path_idx < paths.size(); path_idx++) {
    path& entry = paths[path_idx];
    std::ifstream file (entry.string(), ios::in|ios::binary|ios::ate);

    timespec begin, end;
    clock_gettime(CLOCK_REALTIME, &begin);

    if (file.is_open())
    {
      size = file.tellg();
      memblock = new char [size];
      file.seekg (0, ios::beg);
      file.read (memblock, size);
      file.close();
      if (!already_have_resolution) {
        resolution = readResolution(memblock[0], memblock[1]);
        std::cout << "Resolution: " << resolution  << endl;
        already_have_resolution = true;
      }
      
      std::vector<Color> voxels;
      std::vector<int> count;
      bitset<8> bits_first;
      bitset<8> bits_second;
      bitset<16> bits;
      for (int i = 2; i < size; i+=2) {
        bits_first = bitset<8>(memblock[i]);
        bits_second = bitset<8>(memblock[i+1]);

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

      clock_gettime(CLOCK_REALTIME, &end);
      double time = end.tv_sec - begin.tv_sec + ((end.tv_nsec - begin.tv_nsec) / 1E9);
      sum_time += time;
      if (max_time < time) max_time = time; 
      std::cout << "Decoded " << entry.stem().string() << " -> "
                << time << " s." << std::endl;
      
      int sum = 0;
      for (int cou : count) {
        sum += cou;
      }
      // cout << "voxels size:" << sum << endl;
      assert(sum == resolution*resolution);

      if (write_slices) {
        write_slice_PNG_RLE(resolution, voxels, count, entry);
      }
      
      delete[] memblock;
    } else cout << "Unable to open file";
  }
  std::cout << "Max time of slicing: " << max_time << " s." << endl;
  std::cout << "Total time of slicing: " << sum_time << " s." << endl;
}

void readHEC_RLE_A_8(char* filename, bool write_slices) {
  streampos size;
  char * memblock;
  std::vector<path> paths = std::vector<path>(directory_iterator(filename),boost::filesystem::directory_iterator());
  std::sort(paths.begin(), paths.end());
  bool already_have_resolution = false;
  int resolution = 0;
  double max_time = std::numeric_limits<double>::lowest();
  double sum_time = 0.0;

  for (uint path_idx = 0; path_idx < paths.size(); path_idx++) {
    path& entry = paths[path_idx];
    std::ifstream file (entry.string(), ios::in|ios::binary|ios::ate);

    timespec begin, end;
    clock_gettime(CLOCK_REALTIME, &begin);
    if (file.is_open())
    {
      size = file.tellg();
      memblock = new char [size];
      file.seekg (0, ios::beg);
      file.read (memblock, size);
      file.close();
      if (!already_have_resolution) {
        resolution = readResolution(memblock[0], memblock[1]);
        std::cout << "Resolution: " << resolution  << endl;
        already_have_resolution = true;
      }
      
      std::vector<Color> voxels;
      std::vector<int> count;

      Color order[2] = {W, B};
      Color color = W;
      int n_times_changed = 0;
      int n_gray = 0;
      bitset<8> bits;
      for (int i = 2; i < size; i++) {
        bits = bitset<8>(memblock[i]);
        if (bits.all()) {
          n_times_changed++;
          color = order[n_times_changed % 2];
          continue;
        } else if (bits.to_ulong() == 254) {
          count.push_back(254);
          voxels.push_back(color);
        } else {
          if (bits.any()) {
            count.push_back(bits.to_ulong());
            voxels.push_back(color);
          }
          if (i+1 < size) {
            // don't insert gray at the end
            count.push_back(1);
            voxels.push_back(Color::G);
            n_gray++;
          }
          n_times_changed++;
          color = order[n_times_changed % 2];
        }
      }

      clock_gettime(CLOCK_REALTIME, &end);
      double time = end.tv_sec - begin.tv_sec + ((end.tv_nsec - begin.tv_nsec) / 1E9);
      sum_time += time;
      if (max_time < time) max_time = time; 
      std::cout << "Decoded " << entry.stem().string() << " -> "
                << time << " s." << std::endl;

      int sum = 0;
      for (int cou : count) {
        sum += cou;
      }
      // cout << "voxels size:" << sum << endl;
      // cout << "difference: " << sum - resolution*resolution*resolution << endl;
      // cout << "total grays: " << n_gray << endl;

      assert(sum == resolution*resolution);

      if (write_slices) {
        write_slice_PNG_RLE(resolution, voxels, count, entry);
      }
      delete[] memblock;
    } else cout << "Unable to open file";
  }
  std::cout << "Max time of slicing: " << max_time << " s." << endl;
  std::cout << "Total time of slicing: " << sum_time << " s." << endl;
}

void readHEC_RLE_A_16(char* filename, bool write_slices) {
  streampos size;
  char * memblock;
  std::vector<path> paths = std::vector<path>(directory_iterator(filename),boost::filesystem::directory_iterator());
  std::sort(paths.begin(), paths.end());
  bool already_have_resolution = false;
  int resolution = 0;
  double max_time = std::numeric_limits<double>::lowest();
  double sum_time = 0.0;

  for (uint path_idx = 0; path_idx < paths.size(); path_idx++) {
    path& entry = paths[path_idx];
    std::ifstream file (entry.string(), ios::in|ios::binary|ios::ate);
    
    timespec begin, end;
    clock_gettime(CLOCK_REALTIME, &begin);

    if (file.is_open())
    {
      size = file.tellg();
      memblock = new char [size];
      file.seekg (0, ios::beg);
      file.read (memblock, size);
      file.close();
      if (!already_have_resolution) {
        resolution = readResolution(memblock[0], memblock[1]);
        std::cout << "Resolution: " << resolution  << endl;
        already_have_resolution = true;
      }
      
      std::vector<Color> voxels;
      std::vector<int> count;

      Color order[2] = {W, B};
      Color color = W;
      int n_times_changed = 0;
      int n_gray = 0;
      bitset<16> bits;
      bitset<8> bits_first;
      bitset<8> bits_second;
      for (int i = 2; i < size; i+=2) {
        bits_first = bitset<8>(memblock[i]);
        bits_second = bitset<8>(memblock[i+1]);

        for (uint idx = 0; idx < 8; idx++) {
          bits[idx] = bits_first[idx];
          bits[idx+8] = bits_second[idx];
        }

        if (bits.all()) {
          n_times_changed++;
          color = order[n_times_changed % 2];
          continue;
        } else if (bits.to_ulong() == 65534) {
          count.push_back(65534);
          voxels.push_back(color);
        } else {
          if (bits.any()) {
            count.push_back(bits.to_ulong());
            voxels.push_back(color);
          }
          if (i+2 < size ) {
            // don't insert gray at the end
            count.push_back(1);
            voxels.push_back(Color::G);
            n_gray++;
          }
          n_times_changed++;
          color = order[n_times_changed % 2];
        }

      }

      clock_gettime(CLOCK_REALTIME, &end);
      double time = end.tv_sec - begin.tv_sec + ((end.tv_nsec - begin.tv_nsec) / 1E9);
      sum_time += time;
      if (max_time < time) max_time = time; 
      std::cout << "Decoded " << entry.stem().string() << " -> "
                << time << " s." << std::endl;

      int sum = 0;
      for (int cou : count) {
        sum += cou;
      }
      // cout << "voxels size:" << sum << endl;
      // cout << "difference: " << sum - resolution*resolution*resolution << endl;
      // cout << "total grays: " << n_gray << endl;

      assert(sum == resolution*resolution);

      if (write_slices) {
        write_slice_PNG_RLE(resolution, voxels, count, entry);
      }
      
      delete[] memblock;
    } else cout << "Unable to open file";
  }
  std::cout << "Max time of slicing: " << max_time << " s." << endl;
  std::cout << "Total time of slicing: " << sum_time << " s." << endl;
}

void readHEC_mod(char* filename, bool write_slices) {
  streampos size;
  char * memblock;
  std::vector<path> paths = std::vector<path>(directory_iterator(filename),boost::filesystem::directory_iterator());
  std::sort(paths.begin(), paths.end());
  bool already_have_resolution = false;
  int resolution = 0;
  double max_time = std::numeric_limits<double>::lowest();
  double sum_time = 0.0;

  for (uint path_idx = 0; path_idx < paths.size(); path_idx++) {
    path& entry = paths[path_idx];

    std::ifstream file (entry.string(), ios::in|ios::binary|ios::ate);
    
    timespec begin, end;
    clock_gettime(CLOCK_REALTIME, &begin);

    if (file.is_open())
    {
      size = file.tellg();
      memblock = new char [size];
      file.seekg (0, ios::beg);
      file.read (memblock, size);
      file.close();
      if (!already_have_resolution) {
        resolution = readResolution(memblock[0], memblock[1]);
        std::cout << "Resolution: " << resolution  << endl;
        already_have_resolution = true;
      }
      Color curr_color = W;
      bool needs_to_set_color = true;
      std::vector<Color> voxels;
      bitset<8> bits;
      for (int i = 2; i < size; i++) {
        bits = bitset<8> (memblock[i]);
        if (needs_to_set_color) {
          if (bits[7] == 1) curr_color = G;
          else if (bits[6] == 1) curr_color = B;
          voxels.push_back(curr_color);
        }
        for (int idx = 0; idx < 8; idx+=2) {
          if (needs_to_set_color) {
            needs_to_set_color = false;
            continue;
          }
          if (bits[7-idx] == 0 && bits[6-idx] == 1) {
            if (curr_color == G) curr_color = B;
            else if (curr_color == B) curr_color = G;
            else {
              cout << "1 when W" << endl;
              assert(false);
            }
          } else if (bits[7-idx] == 1 && bits[6-idx] == 0) {
            if (curr_color == G) curr_color = W;
            else if (curr_color == W) curr_color = G;
            else {
              cout << "2 when B" << endl;
              assert(false);
            }
          }
           voxels.push_back(curr_color);
        }
      }

      clock_gettime(CLOCK_REALTIME, &end);
      double time = end.tv_sec - begin.tv_sec + ((end.tv_nsec - begin.tv_nsec) / 1E9);
      sum_time += time;
      if (max_time < time) max_time = time; 
      std::cout << "Decoded " << entry.stem().string() << " -> "
                << time << " s." << std::endl;

      // std::cout << "voxels size:" << voxels.size() << endl;
      assert((int) voxels.size() == resolution*resolution);

      if (write_slices) {
        write_slice_PNG(resolution, voxels, entry);
      }
      delete[] memblock;
    } else std::cout << "Unable to open file";
  }
  std::cout << "Max time of slicing: " << max_time << " s." << endl;
  std::cout << "Total time of slicing: " << sum_time << " s." << endl;
}

int main(int argc, char **argv) {
  
  if (argc != 4) {
    cout << "Usage: ./hec2png /path/to/hec/folder type_of_hec write_slices?" << endl; 
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
  bool write_slices = static_cast<bool>(std::atoi(argv[3]));
  timespec begin, end;
  clock_gettime(CLOCK_REALTIME, &begin);

  switch (type_of_hec) {
    case 0:
      readHEC_normal(argv[1], write_slices);
      break;
    case 1:
      readHEC_RLE_N_8(argv[1], write_slices);
      break;
    case 2:  
      readHEC_RLE_N_16(argv[1], write_slices);
      break;
    case 3:
      readHEC_RLE_A_8(argv[1], write_slices);
      break;
    case 4:
      readHEC_RLE_A_16(argv[1], write_slices);
      break;
    case 5:
      readHEC_mod(argv[1], write_slices);
      break;
    default:
      break;
  }

  clock_gettime(CLOCK_REALTIME, &end);
  std::cout << "Total time of decoding -> "
          << end.tv_sec - begin.tv_sec + ((end.tv_nsec - begin.tv_nsec) / 1E9) << " s." << std::endl;
  
  return 0;
}