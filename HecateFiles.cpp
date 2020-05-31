#include "Grid.h"

#include <bitset>
#include <iostream>
#include <fstream>

void Grid::saveSliceAsHEC(std::vector<Voxel> &voxels, int y) {

  char* slice_memblock;
  std::ofstream bin_file(filename_ + "_" + std::to_string(size_) + 
                              "/Normal/"+ std::to_string(y) + ".hec", std::ios::binary | std::ios::out);
  bin_file.write((char*) &resolution_bits, sizeof(resolution_bits));

  if (bin_file.is_open()) {
    slice_memblock = new char[(size_*size_)/4];
    uint slice_idx = 0;
    for (uint i = 0; i < size_*size_; i+=4) {
      std::bitset<8> bits(0);

      if (voxels[i].color == VoxelColor::BLACK) {
        bits.set(1);
      } else if (voxels[i].color == VoxelColor::GRAY) {
        bits.set(0);
      }

      if (voxels[i+1].color == VoxelColor::BLACK) {
        bits.set(3);
      } else if (voxels[i+1].color == VoxelColor::GRAY) {
        bits.set(2);
      }

      if (voxels[i+2].color == VoxelColor::BLACK) {
        bits.set(5);
      } else if (voxels[i+2].color == VoxelColor::GRAY) {
        bits.set(4);
      }

      if (voxels[i+3].color == VoxelColor::BLACK) {
        bits.set(7);
      } else if (voxels[i+3].color == VoxelColor::GRAY) {
        bits.set(6);
      }

      slice_memblock[slice_idx++] =  static_cast<unsigned char>(bits.to_ulong());
    }

    bin_file.write(slice_memblock, size_*size_/4);

    delete[] slice_memblock;
    bin_file.close();
  } else {
    std::cout << "Unable to open " << filename_ + "_" + std::to_string(size_) + 
                              "/Normal/"+ std::to_string(y) + ".hec" << std::endl;
    assert(false);
  }
}

void Grid::saveSliceAsHEC_RLE_Naive_8b(std::vector<Voxel> &voxels, int y) {

  std::ofstream bin_file(filename_ + "_" + std::to_string(size_) + 
                              "/RLE_n_8/"+ std::to_string(y) + ".hec", std::ios::binary | std::ios::out);
  bin_file.write((char*) &resolution_bits, sizeof(resolution_bits));

  if (bin_file.is_open()) {
    char* slice_memblock;
    const int size_run = 64; // 2^6
    std::vector<std::bitset<8>> runs_to_write;
    std::bitset<8> bits_to_write;
    bool needs_to_set_color = true;
    int curr_runs = 0;
    VoxelColor curr_colors;

    for (uint z = 0; z < size_; z++) { 
      for (uint x = 0; x < size_; x++) {

        if (needs_to_set_color) {
          curr_colors = voxels[z*size_ + x].color;
          needs_to_set_color = false;
        }

        if (voxels[z*size_ + x].color == curr_colors) {
          curr_runs++;
          if (curr_runs >= size_run) {
            bits_to_write.set();
            switch (curr_colors) {
              case VoxelColor::BLACK:
                bits_to_write.set(7, 0);
                break;
              case VoxelColor::GRAY:
                bits_to_write.set(6, 0);
                break;
              default:
                bits_to_write.set(7, 0);
                bits_to_write.set(6, 0);
                break;
            }
            runs_to_write.push_back(bits_to_write);
            bits_to_write.reset();
            needs_to_set_color = true;
            curr_runs = 0;
          }
        } else {
          curr_runs--;
          bits_to_write = std::bitset<8>(curr_runs);
          switch (curr_colors) {
            case VoxelColor::BLACK:
              bits_to_write.set(6, 1);
              break;
            case VoxelColor::GRAY:
              bits_to_write.set(7, 1);
              break;
            default:
              break;
          }
          runs_to_write.push_back(bits_to_write);
          curr_colors = voxels[z*size_ + x].color;
          curr_runs = 1;
        }

      }
    }

    if (curr_runs > 0) {
      bits_to_write = std::bitset<8>(--curr_runs);
      switch (curr_colors) {
        case VoxelColor::BLACK:
          bits_to_write.set(6, 1);
          break;
        case VoxelColor::GRAY:
          bits_to_write.set(7, 1);
          break;
        default:
          break;
      }
      runs_to_write.push_back(bits_to_write);
    } 

    uint size_of_memblock = runs_to_write.size();
    slice_memblock = new char[size_of_memblock];
    uint i = 0;
    for (bitset<8> b : runs_to_write) {
      slice_memblock[i++] = static_cast<unsigned char>(b.to_ulong());
    }
    bin_file.write(slice_memblock, size_of_memblock);
    delete[] slice_memblock;
    bin_file.close();
  } else {
    std::cout << "Unable to open file." << std::endl;
    assert(false);
  }


}

void Grid::saveSliceAsHEC_RLE_Naive_16b(std::vector<Voxel> &voxels, int y) {
  
  std::ofstream bin_file(filename_ + "_" + std::to_string(size_) + 
                              "/RLE_n_16/"+ std::to_string(y) + ".hec", std::ios::binary | std::ios::out);
  bin_file.write((char*) &resolution_bits, sizeof(resolution_bits));
  
  if (bin_file.is_open()) {
    const int size_run = 16384; // 2^14
    std::vector<std::bitset<16>> runs_to_write;
    std::bitset<16> bits_to_write;
    bool needs_to_set_color = true;
    int curr_runs = 0;
    VoxelColor curr_colors;

    for (uint z = 0; z < size_; z++) { 
      for (uint x = 0; x < size_; x++) {

        if (needs_to_set_color) {
          curr_colors = voxels[z*size_ + x].color;
          needs_to_set_color = false;
        }

        if (voxels[z*size_ + x].color == curr_colors) {
          curr_runs++;
          if (curr_runs >= size_run) {
            bits_to_write.set();
            switch (curr_colors) {
              case VoxelColor::BLACK:
                bits_to_write.set(15, 0);
                break;
              case VoxelColor::GRAY:
                bits_to_write.set(14, 0);
                break;
              default:
                bits_to_write.set(15, 0);
                bits_to_write.set(14, 0);
                break;
            }
            runs_to_write.push_back(bits_to_write);
            bits_to_write.reset();
            needs_to_set_color = true;
            curr_runs = 0;
          }
        } else {
          curr_runs--;
          bits_to_write = std::bitset<16>(curr_runs);
          switch (curr_colors) {
            case VoxelColor::BLACK:
              bits_to_write.set(14, 1);
              break;
            case VoxelColor::GRAY:
              bits_to_write.set(15, 1);
              break;
            default:
              break;
          }
          runs_to_write.push_back(bits_to_write);
          curr_colors = voxels[z*size_ + x].color;
          curr_runs = 1;
        }

      }
    }

    if (curr_runs > 0) {
      bits_to_write = std::bitset<16>(--curr_runs);

      assert(bits_to_write[15] == 0 && bits_to_write[14] == 0);
      switch (curr_colors) {
        case VoxelColor::BLACK:
          bits_to_write.set(14, 1);
          break;
        case VoxelColor::GRAY:
          bits_to_write.set(15, 1);
          break;
        default:
          break;
      }
      runs_to_write.push_back(bits_to_write);
    } 

    uint size_of_memblock = runs_to_write.size();
    char16_t* memblock = new char16_t[size_of_memblock];
    int i = 0;
    for (bitset<16> b : runs_to_write) {
      memblock[i++] = static_cast<char16_t>(b.to_ulong());
    }
    bin_file.write((char*) memblock, size_of_memblock * sizeof(*memblock));
    delete[] memblock;
    bin_file.close();
  } else {
    std::cout << "Unable to open file." << std::endl;
    assert(false);
  }
}

void Grid::saveSliceAsHEC_RLE_Alt_8b(std::vector<Voxel> &voxels, int y) {

  std::ofstream bin_file(filename_ + "_" + std::to_string(size_) + 
                              "/RLE_a_8/"+ std::to_string(y) + ".hec", std::ios::binary | std::ios::out);
  bin_file.write((char*) &resolution_bits, sizeof(resolution_bits));

  if (bin_file.is_open()) {
    const int size_run = 254; // 2^8 - 1 - 1 (255 means jumping color)
    std::vector<std::bitset<8>> runs_to_write;
    std::bitset<8> bits_to_write;
    int curr_runs = 0;
    VoxelColor curr_colors = VoxelColor::WHITE;
    
    for (uint z = 0; z < size_; z++) { 
      for (uint x = 0; x < size_; x++) {

        // WORKS
        if (voxels[z*size_ + x].color == curr_colors) {
          curr_runs++;
          if (curr_runs >= size_run) {
            assert(curr_runs == 254);
            bits_to_write = std::bitset<8>(curr_runs);
            runs_to_write.push_back(bits_to_write);
            curr_runs = 0;
          }
        } else if (voxels[z*size_ + x].color == VoxelColor::GRAY){
          bits_to_write = std::bitset<8>(curr_runs);
          runs_to_write.push_back(bits_to_write);
          curr_colors = static_cast<VoxelColor>((curr_colors+1) % 2);
          curr_runs = 0;
        } else {
          if (curr_runs == 0) {
            // Run of unexpected color
            bits_to_write = std::bitset<8>(255);
            runs_to_write.push_back(bits_to_write);
            curr_colors = static_cast<VoxelColor>((curr_colors+1) % 2);
            assert(curr_colors == VoxelColor::BLACK || curr_colors == VoxelColor::WHITE);
            curr_runs++;
          } else {
            cout << "Curr color: " << curr_colors << endl;
            cout << "Actual color: " << voxels[z*size_ + x].color << endl;
            cout << "Run length: " << curr_runs << endl;
            cout << "Slice number: " << y << endl;
            cout << "Row number: " << z << endl;

            // Adjacent Black and white!
            assert(false);
          }
        }

      }
    }

    // if (y == (int) size_-1) {
      bits_to_write = std::bitset<8>(curr_runs);
      runs_to_write.push_back(bits_to_write);
    // } 

    uint size_of_memblock = runs_to_write.size();
    char* memblock = new char[size_of_memblock];
    int i = 0;
    for (bitset<8> b : runs_to_write) {
      memblock[i++] = static_cast<char>(b.to_ulong());
    }
    bin_file.write(memblock, size_of_memblock * sizeof(*memblock));
    delete[] memblock;
    bin_file.close();
  } else {
    std::cout << "Unable to open file." << std::endl;
    assert(false);
  }
}

void Grid::saveSliceAsHEC_RLE_Alt_16b(std::vector<Voxel> &voxels, int y) {
  
  std::ofstream bin_file(filename_ + "_" + std::to_string(size_) + 
                              "/RLE_a_16/"+ std::to_string(y) + ".hec", std::ios::binary | std::ios::out);
  bin_file.write((char*) &resolution_bits, sizeof(resolution_bits));
  
  if (bin_file.is_open()) {
    const int size_run = 65534; // 2^16 - 1 - 1 (65535 means jumping color)
    std::vector<std::bitset<16>> runs_to_write;
    std::bitset<16> bits_to_write;
    int curr_runs = 0;
    VoxelColor curr_colors = VoxelColor::WHITE;

    for (uint z = 0; z < size_; z++) { 
      for (uint x = 0; x < size_; x++) {
        
        // WORKS
        if (voxels[z*size_ + x].color == curr_colors) {
          curr_runs++;
          if (curr_runs >= size_run) {
            assert(curr_runs == 65534);
            bits_to_write = std::bitset<16>(curr_runs);
            runs_to_write.push_back(bits_to_write);
            curr_runs = 0;
          }
        } else if (voxels[z*size_ + x].color == VoxelColor::GRAY){
          bits_to_write = std::bitset<16>(curr_runs);
          runs_to_write.push_back(bits_to_write);
          curr_colors = static_cast<VoxelColor>((curr_colors+1) % 2);
          curr_runs = 0;
        } else {
          if (curr_runs == 0) {
            // Run of unexpected color
            bits_to_write = std::bitset<16>(65535);
            runs_to_write.push_back(bits_to_write);
            curr_colors = static_cast<VoxelColor>((curr_colors+1) % 2);
            assert(curr_colors == VoxelColor::BLACK || curr_colors == VoxelColor::WHITE);
            curr_runs++;
          } else {
            cout << "Curr color: " << curr_colors << endl;
            cout << "Actual color: " << voxels[z*size_ + x].color << endl;
            cout << "Run length: " << curr_runs << endl;
            
            // Adjacent Black and white!
            assert(false);
          }
        }
      }
    }

    // if (y == (int) size_-1) {
      bits_to_write = std::bitset<16>(curr_runs);
      runs_to_write.push_back(bits_to_write);
    // } 

    uint size_of_memblock = runs_to_write.size();
    char16_t* memblock = new char16_t[size_of_memblock];
    int i = 0;
    for (bitset<16> b : runs_to_write) {
      memblock[i++] = static_cast<char16_t>(b.to_ulong());
    }
    bin_file.write((char*) memblock, size_of_memblock * sizeof(*memblock));
    delete[] memblock;
    bin_file.close();
  } else {
    std::cout << "Unable to open file." << std::endl;
    assert(false);
  }
}

void Grid::saveSliceAsHEC_Mod_Enc(std::vector<Voxel> &voxels, int y) {

  std::ofstream bin_file(filename_ + "_" + std::to_string(size_) + 
                              "/Mod/"+ std::to_string(y) + ".hec", std::ios::binary | std::ios::out);
  bin_file.write((char*) &resolution_bits, sizeof(resolution_bits));

  if (bin_file.is_open()) {
    char block;
    VoxelColor curr_color = voxels[0].color;
    bitset<8> encoding;
    
    switch (curr_color) {
      case VoxelColor::BLACK:
        encoding.set(6);      
        break;
      case VoxelColor::GRAY:
        encoding.set(7);
        break;
      default:
        break;
    }
    uint check = 1;
    uint enc_idx = 2;
    uint val;
    bool flag = true;
    for (uint z = 0; z < size_; z++)
    for (uint x = 0; x < size_; x++) {
      if (flag) {
        flag = false;
        continue;
      }
      val = abs(voxels[z*size_ + x].color - curr_color);
      curr_color = voxels[z*size_ + x].color;
      
      if (val == 1) {
        encoding.set(7 - (enc_idx + 1));
      } else if (val == 2) {
        encoding.set(7 - enc_idx);
      }

      enc_idx+=2;
      check++;
      if (enc_idx >= 8) {
        block = static_cast<unsigned char>(encoding.to_ulong());
        bin_file.write(&block, sizeof(block));
        enc_idx = 0;
        encoding.reset();
      }
    }

    bin_file.close();
    // cout << "Check: " << check << " vs size: " << size_*size_ << endl; 
    assert(check == size_*size_);
  } else {
    std::cout << "Unable to open " << filename_ + "_" + std::to_string(size_) + 
                              "/Mod/"+ std::to_string(y) + ".hec" << std::endl;
    assert(false);
  }
}

void Grid::saveSliceAsHEC_Mod_Slice(std::vector<Voxel> &voxels, int y) {
  char* slice_memblock;
  std::ofstream bin_file(filename_ + "_" + std::to_string(size_) + 
                              "/Mod_Slice/"+ std::to_string(y) + ".hec", std::ios::binary | std::ios::out);
  bin_file.write((char*) &resolution_bits, sizeof(resolution_bits));

  if (y == 0) {
    lastSlice = std::vector<Voxel>(voxels);
  
    if (bin_file.is_open()) {
      slice_memblock = new char[(size_*size_)/4];
      uint slice_idx = 0;
      for (uint i = 0; i < size_*size_; i+=4) {
        std::bitset<8> bits(0);

        if (voxels[i].color == VoxelColor::BLACK) {
          bits.set(1);
        } else if (voxels[i].color == VoxelColor::GRAY) {
          bits.set(0);
        }

        if (voxels[i+1].color == VoxelColor::BLACK) {
          bits.set(3);
        } else if (voxels[i+1].color == VoxelColor::GRAY) {
          bits.set(2);
        }

        if (voxels[i+2].color == VoxelColor::BLACK) {
          bits.set(5);
        } else if (voxels[i+2].color == VoxelColor::GRAY) {
          bits.set(4);
        }

        if (voxels[i+3].color == VoxelColor::BLACK) {
          bits.set(7);
        } else if (voxels[i+3].color == VoxelColor::GRAY) {
          bits.set(6);
        }

        slice_memblock[slice_idx++] =  static_cast<unsigned char>(bits.to_ulong());
      }

      bin_file.write(slice_memblock, size_*size_/4);

      delete[] slice_memblock;
      bin_file.close();
    } else {
      std::cout << "Unable to open " << filename_ + "_" + std::to_string(size_) + 
                                "/Mod_Slice/"+ std::to_string(y) + ".hec" << std::endl;
      assert(false);
    }
  } else {
    if (bin_file.is_open()) {
      slice_memblock = new char[(size_*size_)/4];
      uint slice_idx = 0;
      for (uint i = 0; i < size_*size_; i+=4) {
        std::bitset<8> bits(0);

        switch (abs(static_cast<int>(lastSlice[i].color)-static_cast<int>(voxels[i].color))) {
          case 1:
            bits.set(1);
            break;
          case 2:
            bits.set(0);
            break;
        }

        switch (abs(static_cast<int>(lastSlice[i+1].color)-static_cast<int>(voxels[i+1].color))) {
          case 1:
            bits.set(3);
            break;
          case 2:
            bits.set(2);
            break;
        }

        switch (abs(static_cast<int>(lastSlice[i+2].color)-static_cast<int>(voxels[i+2].color))) {
          case 1:
            bits.set(5);
            break;
          case 2:
            bits.set(4);
            break;
        }

        switch (abs(static_cast<int>(lastSlice[i+3].color)-static_cast<int>(voxels[i+3].color))) {
          case 1:
            bits.set(7);
            break;
          case 2:
            bits.set(6);
            break;
        }

        slice_memblock[slice_idx++] =  static_cast<unsigned char>(bits.to_ulong());
      }
      lastSlice = voxels;
      bin_file.write(slice_memblock, size_*size_/4);

      delete[] slice_memblock;
      bin_file.close();
    } else {
      std::cout << "Unable to open " << filename_ + "_" + std::to_string(size_) + 
                                "/Mod_Slice/"+ std::to_string(y) + ".hec" << std::endl;
      assert(false);
    }
  }
}

void Grid::saveSliceAsHEC_Mod_RLE(std::vector<Voxel> &voxels, int y) {
  
  std::vector<VoxelColor> mod_slice(voxels.size());
  if (y == 0) {
    lastSlice_rle = voxels;
    for (uint i = 0; i < voxels.size(); i++) {
      mod_slice[i] = voxels[i].color;
    }
    saveColors_RLE_N_16b(mod_slice, y);
  } else {
      for (uint i = 0; i < voxels.size(); i++) {
        switch (abs(voxels[i].color - lastSlice_rle[i].color))
        {
        case 0:
          mod_slice[i] = WHITE;
          break;
        case 1: 
          mod_slice[i] = BLACK;
          break;
        case 2:
          mod_slice[i] = GRAY;
          break;
        default:
          assert(false);
          break;
        }
      } 
      lastSlice_rle = voxels;
      saveColors_RLE_N_16b(mod_slice, y);
  }
}

void Grid::saveColors_RLE_N_16b(std::vector<VoxelColor> &voxels, int y) {
  
  std::ofstream bin_file(filename_ + "_" + std::to_string(size_) + 
                              "/Mod_RLE/"+ std::to_string(y) + ".hec", std::ios::binary | std::ios::out);
  bin_file.write((char*) &resolution_bits, sizeof(resolution_bits));
  
  if (bin_file.is_open()) {
    const int size_run = 16384; // 2^14
    std::vector<std::bitset<16>> runs_to_write;
    std::bitset<16> bits_to_write;
    bool needs_to_set_color = true;
    int curr_runs = 0;
    VoxelColor curr_colors;

    for (uint z = 0; z < size_; z++) { 
      for (uint x = 0; x < size_; x++) {

        if (needs_to_set_color) {
          curr_colors = voxels[z*size_ + x];
          needs_to_set_color = false;
        }

        if (voxels[z*size_ + x] == curr_colors) {
          curr_runs++;
          if (curr_runs >= size_run) {
            bits_to_write.set();
            switch (curr_colors) {
              case VoxelColor::BLACK:
                bits_to_write.set(15, 0);
                break;
              case VoxelColor::GRAY:
                bits_to_write.set(14, 0);
                break;
              default:
                bits_to_write.set(15, 0);
                bits_to_write.set(14, 0);
                break;
            }
            runs_to_write.push_back(bits_to_write);
            bits_to_write.reset();
            needs_to_set_color = true;
            curr_runs = 0;
          }
        } else {
          curr_runs--;
          bits_to_write = std::bitset<16>(curr_runs);
          switch (curr_colors) {
            case VoxelColor::BLACK:
              bits_to_write.set(14, 1);
              break;
            case VoxelColor::GRAY:
              bits_to_write.set(15, 1);
              break;
            default:
              break;
          }
          runs_to_write.push_back(bits_to_write);
          curr_colors = voxels[z*size_ + x];
          curr_runs = 1;
        }

      }
    }

    if (curr_runs > 0) {
      bits_to_write = std::bitset<16>(--curr_runs);

      assert(bits_to_write[15] == 0 && bits_to_write[14] == 0);
      switch (curr_colors) {
        case VoxelColor::BLACK:
          bits_to_write.set(14, 1);
          break;
        case VoxelColor::GRAY:
          bits_to_write.set(15, 1);
          break;
        default:
          break;
      }
      runs_to_write.push_back(bits_to_write);
    } 

    uint size_of_memblock = runs_to_write.size();
    char16_t* memblock = new char16_t[size_of_memblock];
    int i = 0;
    for (bitset<16> b : runs_to_write) {
      memblock[i++] = static_cast<char16_t>(b.to_ulong());
    }
    bin_file.write((char*) memblock, size_of_memblock * sizeof(*memblock));
    delete[] memblock;
    bin_file.close();
  } else {
    std::cout << "Unable to open file." << std::endl;
    assert(false);
  }
}
