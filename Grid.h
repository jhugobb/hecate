#ifndef GRID_H
#define GRID_H

#include <map>
#include <set>
#include "geometry/Geometry.h"
#include "geometry/2Dgrid.h"

struct ColoringConfiguration {
	bool useNaive, useBox;
	bool calculate_black_white;
	int grid_size;
	int selectedVoxelization;
	double threshold_raycasting;
  bool writePNG;
  bool writePLY;
  bool writeHEC;
  bool writeCSV;
  std::string filename;
};

class Grid {

  enum VoxelColor {
    WHITE,
    BLACK,
    GRAY
  };

  struct Voxel {
    VoxelColor color = VoxelColor::WHITE;
    double x;
    glm::vec3 normal;
    
    bool operator <(const Voxel other ) const {
      return x < other.x;
    }
  };

  public:
    Grid() {};
    Grid(unsigned int size, Geo::BBox space_);
    ~Grid() {}; // TODO
    
    /**
     * Voxelizes the model
     *
     * @param mesh Triangle Mesh of the model to be voxelized
     * @param qt 2Dgrid of the model that subdivided the space
     * @param useNaive boolean that dictates if the naive approach should be used
     * @param useBox
     * @param calculate_b_w boolean that dictates if black and white coloring should be calculated
     * @param threshold_raycasting
     * @param filename name of the file where to save the voxelization
     */ 
    void colorGrid(TriangleMesh* mesh, TwoDGrid* qt, ColoringConfiguration config, std::string filename);


    /**
     * Creates a PLY model with the voxelization
     *
     * @param filename The name of the file to be created
     */ 
    void writePLY(std::string filename);

  private:
    // Number of nodes in a dimension
    uint size_;
    
    // Size of a grid node in any direction
    double node_size;

    // Bounding Box of the model
    Geo::BBox space;

    // Triangles of the mesh_
    std::vector<Triangle*> triangles;

    // Vertices of the mesh_
    std::vector<glm::vec3> vertices;
    
    // Vector of voxels
    // TODO change to only have a single row of voxels at any given point in time
    std::vector<Voxel> elements;

    TriangleMesh* mesh_;

    std::vector<uint> white_runs;
    std::vector<uint> black_runs;
    std::vector<uint> gray_runs;
    
    uint current_run = 0;
    VoxelColor current_run_color = VoxelColor::WHITE;

    int curr_runs[6] = {0,0,0,0,0,0}; 
    VoxelColor curr_colors[6] = {VoxelColor::WHITE,
                                 VoxelColor::WHITE,
                                 VoxelColor::WHITE,
                                 VoxelColor::WHITE,
                                 VoxelColor::WHITE,
                                 VoxelColor::WHITE};
                                 
    bool needs_to_set_color[6] = {true,true,true,true,true,true};
    // used for alt
    int n_times_switched[2] = {0, 0}; 
    bool wroteAll1s[2] = {false, false};

    std::vector<Voxel> lastSlice;
    std::vector<double> similarPercents;

    bool testVoxelGray_Naive(int &representative, glm::vec3 voxel_min, glm::vec3 voxel_max, std::vector<int> candidates);
    
    void testVoxelGray_Box(int z, std::vector<Voxel> &voxels, std::vector<int> &candidates, glm::vec2 &row_coords);

    void writePLY(int x, int y, int z, Voxel& voxel, std::ofstream &out_fobj) const;

    void calculateBlackWhite(int z, 
                             glm::vec2 coords, 
                             node* quad_node, 
                             std::vector<Voxel> &voxels, 
                             double threshold);

    void saveSliceAsPNG(std::vector<Voxel> &voxels, uint y);

    void saveSliceAsHEC(std::vector<Voxel> &voxels, std::ofstream &bin_file);

    void saveSliceAsHEC_RLE_Naive_8b(std::vector<Voxel> &voxels, std::ofstream &bin_file, int y);
    
    void saveSliceAsHEC_RLE_Naive_16b(std::vector<Voxel> &voxels, std::ofstream &bin_file, int y);

    void saveSliceAsHEC_RLE_Alt_8b(std::vector<Voxel> &voxels, std::ofstream &bin_file, int y);
    
    void saveSliceAsHEC_RLE_Alt_16b(std::vector<Voxel> &voxels, std::ofstream &bin_file, int y);

    void calculateStatistics(std::vector<Voxel> &voxels, int y);

    void writeCSV(std::string filename);

};  
#endif