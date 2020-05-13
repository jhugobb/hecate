#include <iostream>
#include <string>
#include <fstream>

#include "../TriangleMesh.h"
#include "../Grid.h"


//Remove console (only works in Visual Studio)
#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

int main(int argc, char **argv)
{
	// Args
	// ./bash_hecate filename grid_size #voxelization_technique 
	// calculateblackwhite? thresholdraycasting writePNG? writePLY? writeCSV?
	
	if (argc != 10) {
		cout << "Usage: ./hecate_bash filename grid_size #voxelization_technique calculateblackwhite? thresholdraycasting writePNG? writePLY? writeCSV? writeHEC?\n";
		return -1;
	} 

	//Filename needs to be models/path/to/model
	std::string filename(argv[1]);
	
	ColoringConfiguration config;

	config.grid_size = std::stoi(argv[2]);

	config.selectedVoxelization = std::stoi(argv[3]);

	config.calculate_black_white = std::stoi(argv[4]) == 1;
	
	config.threshold_raycasting = std::stod(argv[5]);

	config.writePNG = std::stoi(argv[6]) == 1;
	config.writePLY = std::stoi(argv[7]) == 1;
	config.writeCSV = std::stoi(argv[8]) == 1;
	config.writeHEC = std::stoi(argv[9]) == 1;

	// Make sure that output/ exists
	config.filename = "output/" + filename.substr(7,filename.size()-4-7) + "_new";
	std::cout << "filename " <<  config.filename << std::endl;

	ifstream fin;
	
	// Check if file exists
	fin.open(filename);
	if(!fin.is_open())
		return false;
	
	TriangleMesh* mesh = NULL;
	mesh = new TriangleMesh();

	if(mesh->load(filename)) {
		cout << "Triangle mesh " << filename << " loaded. ";
		cout << mesh->getVertices().size() << " vertices. ";
		cout << (mesh->getTriangles().size()) << " faces." << endl;
	}
	else {
		cout << "Could not load mesh." << std::endl;
		return -1;
	}
	Geo::BBox model_bbox;

	for (uint i = 0; i < mesh->getVertices().size(); i++) {
		model_bbox.addPoint(mesh->getVertices()[i]);
	}
	// Compute global bounding box
	glm::vec3 center, size;
	
	center = (model_bbox.minPoint + model_bbox.maxPoint) / 2.0f;
	size = model_bbox.maxPoint + model_bbox.minPoint;

	size = glm::vec3(glm::max(size.x, glm::max(size.y, size.z)));
	std::cout << "min: " << model_bbox.minPoint.x << " " << model_bbox.minPoint.y << " " << model_bbox.minPoint.z << std::endl;	
	std::cout << "max: " << model_bbox.maxPoint.x << " " << model_bbox.maxPoint.y << " " << model_bbox.maxPoint.z << std::endl;		

	std::cout << "=================== Voxelization ===================" << std::endl;	

	timespec begin, end_grid, begin_bt, end_bt, begin_vox, end_vox, end;

	/* =================== 2D Grid =================== */
	std::cout << "Starting Two D Grid Generation" << std::endl;	
	clock_gettime(CLOCK_REALTIME, &begin);

	TwoDGrid* twodgrid = new TwoDGrid(mesh, config.grid_size, model_bbox);
	#pragma omp parallel for
	for (uint i = 0; i < mesh->getTriangles().size(); i++) {
		twodgrid->insert(i);
	}

	clock_gettime(CLOCK_REALTIME, &end_grid);
	std::cout << "Finished Two D Grid in " << end_grid.tv_sec - begin.tv_sec + ((end_grid.tv_nsec - begin.tv_nsec) / 1E9) << " s." << std::endl;
	/* =================== 2D Grid =================== */
	

	/* =================== Voxelization =================== */
	switch (config.selectedVoxelization) {
		case 0: 
			// Naive 
			config.useNaive = true;
			config.useBox = false;
			break;
		case 1: 
			// Box
			config.useNaive = false;
			config.useBox = true;
			break;
		default:
			// Bin Tree
			clock_gettime(CLOCK_REALTIME, &begin_bt);
			twodgrid->buildBinTrees();
			clock_gettime(CLOCK_REALTIME, &end_bt);
			std::cout << "Finished Bin Tree Creation in " << end_bt.tv_sec - begin_bt.tv_sec + ((end_bt.tv_nsec - begin_bt.tv_nsec) / 1E9)<< " s." << std::endl;
			config.useNaive = false;
			config.useBox = false;
			break;
	}
	Grid grid(config.grid_size, model_bbox, mesh);

	// VOXELIZATION
	clock_gettime(CLOCK_REALTIME, &begin_vox);
	
	grid.colorGrid(twodgrid, config, config.filename);

	clock_gettime(CLOCK_REALTIME, &end_vox);

	std::cout << "Finished Voxelization in " << end_vox.tv_sec - begin_vox.tv_sec + ((end_vox.tv_nsec - begin_vox.tv_nsec) / 1E9)<< " s." << std::endl;
	delete twodgrid;
	/* =================== Voxelization =================== */

	clock_gettime(CLOCK_REALTIME, &end);
	
	std::cout << "Finished execution in " << end.tv_sec - begin.tv_sec + ((end.tv_nsec - begin.tv_nsec) / 1E9)<< " s." << std::endl;
	return 0;
}



