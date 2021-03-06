#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <chrono>
#include <time.h>
#define GLM_FORCE_RADIANS
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/matrix_inverse.hpp>
#include "imgui/imgui.h"

#include "Scene.h"

Scene::Scene()
{
	mesh = NULL;
}

Scene::~Scene()
{
	if(mesh != NULL)
	{
		mesh->free();
		delete mesh;
	}
}


void Scene::init()
{
	initShaders();

	currentTime = 0.0f;
	
	camera.init(2.0f);
	
	bLighting = true;
	bWireframe = false;
	config.selectedVoxelization = 1;
	config.grid_size = 512;
	config.writePLY = true;
	config.writePNG = true;
	config.filename = "test";
	// lambda = 1.f;
	// selectedIterativeFunction = 0;
	// selectedGlobalFunction = 0;
	// percentConstraints = 100;
	// constraintWeight = 1.0f;
}

bool Scene::loadScan(const char *filename)
{
	ifstream fin;
	
	// Check if file exists
	fin.open(filename);
	if(!fin.is_open())
		return false;
	
	if(mesh != NULL)
		delete mesh;
	mesh = new TriangleMesh();
	if(mesh->load(filename))
	{
		cout << "Triangle mesh " << filename << " loaded. ";
		cout << mesh->getVertices().size() << " vertices. ";
		cout << (mesh->getTriangles().size()) << " faces." << endl;
	}
	else
		return false;

	// Compute global bounding box
	glm::vec3 bbox[2], center, size;
	
	getBBox(bbox);
	center = (bbox[0] + bbox[1]) / 2.0f;
	size = bbox[1] - bbox[0];
	model_bbox.addPoint(bbox[0]);
	model_bbox.addPoint(bbox[1]);
	size = glm::vec3(glm::max(size.x, glm::max(size.y, size.z)));
	std::cout << "min: " << bbox[0].x << " " << bbox[0].y << " " << bbox[0].z << std::endl;	
	std::cout << "max: " << bbox[1].x << " " << bbox[1].y << " " << bbox[1].z << std::endl;	
	mat = glm::mat4(1.0f);
	mat = glm::scale(mat, glm::vec3(1.0f - 0.1f, 1.0f - 0.1f, 1.0f - 0.1f) / size);
	mat = glm::translate(mat, -center);
	// transform(matrix);
	
	// Save original vertices for reload
	originalVertices.clear();
	copy(mesh->getVertices().begin(), mesh->getVertices().end(), std::back_inserter(originalVertices));
	
	// for (uint i = 0; i < mesh->getVertices().size(); i++) {
	// 	model_bbox.addPoint(mesh->getVertices()[i]);
	// }
	return true;
}

void Scene::update(int deltaTime)
{
	currentTime += deltaTime;
}

void Scene::render()
{
	pointsProgram.use();
	pointsProgram.setUniform1i("bLighting", int(bLighting));
	pointsProgram.setUniform1i("bBlack", false);
	pointsProgram.setUniformMatrix4f("projection", mat *camera.getProjectionMatrix());
	pointsProgram.setUniformMatrix4f("modelview", mat * camera.getModelViewMatrix());
	pointsProgram.setUniformMatrix3f("normalMatrix", glm::inverseTranspose(glm::mat3(mat *camera.getModelViewMatrix())));
	
	if(bWireframe)
	{
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
	}
	if(mesh != NULL)
		mesh->render();
	if(bWireframe)
	{
		glDisable(GL_POLYGON_OFFSET_FILL);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		pointsProgram.setUniform1i("bBlack", true);
		mesh->render();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
}

void Scene::render_gui()
{
	ImGui::Begin("Options", 0, ImGuiWindowFlags_AlwaysAutoResize);
	
	// Model options
	ImGui::Text("Model");
	ImGui::Spacing();
	if(ImGui::Button("Reload"))
		reloadMesh();
	ImGui::Spacing();

	// Render options
	ImGui::Text("Render");
	ImGui::Spacing();
	ImGui::Checkbox("Lighting", &bLighting);
	ImGui::Spacing();
	ImGui::Checkbox("Wireframe", &bWireframe);
	ImGui::Spacing();
	ImGui::Separator();

	// Voxelization Options
	ImGui::Text("Voxelization Options");
	ImGui::Spacing();
	ImGui::RadioButton("Naive Approach", &config.selectedVoxelization, 0); ImGui::SameLine();
	ImGui::RadioButton("Box Approach", &config.selectedVoxelization, 1); ImGui::SameLine();
	ImGui::RadioButton("Bin Tree Approach", &config.selectedVoxelization, 2);
	ImGui::Spacing();
	ImGui::Separator();
	ImGui::Checkbox("Calculate Black and White?", &config.calculate_black_white);
	ImGui::Spacing();
	ImGui::Text("Threshold of Ray Validity: ");
	ImGui::SameLine();
	ImGui::SetNextItemWidth(110);
	ImGui::InputDouble("##Threshold of Ray Validity", &config.threshold_raycasting, 0.0, 0.05);
	config.threshold_raycasting = glm::max(0.0, glm::min(config.threshold_raycasting, 3.0));
	ImGui::Spacing();

	//Output options
	ImGui::Separator();
	ImGui::Text("Voxelization Options");
	ImGui::Checkbox("Write slices as PNGs?", &config.writePNG);
	ImGui::Spacing();
	ImGui::Checkbox("Write grays as PLY?", &config.writePLY);
	ImGui::Spacing();
	ImGui::Checkbox("Write Hecate File?", &config.writeHEC);
	ImGui::Spacing();
	ImGui::Checkbox("Calculate Statistics as CSV?", &config.writeCSV);
	ImGui::Spacing();
	ImGui::Text("Prefix of files written: ");
	ImGui::SameLine();
	ImGui::SetNextItemWidth(110);
	ImGui::InputText("##Filename", &config.filename.at(0), 64);
	ImGui::Separator();

	// Grid parameters
	ImGui::Text("Grid Parameters");
	ImGui::Spacing();

	ImGui::Text("Grid size: ");
	ImGui::SameLine();
	ImGui::SetNextItemWidth(110);
	ImGui::InputInt("##Gridsize", &config.grid_size, 50, 1);
	config.grid_size = glm::max(0, glm::min(config.grid_size, 99999));
	ImGui::Spacing();
	
	if (ImGui::Button("Voxelize")) {

		std::cout << "=================== Voxelization ===================" << std::endl;	

		timespec begin, end_grid, begin_bt, end_bt, begin_vox, end_vox, end;

		/* =================== 2D Grid =================== */
		std::cout << "Starting Two D Grid Generation" << std::endl;	
		clock_gettime(CLOCK_REALTIME, &begin);

		twodgrid = new TwoDGrid(mesh, config.grid_size, model_bbox);
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
				twodgrid->saveMinXs();
				clock_gettime(CLOCK_REALTIME, &end_bt);
				std::cout << "Finished saving min xs in " << end_bt.tv_sec - begin_bt.tv_sec + ((end_bt.tv_nsec - begin.tv_nsec) / 1E9)<< " s." << std::endl;
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
	}

	ImGui::End();
}

Camera &Scene::getCamera()
{
	return camera;
}

void Scene::initShaders()
{
	Shader vShader, fShader;

	// Load shader for point rendering
	vShader.initFromFile(VERTEX_SHADER, "shaders/points.vert");
	if(!vShader.isCompiled())
	{
		cout << "Vertex Shader Error" << endl;
		cout << "" << vShader.log() << endl << endl;
	}
	fShader.initFromFile(FRAGMENT_SHADER, "shaders/points.frag");
	if(!fShader.isCompiled())
	{
		cout << "Fragment Shader Error" << endl;
		cout << "" << fShader.log() << endl << endl;
	}
	pointsProgram.init();
	pointsProgram.addShader(vShader);
	pointsProgram.addShader(fShader);
	pointsProgram.link();
	if(!pointsProgram.isLinked())
	{
		cout << "Shader Linking Error" << endl;
		cout << "" << pointsProgram.log() << endl << endl;
	}
	pointsProgram.bindFragmentOutput("outColor");
	vShader.free();
	fShader.free();
}

void Scene::getBBox(glm::vec3 bbox[2]) const
{
	mesh->getBBox(bbox);
}

void Scene::transform(const glm::mat4 &matrix)
{
	mesh->transform(matrix);
	mesh->sendToOpenGL(pointsProgram);
}

void Scene::reloadMesh()
{
	mesh->getVertices().clear();
	copy(originalVertices.begin(), originalVertices.end(), std::back_inserter(mesh->getVertices()));
	mesh->sendToOpenGL(pointsProgram);
}











