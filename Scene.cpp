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
	config.n_rays_per_voxel = 4;
	config.selectedVoxelization = 1;
	config.grid_size = 64;
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
	
	for (uint i = 0; i < mesh->getVertices().size(); i++) {
		model_bbox.addPoint(mesh->getVertices()[i]);
	}
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
	ImGui::Text("Number of Rays per Voxel Row: ");
	ImGui::SameLine();
	ImGui::SetNextItemWidth(110);
	ImGui::InputInt("##Number of Rays Per Voxel Row", &config.n_rays_per_voxel, 4, 1);
	config.n_rays_per_voxel = glm::max(0, glm::min(config.n_rays_per_voxel, 64));
	ImGui::Separator();
	ImGui::Spacing();

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
		std::cout << "Finished Two D Grid in " << end_grid.tv_sec - begin.tv_sec << " s." << std::endl;
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
				std::cout << "Finished Bin Tree Creation in " << end_bt.tv_sec - begin_bt.tv_sec << " s." << std::endl;
				config.useNaive = false;
				config.useBox = false;
				break;
		}
		Grid grid(config.grid_size, model_bbox);

		// cout << "Model bbox min: " << model_bbox.minPoint.x << " " << model_bbox.minPoint.y << " " << model_bbox.minPoint.z << endl; 
		// cout << "Model bbox max: " << model_bbox.maxPoint.x << " " << model_bbox.maxPoint.y << " " << model_bbox.maxPoint.z << endl; 

		// VOXELIZATION
		clock_gettime(CLOCK_REALTIME, &begin_vox);
		
		grid.colorGrid(mesh, twodgrid, config, "test.ply");

		clock_gettime(CLOCK_REALTIME, &end_vox);

		std::cout << "Finished Voxelization in " << end_vox.tv_sec - begin_vox.tv_sec << " s." << std::endl;
		delete twodgrid;
		/* =================== Voxelization =================== */

		clock_gettime(CLOCK_REALTIME, &end);
		
		std::cout << "Finished execution in " << end.tv_sec - begin.tv_sec << " s." << std::endl;
	}

	// ImGui::Text("#Iterations: ");
	// ImGui::SameLine();
	// ImGui::SetNextItemWidth(75);
	// ImGui::InputInt("##Iterations", &nIterations, 1, 10);
	// nIterations = glm::max(1, glm::min(nIterations, 100));
	// ImGui::Spacing();
	// ImGui::RadioButton("Laplacian", &selectedIterativeFunction, 0); ImGui::SameLine();
	// ImGui::RadioButton("Bilaplacian", &selectedIterativeFunction, 1); ImGui::SameLine();
	// ImGui::RadioButton("Lambda-Nu", &selectedIterativeFunction, 2);
	// ImGui::Spacing();
	// if(ImGui::Button("Apply"))
	// {
	// 	LaplacianSmoothing smoothing;
		
	// 	smoothing.setMesh(mesh);
	// 	switch(selectedIterativeFunction)
	// 	{
	// 	case 0:
	// 		smoothing.iterativeLaplacian(nIterations, lambda);
	// 		break;
	// 	case 1:
	// 		smoothing.iterativeBilaplacian(nIterations, lambda);
	// 		break;
	// 	case 2:
	// 		smoothing.iterativeLambdaNu(nIterations, lambda);
	// 		break;
	// 	}
	// 	mesh->sendToOpenGL(pointsProgram);
	// }
	// ImGui::Spacing();
	// ImGui::Separator();
	// ImGui::Spacing();

	// // Global smoothing
	// ImGui::Text("Global Smoothing");
	// ImGui::Spacing();

	// ImGui::Text("Constraints: ");
	// ImGui::SameLine();
	// ImGui::SetNextItemWidth(110);
	// ImGui::InputInt("##Percentage", &percentConstraints, 1, 10);
	// percentConstraints = glm::max(1, glm::min(percentConstraints, 100));
	// ImGui::Spacing();

	// ImGui::Text("Weight: ");
	// ImGui::SameLine();
	// ImGui::SetNextItemWidth(110);
	// ImGui::InputFloat("##Weight", &constraintWeight, 0.01f, 0.1f, "%.2f");
	// constraintWeight = glm::max(0.f, glm::min(constraintWeight, 1.f));
	// ImGui::Spacing();

	// ImGui::RadioButton("Laplacian##2", &selectedGlobalFunction, 0); ImGui::SameLine();
	// ImGui::RadioButton("Bilaplacian##2", &selectedGlobalFunction, 1);
	// ImGui::Spacing();
	// if(ImGui::Button("Apply##2"))
	// {
	// 	vector<unsigned int> selected;
	// 	vector<bool> constraints;
	// 	LaplacianSmoothing smoothing;
		
	// 	srand(time(NULL));
	// 	selected.resize(mesh->getVertices().size());
	// 	for(unsigned int i=0; i<mesh->getVertices().size(); i++)
	// 		selected[i] = i;
	// 	random_shuffle(selected.begin(), selected.end());
	// 	selected.resize(selected.size() * percentConstraints / 100);
	// 	constraints.resize(mesh->getVertices().size(), false);
	// 	for(unsigned int i=0; i<selected.size(); i++)
	// 		constraints[selected[i]] = true;
	// 	smoothing.setMesh(mesh);
	// 	switch(selectedGlobalFunction)
	// 	{
	// 	case 0:
	// 		smoothing.globalLaplacian(constraints);
	// 		break;
	// 	case 1:
	// 		smoothing.globalBilaplacian(constraints, constraintWeight);
	// 		break;
	// 	}
	// 	mesh->sendToOpenGL(pointsProgram);
	// }

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











