#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include "geometry/SATProjection.h"
#include "TriangleMesh.h"



TriangleMesh::TriangleMesh()
{
	bGLObjsInit = false;
}


void TriangleMesh::init(const vector<glm::vec3> &newVertices, const vector<Triangle> &newTriangles)
{
	copy(newVertices.begin(), newVertices.end(), back_inserter(vertices));
	copy(newTriangles.begin(), newTriangles.end(), back_inserter(triangles));
}

bool TriangleMesh::load(const string &filename)
{
	ifstream fin;
	char line[100];
	int nPoints, nFaces, nVertsInFace, face[4];
	glm::vec3 P;
	
	fin.open(filename.c_str());
	if(!fin.is_open())
		return false;
	
	// Check file format
	fin.getline(line, 100);
	if(strncmp(line, "ply", strlen("ply")))
		return false;
	fin.getline(line, 100);
	if(strncmp(line, "format ascii 1.0", strlen("format ascii 1.0")))
		return false;
	
	// Read # vertices, # faces & skip the rest of the header
	fin.getline(line, 100);
	while(strncmp(line, "end_header", strlen("end_header")) != 0)
	{
		if(strncmp(line, "element vertex", strlen("element vertex")) == 0)
			nPoints = atoi(&line[strlen("element vertex")]);
		else if(strncmp(line, "element face", strlen("element face")) == 0)
			nFaces = atoi(&line[strlen("element face")]);
		fin.getline(line, 100);
	}
	
	// Read vertices
	vertices.clear();
	for(int i=0; i<nPoints; i++)
	{
		fin >> P.x >> P.y >> P.z;
		vertices.push_back(P);
	}
	
	// Read faces
	triangles.clear();
	std::vector<int> tris;
	for(int i=0; i<nFaces; i++)
	{
		fin >> nVertsInFace;
		for(int j=0; j<nVertsInFace; j++)
			fin >> face[j];
		for(int j=1; j<nVertsInFace-1; j++)
		{
			tris.push_back(face[0]);
			tris.push_back(face[j]);
			tris.push_back(face[j+1]);
			Triangle t = Triangle(face[0], face[j], face[j+1]);
			triangles.push_back(t);
		}
	}
		std::cout << "size of tris: " << tris.size() << std::endl;
	return true;
}

void TriangleMesh::getBBox(glm::vec3 box[2]) const
{
	box[0] = vertices[0];
	box[1] = vertices[0];
	for(unsigned int i=1; i<vertices.size(); i++)
	{
		box[0] = glm::min(box[0], vertices[i]);
		box[1] = glm::max(box[1], vertices[i]);
	}
}

void TriangleMesh::transform(const glm::mat4 &matrix)
{
	for(unsigned int i=0; i<vertices.size(); i++)
		vertices[i] = glm::vec3(matrix * glm::vec4(vertices[i], 1.f));
}

void TriangleMesh::sendToOpenGL(ShaderProgram &program)
{
	vector<glm::vec3> *glVertices = new vector<glm::vec3>;
	vector<glm::vec3> *glNormals = new vector<glm::vec3>;
	vector<glm::vec4> *glColors = new vector<glm::vec4>;
	glm::vec3 normal;
	
	for(unsigned int i=0; i<triangles.size(); i++)
	{
		glVertices->push_back(vertices[triangles[i].getV1()]);
		glVertices->push_back(vertices[triangles[i].getV2()]);
		glVertices->push_back(vertices[triangles[i].getV3()]);
		
		normal = glm::cross(vertices[triangles[i].getV2()] - vertices[triangles[i].getV1()], vertices[triangles[i].getV3()] - vertices[triangles[i].getV1()]);
		normal = glm::normalize(normal);
		
		glNormals->push_back(normal);
		glNormals->push_back(normal);
		glNormals->push_back(normal);
	}
	
	// Send data to OpenGL
	if(!bGLObjsInit)
		glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	if(!bGLObjsInit)
	{
		glGenBuffers(1, &vboPosition);
		glGenBuffers(1, &vboNormal);
		glGenBuffers(1, &vboColor);
	}
	bGLObjsInit = true;
	// VBO for vertex positions
	glBindBuffer(GL_ARRAY_BUFFER, vboPosition);
	glBufferData(GL_ARRAY_BUFFER, 3 * glVertices->size() * sizeof(float), &(*glVertices)[0].x, GL_STATIC_DRAW);
	posLocation = program.bindVertexAttribute("position", 3, 0, 0);
	// VBO for vertex normals
	glBindBuffer(GL_ARRAY_BUFFER, vboNormal);
	glBufferData(GL_ARRAY_BUFFER, 3 * glNormals->size() * sizeof(float), &(*glNormals)[0].x, GL_STATIC_DRAW);
	normalLocation = program.bindVertexAttribute("normal", 3, 0, 0);
	// VBO for vertex colors
	// glBindBuffer(GL_ARRAY_BUFFER, vboColor);
	// glBufferData(GL_ARRAY_BUFFER, 4 * glColors->size() * sizeof(float), &(*glColors)[0].r, GL_STATIC_DRAW);
	// colorLocation = program.bindVertexAttribute("color", 4, 0, 0);

	delete glVertices;
	delete glNormals;
	delete glColors;
}

void TriangleMesh::render() const
{
	glBindVertexArray(vao);
	glEnableVertexAttribArray(posLocation);
	glEnableVertexAttribArray(normalLocation);
	// glEnableVertexAttribArray(colorLocation);
	glDrawArrays(GL_TRIANGLES, 0, 3*triangles.size());
}

void TriangleMesh::free()
{
	glDeleteBuffers(1, &vboPosition);
	glDeleteBuffers(1, &vboNormal);
	glDeleteBuffers(1, &vboColor);
	glDeleteVertexArrays(1, &vao);
	
	vertices.clear();
	triangles.clear();
	
	bGLObjsInit = false;
}

// void TriangleMesh::computeSATProjections() {
//     const int num_faces = triangles.size();

//     satProjections.resize(num_faces);

//     #if defined(OPENMP_VOLUMETRIC_RENDERER)
//     #pragma omp parallel for
//     #endif
//     for (int f = 0; f < num_faces; f+=3)
//     {
//         const glm::vec3& tv1o = vertices[triangles[f]];
//         const glm::vec3& tv2o = vertices[triangles[f+1]];
//         const glm::vec3& tv3o = vertices[triangles[f+2]];

//         // tv1 is the origin.
//         glm::vec3 tv1 = glm::vec3{0.f};
//         glm::vec3 tv2 = tv2o-tv1o;
//         glm::vec3 tv3 = tv3o-tv1o;

//         double e01[3];
//         e01[0] = tv2[0] - tv1[0];
//         e01[1] = tv2[1] - tv1[1];
//         e01[2] = tv2[2] - tv1[2];

//         double e12[3];
//         e12[0] = tv3[0] - tv2[0];
//         e12[1] = tv3[1] - tv2[1];
//         e12[2] = tv3[2] - tv2[2];

//         double e20[3];
//         e20[0] = tv1[0] - tv3[0];
//         e20[1] = tv1[1] - tv3[1];
//         e20[2] = tv1[2] - tv3[2];

//         const double axes_x[] = {1, 0, 0,       0,       0,       0,  e01[2],  e12[2],  e20[2], -e01[1], -e12[1], -e20[1]};
//         const double axes_y[] = {0, 1, 0, -e01[2], -e12[2], -e20[2],       0,       0,       0,  e01[0],  e12[0],  e20[0]};
//         const double axes_z[] = {0, 0, 1,  e01[1],  e12[1],  e20[1], -e01[0], -e12[0], -e20[0],       0,       0,       0};

//         SATProjection tproj = satProjections[f];

//         // Bbox.
//         double ptv[3];
//         for (int i = 0; i < 3; ++i)
//         {
//             ptv[0] = tv1[0] * axes_x[i] + tv1[1] * axes_y[i] + tv1[2] * axes_z[i];
//             ptv[1] = tv2[0] * axes_x[i] + tv2[1] * axes_y[i] + tv2[2] * axes_z[i];
//             ptv[2] = tv3[0] * axes_x[i] + tv3[1] * axes_y[i] + tv3[2] * axes_z[i];

//             auto result = std::minmax_element(ptv, ptv+3);
//             tproj.bbox[i][0] = result.first - ptv;
//             tproj.bbox[i][1] = result.second - ptv;
//         }

//         // Cross product axes.
//         for (int i = 3; i < 12; ++i)
//         {
//             ptv[0] = tv1[0] * axes_x[i] + tv1[1] * axes_y[i] + tv1[2] * axes_z[i];
//             ptv[1] = tv2[0] * axes_x[i] + tv2[1] * axes_y[i] + tv2[2] * axes_z[i];
//             ptv[2] = tv3[0] * axes_x[i] + tv3[1] * axes_y[i] + tv3[2] * axes_z[i];

//             auto result = std::minmax_element(ptv, ptv+3);
//             tproj.cross[i-3][0] = result.first - ptv;
//             tproj.cross[i-3][1] = result.second - ptv;
//         }
//     }
// }
