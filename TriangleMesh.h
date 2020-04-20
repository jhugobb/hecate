#ifndef _TRIANGLE_MESH_INCLUDE
#define _TRIANGLE_MESH_INCLUDE


#include <vector>
#include <glm/glm.hpp>
#include "ShaderProgram.h"
#include "geometry/SATProjection.h"
#include "geometry/Triangle.h"

class TriangleMesh
{
public:
	TriangleMesh();

	void init(const vector<glm::vec3> &newVertices, const vector<Triangle*> &newTriangles);
	bool load(const string &filename);

	void getBBox(glm::vec3 box[2]) const;
	void transform(const glm::mat4 &matrix);

	void sendToOpenGL(ShaderProgram &program);
	void render() const;
	void free();
	
	const vector<glm::vec3> &getVertices() const { return vertices; }
	vector<glm::vec3> &getVertices() { return vertices; }
	const vector<Triangle*> &getTriangles() const { return triangles; }
	vector<Triangle*> &getTriangles() { return triangles; }

	
	static int next(int corner);
	static int previous(int corner);

private:
	void computePerVertexNormals();
	void computeCornerTable();
	void computeSATProjections();

private:
	vector<glm::vec3> vertices;
	vector<Triangle*> triangles;

	bool bGLObjsInit;
	GLuint vao;
	GLuint vboPosition, vboNormal, vboColor;
	GLint posLocation, normalLocation, colorLocation;

};


#endif // _TRIANGLE_MESH_INCLUDE



