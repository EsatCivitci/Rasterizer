#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;

class Line {
public:
	Vec3 v0;
	Vec3 v1;
	Color c0;
	Color c1;
	Line(Vec3 v0, Vec3 v1, Color c0, Color c1) : v0(v0), v1(v1), c0(c0), c1(c1) {}

};

Matrix4 getTranslationMatrix(Translation &tr);
Matrix4 getRotationalMatrix(Rotation &tr);
Matrix4 getScalingMatrix(Scaling &s);
Matrix4 getProjectionMatrix(Camera *camera);
bool is_visible(double den, double num, double &t_e, double &t_l);
bool liang_barsky(Line &line);
Color interpolate_colors(Vec3 &target, Vec3 &v0, Vec3 &v1, Color &c0, Color &c1);
Matrix4 getViewportMatrix(Camera *camera);
void rasterize_line(Line &line, Camera *camera, Scene &scene, int triangle_index);
void swap(double &x0, double &x1, double &y0, double &y1, double &z0, double &z1, Color &c0, Color &c1);
double implicitLineEquation(Line line, double x, double y);

inline Vec4 convert3to4(Vec3 *v) {
	return Vec4(v->x, v->y, v->z, 1);
}

inline Vec3 convert4to3(Vec4 *v) {
	return Vec3(v->x, v->y, v->z);
}

inline Vec3 convert4to3_divide(Vec4 *v) {
	return Vec3(v->x/v->t, v->y/v->t, v->z/v->t);
}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

void Scene::initializeDepthBuffer(Camera *camera) {
	if (this->depth_buffer.empty()) {
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<double> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(__DBL_MAX__);
			}

			this->depth_buffer.push_back(rowOfColors);
		}
	}
	else {
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->depth_buffer[i][j] = __DBL_MAX__;
			}
		}
	}

}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}


/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "./magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}


/*
	Transformations, clipping, culling, rasterization are done here.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	// TODO: Implement this function
	Vec3 u = camera->u;
	Vec3 v = camera->v;
	Vec3 w = camera->w;
	Vec3 e = camera->position;

	// Camera modeling transformation
	double M_cam_values[4][4] = {{u.x,u.y,u.z,-(u.x*e.x + u.y*e.y + u.z*e.z)},
								{v.x,v.y,v.z,-(v.x*e.x + v.y*e.y + v.z*e.z)},
								{w.x,w.y,w.z, -(w.x*e.x + w.y*e.y + w.z*e.z)},
								{0,0,0,1}};
	Matrix4 M_cam(M_cam_values);
	// Projection Matrix
	Matrix4 M_projection = getProjectionMatrix(camera);
	Matrix4 M_project_camera = multiplyMatrixWithMatrix(M_projection, M_cam);
	//viewport matrix
	Matrix4 viewport_matrix = getViewportMatrix(camera);

	// Transformation matrix
	double identity_values[4][4] = {{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};

	int mesh_size = this->meshes.size();
	int vertices_size = this->vertices.size();
	vector<vector<Vec4>> transformed_vertices(mesh_size, vector<Vec4>(vertices_size));
	for(int i=0; i<mesh_size; i++) {
		Mesh *current_mesh = this->meshes[i];
		int transformation_size = current_mesh->numberOfTransformations;
		Matrix4 all_transformations_matrix(identity_values);


		// calculating all_transformation_matrix
		for(int j=0; j<transformation_size; j++) {

			if(current_mesh->transformationTypes[j] == 'r') {
				Matrix4 rotation_matrix = getRotationalMatrix(*this->rotations[current_mesh->transformationIds[j]-1]);
				all_transformations_matrix = multiplyMatrixWithMatrix(rotation_matrix, all_transformations_matrix);
			}
			else if(current_mesh->transformationTypes[j] == 't') {
				Matrix4 translation_matrix = getTranslationMatrix(*this->translations[current_mesh->transformationIds[j]-1]);
				all_transformations_matrix = multiplyMatrixWithMatrix(translation_matrix, all_transformations_matrix);
			}
			else if(current_mesh->transformationTypes[j] == 's') {
				Matrix4 scaling_matrix = getScalingMatrix(*this->scalings[current_mesh->transformationIds[j]-1]);
				all_transformations_matrix = multiplyMatrixWithMatrix(scaling_matrix, all_transformations_matrix);
			}
			else {
				return;
			}
		}


		Matrix4 transform_vertex_matrix = multiplyMatrixWithMatrix(M_project_camera, all_transformations_matrix);
		
		// applying transformations to the vertices
		for(int j=0; j<vertices_size; j++) {
			Vec4 curr_vertex = convert3to4(this->vertices[j]);
			transformed_vertices[i][j] = multiplyMatrixWithVec4(transform_vertex_matrix, curr_vertex);
		}

		int triangles_size = current_mesh->triangles.size();
		for(int j=0; j<triangles_size; j++) {
			int v0_id = current_mesh->triangles[j].vertexIds[0]-1;
			int v1_id = current_mesh->triangles[j].vertexIds[1]-1;
			int v2_id = current_mesh->triangles[j].vertexIds[2]-1;
	
			double t0 = transformed_vertices[i][v0_id].t;
			double t1 = transformed_vertices[i][v1_id].t;
			double t2 = transformed_vertices[i][v2_id].t;

			Vec3 v0 = convert4to3(&transformed_vertices[i][v0_id]);
			Vec3 v1 = convert4to3(&transformed_vertices[i][v1_id]);
			Vec3 v2 = convert4to3(&transformed_vertices[i][v2_id]);

			Color c0 = *(this->colorsOfVertices[v0_id]);
			Color c1 = *(this->colorsOfVertices[v1_id]);
			Color c2 = *(this->colorsOfVertices[v2_id]);

			if(this -> cullingEnabled) {
				Vec3 edge1 = subtractVec3(v1, v0);
				Vec3 edge2 = subtractVec3(v2, v0);
				Vec3 normal = normalizeVec3(crossProductVec3(edge1, edge2));

				// PONDER OVER IT
				if(dotProductVec3(normal, v0) < 0) {
					continue;
				}
			}

			// perspective divide
			v0 = multiplyVec3WithScalar(v0, 1/t0);
			v1 = multiplyVec3WithScalar(v1, 1/t1);
			v2 = multiplyVec3WithScalar(v2, 1/t2);

			
			Line line1(v0,v1,c0,c1);
			Line line2(v1,v2,c1,c2);
			Line line3(v2,v0,c2,c0);

			//Wireframe
			if (current_mesh->type == 0){
				// clipping

				// check if there is any problem here
				bool visible1 = liang_barsky(line1);
				if (visible1){
					line1.c0 = interpolate_colors(line1.v0, v0, v1, c0, c1);
					line1.c1 = interpolate_colors(line1.v1, v0, v1, c0, c1);

					//multiplying vertices with viewport matrix
					Vec4 v0_vp_4 = multiplyMatrixWithVec4(viewport_matrix, convert3to4(&line1.v0));
					Vec4 v1_vp_4 = multiplyMatrixWithVec4(viewport_matrix, convert3to4(&line1.v1));

					line1.v0 = convert4to3(&v0_vp_4);
					line1.v1 = convert4to3(&v1_vp_4);

					rasterize_line(line1, camera, *this, j);

				}
				
				bool visible2 = liang_barsky(line2);
				if (visible2){
					line2.c0 = interpolate_colors(line2.v0, v1, v2, c1, c2);
					line2.c1 = interpolate_colors(line2.v1, v1, v2, c1, c2);

					//multiplying vertices with viewport matrix
					Vec4 v0_vp_4 = multiplyMatrixWithVec4(viewport_matrix, convert3to4(&line2.v0));
					Vec4 v1_vp_4 = multiplyMatrixWithVec4(viewport_matrix, convert3to4(&line2.v1));

					line2.v0 = convert4to3(&v0_vp_4);
					line2.v1 = convert4to3(&v1_vp_4);

					rasterize_line(line2, camera, *this, j);
				}
				
				bool visible3 = liang_barsky(line3);
				if (visible3){
					line1.c0 = interpolate_colors(line3.v0, v2, v0, c2, c0);
					line1.c1 = interpolate_colors(line3.v1, v2, v0, c2, c0);

					//multiplying vertices with viewport matrix
					Vec4 v0_vp_4 = multiplyMatrixWithVec4(viewport_matrix, convert3to4(&line3.v0));
					Vec4 v1_vp_4 = multiplyMatrixWithVec4(viewport_matrix, convert3to4(&line3.v1));

					line3.v0 = convert4to3(&v0_vp_4);
					line3.v1 = convert4to3(&v1_vp_4);

					rasterize_line(line3, camera, *this, j);
				}
			}
			//Solid mesh
			else{
				//multiplying vertices with viewport matrix
				Vec4 v0_vp4 = multiplyMatrixWithVec4(viewport_matrix, convert3to4(&v0));
				Vec3 v0_vp = convert4to3(&v0_vp4);

				Vec4 v1_vp4 = multiplyMatrixWithVec4(viewport_matrix, convert3to4(&v1));
				Vec3 v1_vp = convert4to3(&v1_vp4);

				Vec4 v2_vp4 = multiplyMatrixWithVec4(viewport_matrix, convert3to4(&v2));
				Vec3 v2_vp = convert4to3(&v2_vp4);

				line1.v0 = v0_vp;
				line1.v1 = v1_vp;
				
				line2.v0 = v1_vp;
				line2.v1 = v2_vp;

				line3.v0 = v2_vp;
				line3.v1 = v0_vp;

				//finding x_max, x_min, y_max, y_min
				double x_max = max(v0_vp.x,max(v1_vp.x, v2_vp.x));
				double x_min = min(v0_vp.x,min(v1_vp.x, v2_vp.x));
				double y_max = max(v0_vp.y,max(v1_vp.y, v2_vp.y));
				double y_min = min(v0_vp.y, min(v1_vp.y, v2_vp.y));
				
				double alpha_denom = implicitLineEquation(line2, v0_vp.x, v0_vp.y); //f12(x0,y0)
				double beta_denom = implicitLineEquation(line3, v1_vp.x, v1_vp.y); //f20(x1,y1)
				double gama_denom = implicitLineEquation(line1, v2_vp.x, v2_vp.y); //f01(x2,y2)
				
				for (int y = y_min; y < y_max; y++){
					for (int x = x_min; x < x_max; x++){
						if (y >= camera->verRes || y < 0 || x >= camera->horRes || x < 0) continue;


						double alpha = implicitLineEquation(line2, x, y) / alpha_denom;
						double beta = implicitLineEquation(line3, x, y) / beta_denom;
						double gama = implicitLineEquation(line1, x, y) / gama_denom;

						if (alpha >= 0 && beta >= 0 && gama >= 0){
							//depth buffer check
							double depth_value = alpha*v0.z + beta*v1.z + gama*v2.z;
							if (depth_buffer[x][y] > depth_value){
								depth_buffer[x][y] = depth_value;
								Color c =addColors(scaleColor(c0, alpha),addColors(scaleColor(c1, beta), scaleColor(c2, gama)));
								this->image[x][y] = c;
							}

						}

					}
				}
			}

			// WHAT IF Z = 0 IN PERSPECTIVE DIVIDE . THINK ABOUT IT
		}
	}
}

double implicitLineEquation(Line line, double x, double y){
	double x0 = line.v0.x;
	double y0 = line.v0.y;
	double x1 = line.v1.x;
	double y1 = line.v1.y;

	return x*(y0-y1) + y*(x1-x0) + x0*y1 - y0*x1;
}

void rasterize_line(Line &line, Camera *camera, Scene &scene, int triangle_index) {
	double y0 = line.v0.y;
	double y1 = line.v1.y;
	double x0 = line.v0.x;
	double x1 = line.v1.x;
	double z0 = line.v0.z;
	double z1 = line.v1.z;

	Color c0 = line.c0;
	Color c1 = line.c1;
	
	bool is_parallel = false;
	int scale_factor = 1;
	
	// DÜZELT
	if(line.v1.x - line.v0.x == 0) {
		// parallel to y axis
		// according to y
		if(y0 > y1) {
			swap(x0, x1, y0, y1, z0, z1, c0, c1);
		}
		if(x0 > x1) {
			scale_factor = -1;
		}
		
		Color c = c0;
		double x = x0;
		double d = scale_factor*0.5*(y0-y1) + (x1-x0);
		Color dc =  scaleColor(subtractColors(c1,c0),1/(y1-y0));

		for(int i=y0; i<=y1; i++) {
			if((int) x == 208 && i == 300) {
				int asdasd = 12312;
			}
			double alpha_z = (double)sqrt((pow((i-y0),2) + pow(x-x0, 2)) /(pow((y1-y0),2)+pow((x1-x0),2)));
			double interpolated_z = (z1-z0)*alpha_z+z0;


			if(interpolated_z < scene.depth_buffer[(int)x][i]) {
				scene.image[(int)x][i] = c;	
				scene.depth_buffer[(int)x][i] = interpolated_z;
			}
			
			if(scale_factor*d > 0) {
				x += scale_factor;
				d += scale_factor*(y0-y1) + (x1-x0);
			}
			else {
				d += (x1-x0);
			}
			c = addColors(c, dc);
		}
	}
	else {
		double slope = abs((y1 - y0)/(x1 - x0));
		if(slope <= 1) {
			// according to x
			if(x0 > x1) {
				swap(x0, x1, y0, y1, z0, z1, c0, c1);
			}
			if(y0 > y1) {
				scale_factor = -1;
			}

			Color c = c0;
			double y = y0;
			double d = (y0-y1) + scale_factor*0.5*(x1-x0);
			Color dc =  scaleColor(subtractColors(c1,c0),1/(x1-x0));

			for(int i=x0; i<=x1; i++) {
				if((int) i == 208 && (int)y == 300) {
					int asdasd = 12312;
				}
				double alpha_z = (double)sqrt((pow((i-x0),2) + pow(y-y0, 2)) /(pow((y1-y0),2)+pow((x1-x0),2)));
				double interpolated_z = (z1-z0)*alpha_z+z0 ;
				if(interpolated_z < scene.depth_buffer[i][(int)y]) {
					scene.image[i][(int)y] = c;	
					scene.depth_buffer[i][(int)y] = interpolated_z;
				}

				if(scale_factor * d < 0) {
					y += scale_factor;
					d += (y0-y1) + scale_factor*(x1-x0);
				}
				else {
					d += (y0-y1);
				}
				c = addColors(c, dc);
			}
		}
		else if(slope > 1){
			// according to y
			if(y0 > y1) {
				swap(x0, x1, y0, y1, z0, z1, c0, c1);
			}
			if(x0 > x1) {
				scale_factor = -1;
			}
			
			Color c = c0;
			double x = x0;
			double d = scale_factor*0.5*(y0-y1) + (x1-x0);
			Color dc =  scaleColor(subtractColors(c1,c0),1/(y1-y0));

			for(int i=y0; i<=y1; i++) {
				if((int) x == 208 && i == 300) {
					int asdasd = 12312;
				}

				double alpha_z = (double)sqrt((pow((i-y0),2) + pow(x-x0, 2)) /(pow((y1-y0),2)+pow((x1-x0),2)));
				double interpolated_z = (z1-z0)*alpha_z+z0;
				if(interpolated_z < scene.depth_buffer[(int)x][i]) {
					scene.image[(int)x][i] = c;	
					scene.depth_buffer[(int)x][i] = interpolated_z;
				}

				if(scale_factor*d > 0) {
					x += scale_factor;
					d += scale_factor*(y0-y1) + (x1-x0);
				}
				else {
					d += (x1-x0);
				}
				c = addColors(c, dc);
			}
		}
	}
}

void swap(double &x0, double &x1, double &y0, double &y1, double &z0, double &z1, Color &c0, Color &c1) {
	double tmp_x;
	tmp_x = x1;
	x1 = x0;
	x0 = tmp_x;

	double tmp_y;
	tmp_y = y1;
	y1 = y0;
	y0 = tmp_y;

	Color tmp_c0;
	tmp_c0 = c1;
	c1 = c0;
	c0 = tmp_c0;
}

Color interpolate_colors(Vec3 &target, Vec3 &v0, Vec3 &v1, Color &c0, Color &c1) {
	// ROUND THEM LATER
	double alpha;
	if(v1.x - v0.x == 0) {
		// parallel to y axis
		alpha = (target.y - v0.y) / (v1.y - v0.y);
	}
	else {
		double slope = abs((v1.y - v0.y) / (v1.x - v0.x));

		if(slope > 1) {
			// according to y
			alpha = (target.y - v0.y) / (v1.y - v0.y);
		}
		else {
			// according to x
			alpha = (target.x - v0.x) / (v1.x - v0.x);
		}
	}

	return addColors(scaleColor(c0, 1-alpha), scaleColor(c1, alpha));
}

bool liang_barsky(Line &line) {
	double t_e = 0, t_l = 1;
	bool visible = false;
	double x_min = -1, y_min = -1, z_min = -1;
	double x_max = 1, y_max = 1, z_max = 1;
	double d_x = line.v1.x - line.v0.x;
	double d_y = line.v1.y - line.v0.y;
	double d_z = line.v1.z - line.v0.z;

	if(is_visible(d_x, x_min - line.v0.x, t_e, t_l)) { // left
		if(is_visible(-d_x, line.v0.x - x_max, t_e, t_l)) { // right
			if(is_visible(d_y, y_min - line.v0.y, t_e, t_l)) { // bottom
				if(is_visible(-d_y, line.v0.y - y_max, t_e, t_l)) { // top
					if(is_visible(d_z, z_min - line.v0.z, t_e, t_l)) { // front
						if(is_visible(-d_z, line.v0.z - z_max, t_e, t_l)) { // back
							visible = true;
							if(t_l < 1) {
								line.v1.x = line.v0.x + d_x*t_l;
								line.v1.y = line.v0.y + d_y*t_l;
								line.v1.z = line.v0.z + d_z*t_l;
							}
							if(t_e > 0) {
								line.v0.x = line.v0.x + d_x*t_e;
								line.v0.y = line.v0.y + d_y*t_e;
								line.v0.z = line.v0.z + d_z*t_e;
							}
						}
					}
				}
			}
		}
	}

	return visible;
}

bool is_visible(double den, double num, double &t_e, double &t_l) {
	if(den > 0) {
		double t = num/den;
		if(t > t_l) {
			return false;
		}
		if(t > t_e) {
			t_e = t;
		}
	}
	else if(den < 0) {
		double t = num/den;
		if(t < t_e) {
			return false;
		}
		if(t <  t_l) {
			t_l = t;
		}
	}
	else if(num > 0) {
		return false;
	}
	return true;
}

Matrix4 getViewportMatrix(Camera *camera) {
	int n_x = camera->horRes;
	int n_y = camera->verRes;
	
	// normally 3x4 in the slides
	double values[4][4] = {{(double) n_x/2, 0, 0, (double)(n_x-1)/2},
							{0, (double) n_y/2, 0, (double)(n_y-1)/2},
							{0, 0, 0.5, 0.5},
							{0, 0, 0, 1}};
	return Matrix4(values);
}

Matrix4 getProjectionMatrix(Camera *camera) {
	double r = camera->right;
	double l = camera->left;
	double b = camera->bottom;
	double t = camera->top;
	double n = camera->near;
	double f = camera->far;

	Matrix4 projectionMatrix;
	if(camera->projectionType == 0) {
		//ortographic
		double values[4][4] = {{2/(r-l),0,0,-(r+l)/(r-l)},
								{0, 2/(t-b),0,-(t+b)/(t-b)},
								{0,0,-2/(f-n),-(f+n)/(f-n)},
								{0,0,0,1}};

		projectionMatrix = Matrix4(values);
	}
	else if(camera->projectionType == 1) {
		//perspective
		double values[4][4] = {{2*n/(r-l),0,(r+l)/(r-l),0},
								{0, 2*n/(t-b),(t+b)/(t-b),0},
								{0,0,-(f+n)/(f-n),-2*f*n/(f-n)},
								{0,0,-1,0}};

		projectionMatrix = Matrix4(values);
	}
	return projectionMatrix;
}

Matrix4 getTranslationMatrix(Translation &tr) {
	Vec4 v0(1,0,0,tr.tx);
	Vec4 v1(0,1,0,tr.ty);
	Vec4 v2(0,0,1,tr.tz);
	Vec4 v3(0,0,0,1);

	double values[4][4] = {{1,0,0,tr.tx},
		{0,1,0,tr.ty},
		{0,0,1,tr.tz},
		{0,0,0,1}};

	return Matrix4(values);
}

Matrix4 getRotationalMatrix(Rotation &tr) {

	if (tr.ux == 0 && tr.uy == 0 && tr.uz == 0){
			double alpha = tr.angle * M_PI / 180;
			double R_x_values[4][4] = {{1,0,0,0},
				{0, cos(alpha), -sin(alpha), 0},
				{0,sin(alpha),cos(alpha),0},
				{0,0,0,1}};
			Matrix4 R_x(R_x_values);
			return R_x;
	}

	Vec3 u(tr.ux, tr.uy, tr.uz);
	Vec3 v;
	double min_val = min(tr.ux, min(tr.uy, tr.uz));
	if (min_val == tr.ux)
		v = {0, -tr.uz, tr.uy};
	else if (min_val == tr.uy)
		v = {-tr.uz, 0, tr.ux};
	else if (min_val == tr.uy)
		v = {-tr.uy, tr.ux, 0};
	
	Vec3 w = crossProductVec3(u,v);

	u = normalizeVec3(u);
	v = normalizeVec3(v);
	w = normalizeVec3(w);

	double M_values[4][4] = {{u.x, u.y, u.z, 0},
		{v.x, v.y, v.z, 0},
		{w.x, w.y, w.z, 0},
		{0,0,0,1}};
	Matrix4 M(M_values);
				
	double alpha = tr.angle * M_PI / 180;
	double R_x_values[4][4] = {{1,0,0,0},
		{0, cos(alpha), -sin(alpha), 0},
		{0,sin(alpha),cos(alpha),0},
		{0,0,0,1}};
	Matrix4 R_x(R_x_values);

	double M_T_values[4][4] = {{u.x, v.x, w.x, 0},
		{u.y, v.y, w.y, 0},
		{u.z, v.z, w.z, 0},
		{0,0,0,1}};
	Matrix4 M_T(M_T_values);

	Matrix4 rotationMatrix = multiplyMatrixWithMatrix(multiplyMatrixWithMatrix(M_T, R_x), M);
	return rotationMatrix;
}

Matrix4 getScalingMatrix(Scaling &s) {
	double values[4][4] = {{s.sx,0,0,0},
		{0,s.sy,0,0},
		{0,0,s.sz,0},
		{0,0,0,1}};

	return Matrix4(values);
}