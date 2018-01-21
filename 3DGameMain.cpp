#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if defined(__APPLE__)
#include <GLUT/GLUT.h>
#include <OpenGL/gl3.h>
#include <OpenGL/glu.h>
#else

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__)
#include <windows.h>
#endif
#include <GL/glew.h>		 
#include <GL/freeglut.h>	
#endif

#include <string>
#include <vector>
#include <fstream>
#include <algorithm> 

unsigned int windowWidth = 800, windowHeight = 800;
unsigned char keyPressed[256];

// OpenGL major and minor versions
int majorVersion = 3, minorVersion = 0;

void getErrorInfo(unsigned int handle) 
{
	int logLen;
	glGetShaderiv(handle, GL_INFO_LOG_LENGTH, &logLen);
	if (logLen > 0) 
	{
		char * log = new char[logLen];
		int written;
		glGetShaderInfoLog(handle, logLen, &written, log);
		printf("Shader log:\n%s", log);
		delete log;
	}
}

// check if shader could be compiled
void checkShader(unsigned int shader, char * message) 
{
	int OK;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &OK);
	if (!OK) 
	{
		printf("%s!\n", message);
		getErrorInfo(shader);
	}
}

// check if shader could be linked
void checkLinking(unsigned int program) 
{
	int OK;
	glGetProgramiv(program, GL_LINK_STATUS, &OK);
	if (!OK) 
	{
		printf("Failed to link shader program!\n");
		getErrorInfo(program);
	}
}
//

// shader program for rendering textured meshes

const char *vertexSource0 = "\n\
#version 410 \n\
precision highp float; \n\
in vec3 vertexPosition; \n\
in vec2 vertexTexCoord; \n\
in vec3 vertexNormal; \n\
uniform mat4 M, MInv, MVP; \n\
uniform vec3 wEye; \n\
uniform vec4 worldLightPosition; \n\
out vec2 texCoord; \n\
out vec3 worldNormal; \n\
out vec3 worldView; \n\
out vec3 worldLight; \n\
void main() { \n\
texCoord = vertexTexCoord; \n\
vec4 worldPosition = \n\
vec4(vertexPosition, 1) * M; \n\
worldLight  = \n\
worldLightPosition.xyz * worldPosition.w - worldPosition.xyz * worldLightPosition.w; \n\
worldView = wEye - worldPosition.xyz; \n\
worldNormal = (MInv * vec4(vertexNormal, 0.0)).xyz; \n\
gl_Position = vec4(vertexPosition, 1) * MVP; \n\
} \n\
";

const char *fragmentSource0 = "\n\
#version 410 \n\
precision highp float; \n\
uniform sampler2D samplerUnit; \n\
uniform vec3 La, Le; \n\
uniform vec3 ka, kd, ks; \n\
uniform float shininess; \n\
uniform samplerCube environmentMap;\n\
in vec2 texCoord; \n\
in vec3 worldNormal; \n\
in vec3 worldView; \n\
in vec3 worldLight; \n\
out vec4 fragmentColor; \n\
void main() {\n\
    vec3 N = normalize(worldNormal);\n\
    vec3 V = normalize(worldView);\n\
    vec3 L = normalize(worldLight);\n\
    vec3 H = normalize((V + L) / 2);\n\
    vec3 R = N * dot(N, V) * 2 - V;\n\
    vec3 texel = texture(samplerUnit, texCoord).xyz;\n\
    vec3 texel1 = texture(environmentMap, R).xyz;\n\
    vec3 color =\n\
    La * ka +\n\
    Le * kd * texel* max(0.0, dot(L, N)) +\n\
    Le * ks * pow(max(0.0, dot(H, N)), shininess);\n\
    fragmentColor = vec4(color * 0.9 + texel1 * 0.1, 1);\n\
} \n\
";

// shader program for rendering the ground as an infinite quad

const char *vertexSource1 = "\n\
#version 410 \n\
precision highp float; \n\
\n\
in vec4 vertexPosition; \n\
in vec2 vertexTexCoord; \n\
in vec3 vertexNormal; \n\
uniform mat4 M, MInv, MVP; \n\
\n\
out vec2 texCoord; \n\
out vec4 worldPosition; \n\
out vec3 worldNormal; \n\
\n\
void main() { \n\
texCoord = vertexTexCoord; \n\
worldPosition = vertexPosition * M; \n\
worldNormal = (MInv * vec4(vertexNormal, 0.0)).xyz; \n\
gl_Position = vertexPosition * MVP; \n\
} \n\
";


const char *fragmentSource1 = "\n\
#version 410 \n\
precision highp float; \n\
uniform sampler2D samplerUnit; \n\
uniform vec3 La, Le; \n\
uniform vec3 ka, kd, ks; \n\
uniform float shininess; \n\
uniform vec3 wEye; \n\
uniform vec4 worldLightPosition; \n\
in vec2 texCoord; \n\
in vec4 worldPosition; \n\
in vec3 worldNormal; \n\
out vec4 fragmentColor; \n\
void main() { \n\
vec3 N = normalize(worldNormal); \n\
vec3 V = normalize(wEye * worldPosition.w - worldPosition.xyz);\n\
vec3 L = normalize(worldLightPosition.xyz * worldPosition.w - worldPosition.xyz * worldLightPosition.w);\n\
vec3 H = normalize(V + L); \n\
vec2 position = worldPosition.xz / worldPosition.w; \n\
vec2 tex = position.xy - floor(position.xy); \n\
vec3 texel = texture(samplerUnit, tex).xyz; \n\
vec3 color = vec3(10, 10, 0);\n\
if((position.x + position.y)>0.1){\n\
color = vec3(10, 0, 0); \n\
}\n\
if((position.x + position.y) < 0){\n\
color = vec3(10, 10, 0); \n\
}\n\
if((position.x + position.y) > 0.2){\n\
color = vec3(10, 10, 10); \n\
color = (La * ka + Le * kd * texel* max(0.0, dot(L, N)) + Le * ks * pow(max(0.0, dot(H, N)), shininess)); \n\
}\n\
fragmentColor = vec4(color, 1); \n\
} \n\
";


// shader program for rendering plane-projected shadows

const char *vertexSource2 = " \n\
#version 410 \n\
precision highp float; \n\
in vec3 vertexPosition; \n\
in vec2 vertexTexCoord; \n\
in vec3 vertexNormal; \n\
uniform mat4 M, MVP; \n\
uniform vec4 worldLightPosition; \n\
void main() { \n\
vec4 p = vec4(vertexPosition, 1) * M; \n\
vec3 s; \n\
s.y = -0.999; \n\
s.x = (p.x - worldLightPosition.x) / (p.y - worldLightPosition.y) * (s.y - worldLightPosition.y) + worldLightPosition.x; \n\
s.z = (p.z - worldLightPosition.z) / (p.y - worldLightPosition.y) * (s.y - worldLightPosition.y) + worldLightPosition.z; \n\
gl_Position = vec4(s, 1) * MVP; \n\
} \n\
";

const char *fragmentSource2 = " \n\
#version 410 \n\
precision highp float; \n\
uniform sampler2D samplerUnit; \n\
uniform samplerCube environmentMap; \n\
out vec4 fragmentColor; \n\
void main() { \n\
fragmentColor = vec4(0.1, 0.2, 0.1, 1); \n\
} \n\
";

///
//Shader for environment map
const char *vertexSource3 = " \n\
#version 410\n\
precision highp float;\n\
in vec4 vertexPosition;\n\
uniform mat4 viewDirMatrix;\n\
out vec4 worldPosition;\n\
out vec3 viewDir;\n\
void main(){\n\
worldPosition = vertexPosition;\n\
viewDir = (vertexPosition * viewDirMatrix).xyz;\n\
gl_Position = vertexPosition;\n\
gl_Position.z = 0.999999;\n\
}\n\
";

const char *fragmentSource3 = "\n\
#version 410\n\
precision highp float;\n\
uniform sampler2D samplerUnit;\n\
uniform samplerCube  environmentMap;\n\
in vec4 worldPosition;\n\
in vec3 viewDir;\n\
out vec4 fragmentColor;\n\
void main(){\n\
vec3 texel = texture(environmentMap, viewDir).xyz;\n\
fragmentColor = vec4(texel, 1);\n\
}\n\
";
///




////////////////////////////

struct vec2
{
float x, y;

	vec2(float x = 0.0, float y = 0.0) : x(x), y(y) {}

	static vec2 random() { return vec2(((float)rand() / RAND_MAX) * 2 - 1, ((float)rand() / RAND_MAX) * 2 - 1); }

	vec2 operator+(const vec2& v) { return vec2(x + v.x, y + v.y); }

	vec2 operator-(const vec2& v) { return vec2(x - v.x, y - v.y); }

	vec2 operator*(float s) { return vec2(x * s, y * s); }

	vec2 operator/(float s) { return vec2(x / s, y / s); }

	float length() { return sqrt(x * x + y * y); }

	vec2 normalize() { return *this / length(); }
};


struct vec3 
{
	float x, y, z;

	vec3(float x = 0.0, float y = 0.0, float z = 0.0) : x(x), y(y), z(z) {}

	static vec3 random() { return vec3(((float)rand() / RAND_MAX) * 2 - 1, ((float)rand() / RAND_MAX) * 2 - 1, ((float)rand() / RAND_MAX) * 2 - 1); }

	vec3 operator+(const vec3& v) { return vec3(x + v.x, y + v.y, z + v.z); }

	vec3 operator-(const vec3& v) { return vec3(x - v.x, y - v.y, z - v.z); }

	vec3 operator*(float s) { return vec3(x * s, y * s, z * s); }

	vec3 operator/(float s) { return vec3(x / s, y / s, z / s); }
    
    vec3 operator*=(const float f)
    {
        return vec3(x *= f, y *= f, z *=f);
    }
    
    vec3 operator+=(const float f)
    {
        return vec3(x += f, y += f, z +=f);
    }
    
    vec3 operator-=(const float f)
    {
        return vec3(x -= f, y -= f, z -=f);
    }
    
    vec3 operator*=(const vec3& v)
    {
        return vec3(x *= v.x, y *= v.y, z *= v.z);
    }
    
    vec3 operator+=(const vec3& v)
    {
        return vec3(x += v.x, y += v.y, z += v.z);
    }

    vec3 operator-=(const vec3& v)
    {
        return vec3(x -= v.x, y -= v.y, z -= v.z);
    }
    
	float length() { return sqrt(x * x + y * y + z * z); }

	vec3 normalize() { return *this / length(); }

	void print() { printf("%f \t %f \t %f \n", x, y, z); }
};

vec3 cross(const vec3& a, const vec3& b)
{
	return vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x );
}


// row-major matrix 4x4
struct mat4 
{
	float m[4][4];
public:
	mat4() {}
	mat4(float m00, float m01, float m02, float m03,
		float m10, float m11, float m12, float m13,
		float m20, float m21, float m22, float m23,
		float m30, float m31, float m32, float m33) 
	{
		m[0][0] = m00; m[0][1] = m01; m[0][2] = m02; m[0][3] = m03;
		m[1][0] = m10; m[1][1] = m11; m[1][2] = m12; m[1][3] = m13;
		m[2][0] = m20; m[2][1] = m21; m[2][2] = m22; m[2][3] = m23;
		m[3][0] = m30; m[3][1] = m31; m[3][2] = m32; m[3][3] = m33;
	}

	mat4 operator*(const mat4& right) 
	{
		mat4 result;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				result.m[i][j] = 0;
				for (int k = 0; k < 4; k++) result.m[i][j] += m[i][k] * right.m[k][j];
			}
		}
		return result;
	}

	operator float*() { return &m[0][0]; }
};

// 3D point in homogeneous coordinates
struct vec4 
{
	float v[4];

	vec4(float x = 0, float y = 0, float z = 0, float w = 1) 
	{
		v[0] = x; v[1] = y; v[2] = z; v[3] = w;
	}

	vec4 operator*(const mat4& mat) 
	{
		vec4 result;
		for (int j = 0; j < 4; j++) 
		{
			result.v[j] = 0;
			for (int i = 0; i < 4; i++) result.v[j] += v[i] * mat.m[i][j];
		}
		return result;
	}

	vec4 operator+(const vec4& vec) 
	{
		vec4 result(v[0] + vec.v[0], v[1] + vec.v[1], v[2] + vec.v[2], v[3] + vec.v[3]);
		return result;
	}
};



class   Mesh
{
	struct  Face
	{
		int       positionIndices[4];
		int       normalIndices[4];
		int       texcoordIndices[4];
		bool      isQuad;
	};

	std::vector<std::string*>	rows;
	std::vector<vec3*>		positions;
	std::vector<std::vector<Face*> >          submeshFaces;
	std::vector<vec3*>		normals;
	std::vector<vec2*>		texcoords;

	unsigned int vao;	
	int nTriangles;

public:
	Mesh(const char *filename);
	~Mesh();

	void DrawModel();
};


Mesh::Mesh(const char *filename)
{
	std::fstream file(filename); 
	if(!file.is_open())       
	{
		return;
	}

	char buffer[256];
	while(!file.eof())
	{
		file.getline(buffer,256);
		rows.push_back(new std::string(buffer));
	}

	submeshFaces.push_back(std::vector<Face*>());
	std::vector<Face*>* faces = &submeshFaces.at(submeshFaces.size()-1);

	for(int i = 0; i < rows.size(); i++)
	{
		if(rows[i]->empty() || (*rows[i])[0] == '#') 
			continue;      
		else if((*rows[i])[0] == 'v' && (*rows[i])[1] == ' ')
		{
			float tmpx,tmpy,tmpz;
			sscanf(rows[i]->c_str(), "v %f %f %f" ,&tmpx,&tmpy,&tmpz);      
			positions.push_back(new vec3(tmpx,tmpy,tmpz));  
		}
		else if((*rows[i])[0] == 'v' && (*rows[i])[1] == 'n')    
		{
			float tmpx,tmpy,tmpz;   
			sscanf(rows[i]->c_str(), "vn %f %f %f" ,&tmpx,&tmpy,&tmpz);
			normals.push_back(new vec3(tmpx,tmpy,tmpz));     
		}
		else if((*rows[i])[0] == 'v' && (*rows[i])[1] == 't')
		{
			float tmpx,tmpy;
			sscanf(rows[i]->c_str(), "vt %f %f" ,&tmpx,&tmpy);
			texcoords.push_back(new vec2(tmpx,tmpy));     
		}
		else if((*rows[i])[0] == 'f')  
		{
			if(count(rows[i]->begin(),rows[i]->end(), ' ') == 3)
			{
				Face* f = new Face();
				f->isQuad = false;
				sscanf(rows[i]->c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d",
					&f->positionIndices[0], &f->texcoordIndices[0], &f->normalIndices[0],
					&f->positionIndices[1], &f->texcoordIndices[1], &f->normalIndices[1],
					&f->positionIndices[2], &f->texcoordIndices[2], &f->normalIndices[2]);
				faces->push_back(f);
			}
			else
			{
				Face* f = new Face();
				f->isQuad = true;
				sscanf(rows[i]->c_str(), "f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d", 
					&f->positionIndices[0], &f->texcoordIndices[0], &f->normalIndices[0],
					&f->positionIndices[1], &f->texcoordIndices[1], &f->normalIndices[1],
					&f->positionIndices[2], &f->texcoordIndices[2], &f->normalIndices[2],
					&f->positionIndices[3], &f->texcoordIndices[3], &f->normalIndices[3]);
				faces->push_back(f);   
			}
		}
		else if((*rows[i])[0] == 'g')
		{
			if(faces->size() > 0)
			{
				submeshFaces.push_back(std::vector<Face*>());
				faces = &submeshFaces.at(submeshFaces.size()-1);
			}
		}
	}
	
	int numberOfTriangles = 0;
	for(int iSubmesh=0; iSubmesh<submeshFaces.size(); iSubmesh++)
	{
		std::vector<Face*>& faces = submeshFaces.at(iSubmesh);

		for(int i=0;i<faces.size();i++)
		{
			if(faces[i]->isQuad) numberOfTriangles += 2;
			else numberOfTriangles += 1;
		}
	}

	nTriangles = numberOfTriangles;
	
	float *vertexCoords = new float[numberOfTriangles * 9];
	float *vertexTexCoords = new float[numberOfTriangles * 6];
	float *vertexNormalCoords = new float[numberOfTriangles * 9];


	int triangleIndex = 0;
	for(int iSubmesh=0; iSubmesh<submeshFaces.size(); iSubmesh++)
	{
		std::vector<Face*>& faces = submeshFaces.at(iSubmesh);

		for(int i=0;i<faces.size();i++)
		{
			if(faces[i]->isQuad) 
			{
				vertexTexCoords[triangleIndex * 6] =     texcoords[faces[i]->texcoordIndices[0]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 1] = 1-texcoords[faces[i]->texcoordIndices[0]-1]->y;
				
				vertexTexCoords[triangleIndex * 6 + 2] = texcoords[faces[i]->texcoordIndices[1]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 3] = 1-texcoords[faces[i]->texcoordIndices[1]-1]->y;

				vertexTexCoords[triangleIndex * 6 + 4] = texcoords[faces[i]->texcoordIndices[2]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 5] = 1-texcoords[faces[i]->texcoordIndices[2]-1]->y;


				vertexCoords[triangleIndex * 9] =     positions[faces[i]->positionIndices[0]-1]->x;
				vertexCoords[triangleIndex * 9 + 1] = positions[faces[i]->positionIndices[0]-1]->y;
				vertexCoords[triangleIndex * 9 + 2] = positions[faces[i]->positionIndices[0]-1]->z;

				vertexCoords[triangleIndex * 9 + 3] = positions[faces[i]->positionIndices[1]-1]->x;
				vertexCoords[triangleIndex * 9 + 4] = positions[faces[i]->positionIndices[1]-1]->y;
				vertexCoords[triangleIndex * 9 + 5] = positions[faces[i]->positionIndices[1]-1]->z;

				vertexCoords[triangleIndex * 9 + 6] = positions[faces[i]->positionIndices[2]-1]->x;
				vertexCoords[triangleIndex * 9 + 7] = positions[faces[i]->positionIndices[2]-1]->y;
				vertexCoords[triangleIndex * 9 + 8] = positions[faces[i]->positionIndices[2]-1]->z;


				vertexNormalCoords[triangleIndex * 9] =     normals[faces[i]->normalIndices[0]-1]->x;    
				vertexNormalCoords[triangleIndex * 9 + 1] = normals[faces[i]->normalIndices[0]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 2] = normals[faces[i]->normalIndices[0]-1]->z;

				vertexNormalCoords[triangleIndex * 9 + 3] = normals[faces[i]->normalIndices[1]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 4] = normals[faces[i]->normalIndices[1]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 5] = normals[faces[i]->normalIndices[1]-1]->z;

				vertexNormalCoords[triangleIndex * 9 + 6] = normals[faces[i]->normalIndices[2]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 7] = normals[faces[i]->normalIndices[2]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 8] = normals[faces[i]->normalIndices[2]-1]->z;
								
				triangleIndex++;


				vertexTexCoords[triangleIndex * 6] =     texcoords[faces[i]->texcoordIndices[1]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 1] = 1-texcoords[faces[i]->texcoordIndices[1]-1]->y;
				
				vertexTexCoords[triangleIndex * 6 + 2] = texcoords[faces[i]->texcoordIndices[2]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 3] = 1-texcoords[faces[i]->texcoordIndices[2]-1]->y;

				vertexTexCoords[triangleIndex * 6 + 4] = texcoords[faces[i]->texcoordIndices[3]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 5] = 1-texcoords[faces[i]->texcoordIndices[3]-1]->y;


				vertexCoords[triangleIndex * 9] =     positions[faces[i]->positionIndices[1]-1]->x;
				vertexCoords[triangleIndex * 9 + 1] = positions[faces[i]->positionIndices[1]-1]->y;
				vertexCoords[triangleIndex * 9 + 2] = positions[faces[i]->positionIndices[1]-1]->z;

				vertexCoords[triangleIndex * 9 + 3] = positions[faces[i]->positionIndices[2]-1]->x;
				vertexCoords[triangleIndex * 9 + 4] = positions[faces[i]->positionIndices[2]-1]->y;
				vertexCoords[triangleIndex * 9 + 5] = positions[faces[i]->positionIndices[2]-1]->z;

				vertexCoords[triangleIndex * 9 + 6] = positions[faces[i]->positionIndices[3]-1]->x;
				vertexCoords[triangleIndex * 9 + 7] = positions[faces[i]->positionIndices[3]-1]->y;
				vertexCoords[triangleIndex * 9 + 8] = positions[faces[i]->positionIndices[3]-1]->z;


				vertexNormalCoords[triangleIndex * 9] =     normals[faces[i]->normalIndices[1]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 1] = normals[faces[i]->normalIndices[1]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 2] = normals[faces[i]->normalIndices[1]-1]->z;

				vertexNormalCoords[triangleIndex * 9 + 3] = normals[faces[i]->normalIndices[2]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 4] = normals[faces[i]->normalIndices[2]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 5] = normals[faces[i]->normalIndices[2]-1]->z;

				vertexNormalCoords[triangleIndex * 9 + 6] = normals[faces[i]->normalIndices[3]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 7] = normals[faces[i]->normalIndices[3]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 8] = normals[faces[i]->normalIndices[3]-1]->z;

				triangleIndex++;
			}
			else 
			{
				vertexTexCoords[triangleIndex * 6] =     texcoords[faces[i]->texcoordIndices[0]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 1] = 1-texcoords[faces[i]->texcoordIndices[0]-1]->y;
				
				vertexTexCoords[triangleIndex * 6 + 2] = texcoords[faces[i]->texcoordIndices[1]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 3] = 1-texcoords[faces[i]->texcoordIndices[1]-1]->y;

				vertexTexCoords[triangleIndex * 6 + 4] = texcoords[faces[i]->texcoordIndices[2]-1]->x;
				vertexTexCoords[triangleIndex * 6 + 5] = 1-texcoords[faces[i]->texcoordIndices[2]-1]->y;

				vertexCoords[triangleIndex * 9] =     positions[faces[i]->positionIndices[0]-1]->x;
				vertexCoords[triangleIndex * 9 + 1] = positions[faces[i]->positionIndices[0]-1]->y;
				vertexCoords[triangleIndex * 9 + 2] = positions[faces[i]->positionIndices[0]-1]->z;

				vertexCoords[triangleIndex * 9 + 3] = positions[faces[i]->positionIndices[1]-1]->x;
				vertexCoords[triangleIndex * 9 + 4] = positions[faces[i]->positionIndices[1]-1]->y;
				vertexCoords[triangleIndex * 9 + 5] = positions[faces[i]->positionIndices[1]-1]->z;

				vertexCoords[triangleIndex * 9 + 6] = positions[faces[i]->positionIndices[2]-1]->x;
				vertexCoords[triangleIndex * 9 + 7] = positions[faces[i]->positionIndices[2]-1]->y;
				vertexCoords[triangleIndex * 9 + 8] = positions[faces[i]->positionIndices[2]-1]->z;


				vertexNormalCoords[triangleIndex * 9] =     normals[faces[i]->normalIndices[0]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 1] = normals[faces[i]->normalIndices[0]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 2] = normals[faces[i]->normalIndices[0]-1]->z;

				vertexNormalCoords[triangleIndex * 9 + 3] = normals[faces[i]->normalIndices[1]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 4] = normals[faces[i]->normalIndices[1]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 5] = normals[faces[i]->normalIndices[1]-1]->z;

				vertexNormalCoords[triangleIndex * 9 + 6] = normals[faces[i]->normalIndices[2]-1]->x;
				vertexNormalCoords[triangleIndex * 9 + 7] = normals[faces[i]->normalIndices[2]-1]->y;
				vertexNormalCoords[triangleIndex * 9 + 8] = normals[faces[i]->normalIndices[2]-1]->z;

				triangleIndex++;
			}
		}
	}

	glGenVertexArrays(1, &vao);	
	glBindVertexArray(vao);		

	unsigned int vbo[3];		
	glGenBuffers(3, &vbo[0]);	

	glBindBuffer(GL_ARRAY_BUFFER, vbo[0]); 
	glBufferData(GL_ARRAY_BUFFER, nTriangles * 9 * sizeof(float), vertexCoords, GL_STATIC_DRAW);	    
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);     

	glBindBuffer(GL_ARRAY_BUFFER, vbo[1]); 
	glBufferData(GL_ARRAY_BUFFER, nTriangles * 6 * sizeof(float), vertexTexCoords, GL_STATIC_DRAW);	
	glEnableVertexAttribArray(1);  
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, NULL); 
	
	glBindBuffer(GL_ARRAY_BUFFER, vbo[2]); 
	glBufferData(GL_ARRAY_BUFFER, nTriangles * 9 * sizeof(float), vertexNormalCoords, GL_STATIC_DRAW);	    
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, NULL);     
	
	delete vertexCoords;
	delete vertexTexCoords;
	delete vertexNormalCoords;
}


void Mesh::DrawModel()
{
	glBindVertexArray(vao); 
	glDrawArrays(GL_TRIANGLES, 0, nTriangles * 3);	
}


Mesh::~Mesh()
{
	for(unsigned int i = 0; i < rows.size(); i++) delete rows[i];   
	for(unsigned int i = 0; i < positions.size(); i++) delete positions[i];
	for(unsigned int i = 0; i < submeshFaces.size(); i++)
		for(unsigned int j = 0; j < submeshFaces.at(i).size(); j++)
			delete submeshFaces.at(i).at(j); 
	for(unsigned int i = 0; i < normals.size(); i++) delete normals[i];
	for(unsigned int i = 0; i < texcoords.size(); i++) delete texcoords[i];
}


///
class Material
{
    vec3 ka, kd, ks;
    float shininess;
    
public:
    Material(vec3 ka, vec3 kd, vec3 ks, float shininess) : ka(ka), kd(kd), ks(ks), shininess(shininess) {}
    
    void SetMaterial(unsigned int shader)
    {
        
        int location = glGetUniformLocation(shader, "ka");
        if(location >= 0) glUniform3fv(location, 1, &(ka.x));
        else printf("uniform ka cant be set \n");
        
        
        location = glGetUniformLocation(shader, "kd");
        if(location >= 0) glUniform3fv(location, 1, &(kd.x));
        else printf("uniform kd cant be set \n");
        
        
        location = glGetUniformLocation(shader, "ks");
        if(location >= 0) glUniform3fv(location, 1, &(ks.x));
        else printf("uniform ks cant be set \n");
        
        
        location = glGetUniformLocation(shader, "shininess");
        if(location >= 0) glUniform1f(location, shininess);
        else printf("uniform shininess cant be set \n");
    }
    
    
};

class Light
{
vec3 La, Le;
    vec4 worldLightPosition;
    
    
public:
    Light(vec3 La, vec3 Le, vec4 worldLightPosition) : La(La), Le(Le), worldLightPosition(worldLightPosition){
        
    }
    
    
    void SetLightPosition(unsigned int shader){
        int location = glGetUniformLocation(shader, "worldLightPosition");
        if(location >= 0) glUniform4fv(location, 1, &worldLightPosition.v[0]);
        else printf("uniform worldLightPositionn cant be set \n");
    }
    
    void SetLight(unsigned int shader)
    {
        
        int location = glGetUniformLocation(shader, "La");
        if(location >= 0) glUniform3fv(location, 1, &(La.x));
        else printf("uniform La cant be set \n");
        
        
        location = glGetUniformLocation(shader, "Le");
        if(location >= 0) glUniform3fv(location, 1, &(Le.x));
        else printf("uniform Le cant be set \n");
        
        location = glGetUniformLocation(shader, "worldLightPosition");
        if(location >= 0) glUniform4fv(location, 1, &(worldLightPosition.v[0]));
        else printf("uniform worldLightPosition cant be set \n");
        
    }
    
    void SetPointLightSource(vec3 pos)
    {
        worldLightPosition.v[0] = pos.x;
        worldLightPosition.v[1] = pos.y;
        worldLightPosition.v[2] = pos.z;
        worldLightPosition.v[3] = 1;
    }
    
    
    void SetDirectionalLightSource(vec3 pos)
    {
        worldLightPosition.v[0] = pos.x;
        worldLightPosition.v[1] = pos.y;
        worldLightPosition.v[2] = pos.z;
        worldLightPosition.v[3] = 0;
    }
    
};

Light light(vec3(1, 1, 1), vec3(1, 1, 1), vec4(0.1, 0.1, 0.1, 0.0));


class Camera {
   vec3  wEye, wLookat, wVup, velocity, position;
   float fov, asp, fp, bp, angularVelocity, orientation;

public:
	Camera()
	{
		wEye = vec3(0.0, 0.0, 1.0);
		wLookat = vec3(0.0, 0.0, 0.0);
		wVup = vec3(0.0, 1.0, 0.0);
		fov = M_PI / 2.0; asp = 1.0; fp = 0.01; bp = 10.0;
        light.SetPointLightSource(wEye);
	}
	
    void SetAspectRatio(float a) { asp = a; }
    
    void SetEyePosition(unsigned int shader) {
        int location = glGetUniformLocation(shader, "wEye");
        if(location >= 0) glUniform3fv(location, 1, &(wEye.x));
        else printf("uniform wEye cant be set \n");

    
    }
    
    vec3 GetEyePosition() {
        return wEye;
        
    }
    
    void Control(){
        velocity = vec3(0,0,0);
        angularVelocity = 0;
        float alpha = M_PI / 18000.0;
        
        
        if(keyPressed['w']){
            velocity.z -= 0.008;
            light.SetPointLightSource(wEye);
            
        }
        if(keyPressed['s']){
            velocity.z += 0.008;
            light.SetPointLightSource(wEye);
        }
        
        if(keyPressed['z']){
            velocity.x -= 0.008;
            light.SetPointLightSource(wEye);
            
        }
        if(keyPressed['c']){
            velocity.x += 0.008;
            light.SetPointLightSource(wEye);
        }
        
        if(keyPressed['d']){
            vec3 ahead = wLookat - wEye;
            float length = ahead.length();
            vec3 right = cross(ahead, wVup);
            ahead = ahead.normalize()*cos(alpha) + right * sin(alpha);
            wLookat = wEye + ahead * length;
            angularVelocity += .060;
            light.SetPointLightSource(wEye);
        }
        if(keyPressed['a']){
            vec3 ahead = wLookat - wEye;
            float length = ahead.length();
            vec3 right = cross(ahead, wVup);
            ahead = ahead.normalize()*cos(alpha) - right * sin(alpha);
            wLookat = wEye + ahead * length;
            angularVelocity += -.060;
            light.SetPointLightSource(wEye);
        }
    }
    
    void Move(){
        wEye += velocity;
        wLookat += velocity;
        wLookat.x += angularVelocity;
    }

	mat4 GetViewMatrix() // view matrix 
	{ 
		vec3 w = (wEye - wLookat).normalize();
		vec3 u = cross(wVup, w).normalize();
		vec3 v = cross(w, u);
	
		return  
			mat4(	
				1.0f,    0.0f,    0.0f,    0.0f,
				0.0f,    1.0f,    0.0f,    0.0f,
				0.0f,    0.0f,    1.0f,    0.0f,
				-wEye.x, -wEye.y, -wEye.z, 1.0f ) *
			mat4(	
				u.x,  v.x,  w.x,  0.0f,
				u.y,  v.y,  w.y,  0.0f,
				u.z,  v.z,  w.z,  0.0f,
				0.0f, 0.0f, 0.0f, 1.0f );
   }
    
    
    mat4 GetInverseViewMatrix() {
        vec3 w = (wEye - wLookat).normalize();
        vec3 u = cross(wVup, w).normalize();
        vec3 v = cross(w, u);
        return mat4(
                    u.x,  u.y,  u.z,  0.0f,
                    v.x,  v.y,  v.z,  0.0f,
                    w.x,  w.y,  w.z,  0.0f,
                    0.0f, 0.0f, 0.0f, 1.0f );
    }
    
    mat4 GetProjectionMatrix() // projection matrix
    {
        float sy = 1/tan(fov/2);
        return mat4(
                    sy/asp, 0.0f,  0.0f,               0.0f,
                    0.0f,   sy,    0.0f,               0.0f,
                    0.0f,   0.0f, -(fp+bp)/(bp - fp), -1.0f,
                    0.0f,   0.0f, -2*fp*bp/(bp - fp),  0.0f);
    }
    
    mat4 GetInverseProjectionMatrix() {
        float sy = 1/tan(fov/2);
        return mat4(
                    asp/sy,		0.0f,		0.0f,		0.0f,
                    0.0f,		1.0/sy,		0.0f,		0.0f,
                    0.0f,		0.0f,		0.0f,		(fp - bp) / 2.0 / fp / bp,
                    0.0f,		0.0f,		-1.0f,		(fp + bp) / 2.0 / fp / bp);
    }
};

Camera camera;


unsigned int shaderProgram0;
unsigned int shaderProgram1;
unsigned int shaderProgram2;
unsigned int shaderProgram3;

extern "C" unsigned char* stbi_load(char const *filename, int *x, int *y, int *comp, int req_comp);


class TextureCube
{
    unsigned int textureId;
public:
    TextureCube(
                const std::string& inputFileName0, const std::string& inputFileName1, const std::string& inputFileName2,
                const std::string& inputFileName3, const std::string& inputFileName4, const std::string& inputFileName5)
    {
        unsigned char* data[6]; int width[6]; int height[6]; int nComponents[6]; std::string filename[6];
        filename[0] = inputFileName0; filename[1] = inputFileName1; filename[2] = inputFileName2;
        filename[3] = inputFileName3; filename[4] = inputFileName4; filename[5] = inputFileName5;
        for(int i = 0; i < 6; i++) {
            data[i] = stbi_load(filename[i].c_str(), &width[i], &height[i], &nComponents[i], 0);
            if(data == NULL) return;
        }
        glGenTextures(1, &textureId); glBindTexture(GL_TEXTURE_CUBE_MAP, textureId);
        for(int i = 0; i < 6; i++) {
            if(nComponents[i] == 4) glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGBA,
                                                 width[i], height[i], 0, GL_RGBA, GL_UNSIGNED_BYTE, data[i]);
            if(nComponents[i] == 3) glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB,
                                                 width[i], height[i], 0, GL_RGB, GL_UNSIGNED_BYTE, data[i]);
        }
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        for(int i = 0; i < 6; i++)  delete data[i];
    }
    
    
    void Bind(unsigned int shader)
    {
        int samplerCube = 1;
        int location = glGetUniformLocation(shader, "environmentMap");
        glUniform1i(location, samplerCube);
        glActiveTexture(GL_TEXTURE0 + samplerCube);
        glBindTexture(GL_TEXTURE_CUBE_MAP, textureId);
    }
};

TextureCube *environmentMap = 0;




class Object {
protected:
    vec3 position, scaling;
    float orientation;

	vec3 velocity, acceleration;
	float angularVelocity;

	unsigned int shader;
    unsigned int shaderShadow;
    
    
    Material *material;

public:
	Object(unsigned int sp, Material *m) : scaling(1.0, 1.0, 1.0), position(0.0, 0.0, 0.0), orientation(0.0), angularVelocity(0.0), shader(sp), material(m) { }

    
	virtual void SetTransform(unsigned int sp)
	{
		mat4 scale(
			scaling.x, 0, 0, 0,
			0, scaling.y, 0, 0,
			0, 0, scaling.z, 0,
			0, 0, 0, 1); 

		mat4 scaleInv(
			1.0 / scaling.x, 0, 0, 0,
			0, 1.0 / scaling.y, 0, 0,
			0, 0, 1.0 / scaling.z, 0,
			0, 0, 0, 1); 

		float alpha = orientation / 180 * M_PI;

		mat4 rotate(
			cos(alpha), 0, sin(alpha), 0,
			0, 1, 0, 0,
			-sin(alpha), 0, cos(alpha), 0,
			0, 0, 0, 1); 

        
        ///
        //wheel rotate mat
        // M1M, frst apply wheel  matrix then general matrix (apply sequence of modeling trans to parts of car)
        // I CHANGED THE GENERAL ROTATION TO BE AROUND THE X AXIS INSTEAD OF Y (original version)
        
        float alphaRot = orientation / 180 * M_PI;
        
        mat4 rotateWheel(
                    1, 0, 0, 0,
                    0, cos(alphaRot), sin(alphaRot), 0,
                    0, -sin(alphaRot), cos(alphaRot), 0,
                    0, 0, 0, 1);
  
		mat4 rotateInv(
			cos(alpha), 0, -sin(alpha), 0,
			0, 1, 0, 0,
			sin(alpha), 0, cos(alpha), 0,
			0, 0, 0, 1);
        
        
        mat4 rotateWheelInv(
                    1, 0, 0, 0,
                    0, cos(alpha), -sin(alpha), 0,
                    0, sin(alpha), cos(alpha), 0,
                    0, 0, 0, 1);

		mat4 translate(
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			position.x, position.y, position.z, 1); 
				
		mat4 translateInv(
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			-position.x, -position.y, -position.z, 1); 
		
		mat4 M = scale * rotate * translate;
        
          
        
        int location = glGetUniformLocation(shader, "M");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, M);
        else printf("uniform M cannot be set\n");
        
	mat4 MInv = translateInv * rotateInv * scaleInv;

	location = glGetUniformLocation(shader, "MInv");
	if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, MInv); 
	else printf("uniform MInv cannot be set\n");

        
	mat4 MVP = M * camera.GetViewMatrix() * camera.GetProjectionMatrix();



	location = glGetUniformLocation(shader, "MVP");
	if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, MVP); 
	else printf("uniform MVP cannot be set\n");
        
    
	}
    
    mat4 GetModelM()
    {
        mat4 scale(
                   scaling.x, 0, 0, 0,
                   0, scaling.y, 0, 0,
                   0, 0, scaling.z, 0,
                   0, 0, 0, 1);
        
        mat4 scaleInv(
                      1.0 / scaling.x, 0, 0, 0,
                      0, 1.0 / scaling.y, 0, 0,
                      0, 0, 1.0 / scaling.z, 0,
                      0, 0, 0, 1);
        
        float alpha = orientation / 180 * M_PI;
        
        mat4 rotate(
                    cos(alpha), 0, sin(alpha), 0,
                    0, 1, 0, 0,
                    -sin(alpha), 0, cos(alpha), 0,
                    0, 0, 0, 1);
        
        
        mat4 rotateInv(
                       cos(alpha), 0, -sin(alpha), 0,
                       0, 1, 0, 0,
                       sin(alpha), 0, cos(alpha), 0,
                       0, 0, 0, 1);
        
        
        mat4 translate(
                       1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, 1, 0,
                       position.x, position.y, position.z, 1); 
        
        mat4 translateInv(
                          1, 0, 0, 0,
                          0, 1, 0, 0,
                          0, 0, 1, 0,
                          -position.x, -position.y, -position.z, 1); 
        
        
        mat4 M = scale * rotate * translate;
        return M;
        
    }
    
    
    mat4 GetModelMInv()
    {
        mat4 scale(
                   scaling.x, 0, 0, 0,
                   0, scaling.y, 0, 0,
                   0, 0, scaling.z, 0,
                   0, 0, 0, 1);
        
        mat4 scaleInv(
                      1.0 / scaling.x, 0, 0, 0,
                      0, 1.0 / scaling.y, 0, 0,
                      0, 0, 1.0 / scaling.z, 0,
                      0, 0, 0, 1);
        
        float alpha = orientation / 180 * M_PI;
        
        mat4 rotate(
                    cos(alpha), 0, sin(alpha), 0,
                    0, 1, 0, 0,
                    -sin(alpha), 0, cos(alpha), 0,
                    0, 0, 0, 1);
        
        
        mat4 rotateInv(
                       cos(alpha), 0, -sin(alpha), 0,
                       0, 1, 0, 0,
                       sin(alpha), 0, cos(alpha), 0,
                       0, 0, 0, 1);
        
        
        mat4 translate(
                       1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, 1, 0,
                       position.x, position.y, position.z, 1);
        
        mat4 translateInv(
                          1, 0, 0, 0,
                          0, 1, 0, 0,
                          0, 0, 1, 0,
                          -position.x, -position.y, -position.z, 1);
        
        mat4 MInv = translateInv * rotateInv * scaleInv;
        return MInv;
        
        
        
    }

    void SetTransform2(unsigned int sp){
        mat4 scale(
                   scaling.x, 0, 0, 0,
                   0, scaling.y, 0, 0,
                   0, 0, scaling.z, 0,
                   0, 0, 0, 1);
        
        
        float alpha = orientation / 180 * M_PI;
        
        mat4 rotate(
                    cos(alpha), 0, sin(alpha), 0,
                    0, 1, 0, 0,
                    -sin(alpha), 0, cos(alpha), 0,
                    0, 0, 0, 1);
        
        
        
        
        mat4 translate(
                       1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, 1, 0,
                       position.x, position.y, position.z, 1);
        
        
        mat4 M = scale * rotate * translate;
        
        
        
        int location = glGetUniformLocation(sp, "M");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, M);
        else printf("uniform M cannottt be set\n");
        
        
        mat4 MVP = camera.GetViewMatrix() * camera.GetProjectionMatrix();
        
        
        
        location = glGetUniformLocation(sp, "MVP");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, MVP); 
        else printf("uniform MVP cannottt be set\n");
        
    }

    
	virtual void Draw()
	{
		glUseProgram(shader);
        if(environmentMap) environmentMap->Bind(shader);
		SetTransform(shader);
        light.SetPointLightSource(camera.GetEyePosition());
        light.SetLight(shader);
        camera.SetEyePosition(shader);
        material->SetMaterial(shader);
		DrawModel();
	}
    
    void setPos(vec3 newPos)
    {
        position = newPos;
    }
    
    void setPosX(int newX)
    {
        position.x = newX;
    }
    
    void setPosY(int newY)
    {
        position.y = newY;
    }
    
    void setPosZ(int newZ)
    {
        position.z = newZ;
    }
    
    vec3 getPos()
    {
        return position;
    }
    
    
    virtual void DrawShadow()
    {
        glUseProgram(shaderProgram2);
        SetTransform2(shaderProgram2);
        light.SetPointLightSource(vec3(0.0, 10.0, 3.0));
        light.SetLightPosition(shaderProgram2);
        DrawModel();
    }
    
    virtual void DrawModel() = 0;
    
    Object* Scale(const vec3& s) { scaling *= s; return this; }
    
    Object* Rotate(float angle) { orientation += angle; return this; }
    
    Object* Translate(const vec3& t) { position += t; return this; }


	virtual void Move(float dt)
	{
		position = position + velocity * dt;
        
        velocity = velocity + acceleration * dt;

		orientation = orientation + angularVelocity * dt;
	}

	virtual void Control()
	{
        
	}
};



class Texture
{
	unsigned int textureId;

public:
	Texture(const std::string& inputFileName)
	{
		unsigned char* data;
		int width; int height; int nComponents = 4;

		data = stbi_load(inputFileName.c_str(), &width, &height, &nComponents, 0);

		if(data == NULL) 
		{ 
			return;
		}

		glGenTextures(1, &textureId); 
		glBindTexture(GL_TEXTURE_2D, textureId); 
		
		if(nComponents == 4) glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
		if(nComponents == 3) glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 

		delete data; 
	}

	void Bind(unsigned int shader)
	{
		int samplerUnit = 0; 
		int location = glGetUniformLocation(shader, "samplerUnit");
		glUniform1i(location, samplerUnit);
		glActiveTexture(GL_TEXTURE0 + samplerUnit); 
		glBindTexture(GL_TEXTURE_2D, textureId);
	}
};



class TexturedQuad : public Object {
    
    Texture *texture;
    unsigned int vao;
    
public:
    TexturedQuad(Material *m, Texture* t, unsigned int sp = shaderProgram1) : Object(sp, m), texture(t)
    {
        // similar to the default constructor of Quad, but vertex attribute array 1 has to
        // store the texture coordinates rather than the color attributes
        glGenVertexArrays(1, &vao);// create 1 vertex array object
        glBindVertexArray(vao);// make it active
        
        unsigned int vbo[3];// vertex buffer objects
        glGenBuffers(3, &vbo[0]);// generate 2 vertex buffer objects
        
        // vertex coordinates: vbo[0] -> Attrib Array 0 -> vertexPosition of the vertex shader
        glBindBuffer(GL_ARRAY_BUFFER, vbo[0]); // make it active, it is an array
        static float vertexCoords[] =  {-10, -0.7, -10, -10, -0.7, 10, 10, -0.7, -5, 10, -0.7, 10};// vertex data on the CPU
                
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertexCoords), vertexCoords, GL_STATIC_DRAW);
        // map Attribute Array 0 to the current bound vertex buffer (vbo[0])
        glEnableVertexAttribArray(0);
        // data organization of Attribute Array 0
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
        
        // vertex colors: vbo[1] -> Attrib Array 1 -> vertexColor of the vertex shader
        glBindBuffer(GL_ARRAY_BUFFER, vbo[1]); // make it active, it is an array
        static float texCoord[] = {0, 0, 5, 0, 0, 5, 5, 5};// texture coordinates
        glBufferData(GL_ARRAY_BUFFER, sizeof(texCoord), texCoord, GL_STATIC_DRAW);
            glEnableVertexAttribArray(1);  // Vertex position
        // data organization of Attribute Array 1
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, NULL);
        
        glBindBuffer(GL_ARRAY_BUFFER, vbo[2]); // make it active, it is an array
        static float vertexNormal[] = {0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0};// texture coordinates
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertexNormal), vertexNormal, GL_STATIC_DRAW);      glEnableVertexAttribArray(2);  // Vertex position
        // data organization of Attribute Array 1
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    }
    
    void DrawModel()
    {
        texture->Bind(shader);

        glEnable(GL_DEPTH_TEST);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);
        glDisable(GL_DEPTH_TEST);
    }
    
    
    
};


class FullQuad : public Object {
    
    Texture *texture;
    unsigned int vao;
    
public:
    FullQuad(Material *m, Texture* t, unsigned int sp = shaderProgram3) : Object(sp, m), texture(t)
    {
        // similar to the default constructor of Quad, but vertex attribute array 1 has to
        // store the texture coordinates rather than the color attributes
        glGenVertexArrays(1, &vao);// create 1 vertex array object
        glBindVertexArray(vao);// make it active
        
        unsigned int vbo;// vertex buffer objects
        glGenBuffers(1, &vbo);// generate 2 vertex buffer objects
        
        // vertex coordinates: vbo[0] -> Attrib Array 0 -> vertexPosition of the vertex shader
        glBindBuffer(GL_ARRAY_BUFFER, vbo); // make it active, it is an array
        
        static float vertexCoords[] = {-1.0, -1.0, 0.0, 1.0,
                                        1.0, -1.0, 0.0, 1.0,
                                        -1.0, 1.0, 0.0, 1.0,
                                        1.0, 1.0, 0.0, 1.0};
        
        
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertexCoords), vertexCoords, GL_STATIC_DRAW);
        // map Attribute Array 0 to the current bound vertex buffer (vbo[0])
        glEnableVertexAttribArray(0);
        // data organization of Attribute Array 0
        glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, NULL);

    }
    
    void DrawModel()
    {
        texture->Bind(shader);
        glEnable(GL_DEPTH_TEST);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);
        glDisable(GL_DEPTH_TEST);
        
    }
    
    void Draw(){
        glUseProgram(shader);
        if (environmentMap) environmentMap->Bind(shader);
        SetTransform(shader);
        DrawModel();
    }
    
    void SetTransform(unsigned int sp)
    {
        
        mat4 viewDirMatrix =
        camera.GetInverseProjectionMatrix() * camera.GetInverseViewMatrix();
        
        int location = glGetUniformLocation(sp, "viewDirMatrix");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, viewDirMatrix);
        else printf("uniform viewDirMatrix cannot be set\n");

    }
    
    void DrawShadow() {}
    
};



class MeshInstance : public Object
{
    Texture *texture;
	Mesh * mesh;

public:

	MeshInstance(Texture* t, Mesh* m, Material* mt, unsigned int sp = shaderProgram0) : Object(sp, mt), texture(t), mesh(m)
	{
		position = vec3(0.0, -0.9, 0.0);
	}

	virtual void DrawModel()
	{
		texture->Bind(shader);
	
		glEnable(GL_DEPTH_TEST);
		mesh->DrawModel();
		glDisable(GL_DEPTH_TEST);
	}
};

class Tigger : public MeshInstance
{
    Texture *texture;
    Mesh * mesh;
    vec3  wEye, wLookat, wVup;
    float fov, asp, fp, bp;
    
public:
    
    Tigger(Texture* t, Mesh* m, Material* mt, unsigned int sp = shaderProgram0) : MeshInstance(t, m, mt, sp), texture(t), mesh(m)
    {
        scaling = vec3(0.03, 0.03, 0.03);
        position = vec3(0.8, -0.9, -0.4);
        angularVelocity = 0;
        velocity = vec3(0, 0, 0);
        
    }
    
    void DrawModel()
    {
        texture->Bind(shader);
        
        glEnable(GL_DEPTH_TEST);
        mesh->DrawModel();
        glDisable(GL_DEPTH_TEST);
    }
    
    void Control(){
        acceleration = vec3(0, 0.0, 0);
        angularVelocity = 0.0;
        
        if(keyPressed[';']){
            acceleration += vec3(0.8, 0, 0);
        }
        if(keyPressed['k']){
            acceleration += vec3(-0.8, 0, 0);
        }
        if(keyPressed['i']){
            angularVelocity +=  100;
        }
        if(keyPressed['p']){
            angularVelocity -=  100;
        }
        
        if(keyPressed['o']){
            acceleration += vec3(0, 0, -0.6);
        }
        if(keyPressed['l']){
            acceleration += vec3(0.0, 0, 0.6);
        }
    }
    
};

Tigger *tigger;


class Car : public MeshInstance
{
    Texture *texture;
    Mesh * mesh;
    vec3  wEye, wLookat, wVup;
    float fov, asp, fp, bp;
    
public:
    
    Car(Texture* t, Mesh* m, Material* mt, unsigned int sp = shaderProgram0) : MeshInstance(t, m, mt, sp), texture(t), mesh(m)
    {
        scaling = vec3(0.05, 0.05, 0.05);
        position = vec3(2, -0.6, -4);
        angularVelocity = 0;
        velocity = vec3(0, 0, 0);
    }
    
    void DrawModel()
    {
        texture->Bind(shader);
        
        glEnable(GL_DEPTH_TEST);
        mesh->DrawModel();
        glDisable(GL_DEPTH_TEST);
    }
    
};


Car *car;

class Wheel : public MeshInstance
{
    Texture *texture;
    Mesh * mesh;
    
public:
    
    Wheel(Texture* t, Mesh* m, Material* mt, unsigned int sp = shaderProgram0) : MeshInstance(t, m, mt, sp), texture(t), mesh(m)
    {
        scaling = vec3(0.05, 0.05, 0.05);
        position = vec3(0.0, -0.9, 0.0);
        angularVelocity = 580;
        
    }
    
    void DrawModel()
    {
        texture->Bind(shader);
        
        glEnable(GL_DEPTH_TEST);
        mesh->DrawModel();
        glDisable(GL_DEPTH_TEST);
    }
    
    
    
    void SetTransform(unsigned int sp)
    {
        mat4 scale(
                   scaling.x, 0, 0, 0,
                   0, scaling.y, 0, 0,
                   0, 0, scaling.z, 0,
                   0, 0, 0, 1);
        
        mat4 scaleInv(
                      1.0 / scaling.x, 0, 0, 0,
                      0, 1.0 / scaling.y, 0, 0,
                      0, 0, 1.0 / scaling.z, 0,
                      0, 0, 0, 1);
        
        float alpha = orientation / 180 * M_PI;
        
        mat4 rotate(
                    cos(alpha), 0, sin(alpha), 0,
                    0, 1, 0, 0,
                    -sin(alpha), 0, cos(alpha), 0,
                    0, 0, 0, 1);

        
        float alphaRot = orientation / 180 * M_PI;
        
        mat4 rotateWheel(
                         1, 0, 0, 0,
                         0, cos(alphaRot), sin(alphaRot), 0,
                         0, -sin(alphaRot), cos(alphaRot), 0,
                         0, 0, 0, 1);
     
        mat4 rotateInv(
                       cos(alpha), 0, -sin(alpha), 0,
                       0, 1, 0, 0,
                       sin(alpha), 0, cos(alpha), 0,
                       0, 0, 0, 1);
        
        
        mat4 rotateWheelInv(
                            1, 0, 0, 0,
                            0, cos(alpha), -sin(alpha), 0,
                            0, sin(alpha), cos(alpha), 0,
                            0, 0, 0, 1);
        
        mat4 translate(
                       1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, 1, 0,
                       position.x, position.y, position.z, 1);
        
        mat4 translateInv(
                          1, 0, 0, 0,
                          0, 1, 0, 0,
                          0, 0, 1, 0,
                          -position.x, -position.y, -position.z, 1);
        
        mat4 M = scale * rotate * translate;
        
        
        mat4 MWheel = scale * rotateWheel * translate;
        
        
        int location = glGetUniformLocation(sp, "M");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, M);
        else printf("uniform M cannot be set\n");
        
        mat4 MInv = translateInv * rotateWheelInv * scaleInv;
        
        location = glGetUniformLocation(sp, "MInv");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, MInv);
        else printf("uniform MInv cannot be set\n");
        
        
        mat4 MVP = MWheel * camera.GetViewMatrix() * camera.GetProjectionMatrix();
        
        
        location = glGetUniformLocation(sp, "MVP");
        if (location >= 0) glUniformMatrix4fv(location, 1, GL_TRUE, MVP);
        else printf("uniform MVP cannot be set\n");
        
    }
    
    mat4 GetModelM()
    {
        mat4 scale(
                   scaling.x, 0, 0, 0,
                   0, scaling.y, 0, 0,
                   0, 0, scaling.z, 0,
                   0, 0, 0, 1);
        
        mat4 scaleInv(
                      1.0 / scaling.x, 0, 0, 0,
                      0, 1.0 / scaling.y, 0, 0,
                      0, 0, 1.0 / scaling.z, 0,
                      0, 0, 0, 1);
        
        
        float alphaRot = orientation / 180 * M_PI;
        
        mat4 rotateWheel(
                         1, 0, 0, 0,
                         0, cos(alphaRot), sin(alphaRot), 0,
                         0, -sin(alphaRot), cos(alphaRot), 0,
                         0, 0, 0, 1);
        
        mat4 rotateWheelInv(
                            1, 0, 0, 0,
                            0, cos(alphaRot), -sin(alphaRot), 0,
                            0, sin(alphaRot), cos(alphaRot), 0,
                            0, 0, 0, 1);
        
        mat4 translate(
                       1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, 1, 0,
                       position.x, position.y, position.z, 1);
        
        mat4 translateInv(
                          1, 0, 0, 0,
                          0, 1, 0, 0,
                          0, 0, 1, 0,
                          -position.x, -position.y, -position.z, 1);
        
        
        mat4 M = scale * rotateWheel * translate;
        return M;
        
    }
    
    
    mat4 GetModelMInv()
    {
        mat4 scale(
                   scaling.x, 0, 0, 0,
                   0, scaling.y, 0, 0,
                   0, 0, scaling.z, 0,
                   0, 0, 0, 1);
        
        mat4 scaleInv(
                      1.0 / scaling.x, 0, 0, 0,
                      0, 1.0 / scaling.y, 0, 0,
                      0, 0, 1.0 / scaling.z, 0,
                      0, 0, 0, 1);
        
        
        float alphaRot = orientation / 180 * M_PI;
        
        mat4 rotateWheel(
                         1, 0, 0, 0,
                         0, cos(alphaRot), sin(alphaRot), 0,
                         0, -sin(alphaRot), cos(alphaRot), 0,
                         0, 0, 0, 1);
        
        mat4 rotateWheelInv(
                            1, 0, 0, 0,
                            0, cos(alphaRot), -sin(alphaRot), 0,
                            0, sin(alphaRot), cos(alphaRot), 0,
                            0, 0, 0, 1);
        
        mat4 translate(
                       1, 0, 0, 0,
                       0, 1, 0, 0,
                       0, 0, 1, 0,
                       position.x, position.y, position.z, 1);
        
        mat4 translateInv(
                          1, 0, 0, 0,
                          0, 1, 0, 0,
                          0, 0, 1, 0,
                          -position.x, -position.y, -position.z, 1);
        
        mat4 MInv = translateInv * rotateWheelInv * scaleInv;
        return MInv;
        
    }
    
};


class Tree : public MeshInstance
{
    Texture *texture;
    Mesh * mesh;
    
public:
    
    Tree(Texture* t, Mesh* m,  Material* mt, unsigned int sp = shaderProgram0) : MeshInstance(t, m, mt, sp), texture(t), mesh(m)
    {
        scaling = vec3(0.01, 0.01, 0.01);
        position = vec3(-2, -0.6, 0);
        velocity = vec3(0, 0, 0);
    }
    
    
    
    void DrawModel()
    {
        texture->Bind(shader);
        
        glEnable(GL_DEPTH_TEST);
        mesh->DrawModel();
        glDisable(GL_DEPTH_TEST);
    }
    
};

class Ground : public TexturedQuad
{
    
public:
    Ground(Material *m, Texture* t, unsigned int sp = shaderProgram1) : TexturedQuad(m, t, sp)
    {
        scaling = vec3(0.1, 0.1, 0.1);
        position = vec3(-0.5, 0.5, -0.9);

    }

    
    
};


Tree *tree;
Wheel *wheelFL;
Wheel *wheelFR;
Wheel *wheelBL;
Wheel *wheelBR;

class Scene
{
	std::vector<Texture*> textures;
	std::vector<Mesh*> meshes;
    std::vector<Object*> objects;
    std::vector<Material*> materials;

public:
	Scene() 
	{ 
		
	}

	void Initialize()
	{
        environmentMap = new TextureCube("img/posx.jpg", "img/negx.jpg", "img/posy.jpg", "img/negy.jpg", "img/posz.jpg", "img/negz.jpg");
        
        textures.push_back(new Texture("img/chevy.png"));
        textures.push_back(new Texture("img/tree.png"));
        
        textures.push_back(new Texture("img/heliait.png"));
        textures.push_back(new Texture("img/tigger.png"));

        meshes.push_back(new Mesh("img/chassis.obj"));
        meshes.push_back(new Mesh("img/tree.obj"));
        meshes.push_back(new Mesh("img/wheel.obj"));
        meshes.push_back(new Mesh("img/tigger.obj"));

        
        materials.push_back(new Material(vec3(0.1, 0.1, 0.1), vec3(0.9, 0.9, 0.9), vec3(0.0, 0.0, 0.0), 0.0));
        materials.push_back(new Material(vec3(0.1, 0.1, 0.1), vec3(0.6, 0.6, 0.6), vec3(0.3, 0.3, 0.3), 200));
  
	objects.push_back(car = new Car(textures[0], meshes[0], materials[1]));
        objects.push_back(tree = new Tree(textures[1], meshes[1], materials[0]));
        objects.push_back(tigger = new Tigger(textures[3], meshes[3], materials[0]));
        objects.push_back(new Ground(materials[1], textures[2]));
        objects.push_back(wheelFL = new Wheel(textures[0], meshes[2], materials[1]));
        objects.push_back(wheelFR = new Wheel(textures[0], meshes[2], materials[1]));
        objects.push_back(wheelBL = new Wheel(textures[0], meshes[2], materials[1]));
        objects.push_back(wheelBR = new Wheel(textures[0], meshes[2], materials[1]));
        objects.push_back(new FullQuad(materials[1], textures[2]));
        
    }

	~Scene()
	{
	for(int i = 0; i < textures.size(); i++) delete textures[i];
	for(int i = 0; i < meshes.size(); i++) delete meshes[i];
	for(int i = 0; i < objects.size(); i++) delete objects[i];
        delete environmentMap;
	}

	void Draw()
	{
        
        for(int i = 0; i < objects.size(); i++) objects[i]->DrawShadow();
        
        int carX = car->getPos().x;
        int carY = car->getPos().y;
        int carZ = car->getPos().z;
        
        wheelFL->setPos(vec3(carX-0.3, carY-0.8, carZ+0.7));
        wheelFR->setPos(vec3(carX+0.3, carY-0.8, carZ+0.7));
        wheelBL->setPos(vec3(carX-0.3, carY-0.8, carZ-0.55));
        wheelBR->setPos(vec3(carX+0.3, carY-0.8, carZ-0.55));
        
        
	for(int i = 0; i < objects.size(); i++) objects[i]->Draw();
        
	}
    
    	void DrawShadow(){
        	for(int i = 0; i < objects.size(); i++) objects[i]->DrawShadow();
    	}

	void Move(float dt)
	{
		for(int i = 0; i < objects.size(); i++) objects[i]->Move(dt);
        camera.Move();
	}

	void Control()
	{
		for(int i = 0; i < objects.size(); i++) objects[i]->Control();
        camera.Control();
	}
    
};


Scene scene;


void onInitialization() {
	glViewport(0, 0, windowWidth, windowHeight);

	// Create vertex shader from string
	unsigned int vertexShader0 = glCreateShader(GL_VERTEX_SHADER);
	if (!vertexShader0) { printf("Error in vertex shader 0 creation\n"); exit(1); }

	glShaderSource(vertexShader0, 1, &vertexSource0, NULL);
	glCompileShader(vertexShader0);
	checkShader(vertexShader0, "Vertex shader 0 error");

	// Create fragment shader from string
	unsigned int fragmentShader0 = glCreateShader(GL_FRAGMENT_SHADER);
	if (!fragmentShader0) { printf("Error in fragment shader 0 creation\n"); exit(1); }

	glShaderSource(fragmentShader0, 1, &fragmentSource0, NULL);
	glCompileShader(fragmentShader0);
	checkShader(fragmentShader0, "Fragment shader 0 error");

	// Attach shaders to a single program
	shaderProgram0 = glCreateProgram();
	if (!shaderProgram0) { printf("Error in shader program 0 creation\n"); exit(1); }

	glAttachShader(shaderProgram0, vertexShader0);
	glAttachShader(shaderProgram0, fragmentShader0);

	// Connect Attrib Arrays to input variables of the vertex shader
	glBindAttribLocation(shaderProgram0, 0, "vertexPosition");	// vertexPosition gets values from Attrib Array 0
	glBindAttribLocation(shaderProgram0, 1, "vertexTexCoord");  // vertexTexCoord gets values from Attrib Array 1
	glBindAttribLocation(shaderProgram0, 2, "vertexNormal");  // vertexNormalCoord gets values from Attrib Array 2

	// Connect the fragmentColor to the frame buffer memory
	glBindFragDataLocation(shaderProgram0, 0, "fragmentColor");	// fragmentColor goes to the frame buffer memory

	// Program packaging
	glLinkProgram(shaderProgram0);
	checkLinking(shaderProgram0);
    
    
    ///
    // Create vertex shader from string
    unsigned int vertexShader1 = glCreateShader(GL_VERTEX_SHADER);
    if (!vertexShader1) { printf("Error in vertex shader 1 creation\n"); exit(1); }
    
    glShaderSource(vertexShader1, 1, &vertexSource1, NULL);
    glCompileShader(vertexShader1);
    checkShader(vertexShader1, "Vertex shader 1 error");
    
    // Create fragment shader from string
    unsigned int fragmentShader1 = glCreateShader(GL_FRAGMENT_SHADER);
    if (!fragmentShader1) { printf("Error in fragment shader 1 creation\n"); exit(1); }
    
    glShaderSource(fragmentShader1, 1, &fragmentSource1, NULL);
    glCompileShader(fragmentShader1);
    checkShader(fragmentShader1, "Fragment shader 1 error");
    
    // Attach shaders to a single program
    shaderProgram1 = glCreateProgram();
    if (!shaderProgram1) { printf("Error in shader program 1 creation\n"); exit(1); }
    
    glAttachShader(shaderProgram1, vertexShader1);
    glAttachShader(shaderProgram1, fragmentShader1);
    
    // Connect Attrib Arrays to input variables of the vertex shader
    glBindAttribLocation(shaderProgram1, 0, "vertexPosition");	// vertexPosition gets values from Attrib Array 0
    glBindAttribLocation(shaderProgram1, 1, "vertexTexCoord");  // vertexTexCoord gets values from Attrib Array 1
    glBindAttribLocation(shaderProgram1, 2, "vertexNormal");  // vertexNormalCoord gets values from Attrib Array 2
    
    // Connect the fragmentColor to the frame buffer memory
    glBindFragDataLocation(shaderProgram1, 0, "fragmentColor");	// fragmentColor goes to the frame buffer memory
    
    // Program packaging
    glLinkProgram(shaderProgram1);
    checkLinking(shaderProgram1);
    ///
	
    // Create vertex shader from string
    unsigned int vertexShader2 = glCreateShader(GL_VERTEX_SHADER);
    if (!vertexShader2) { printf("Error in vertex shader 2 creation\n"); exit(1); }
    
    glShaderSource(vertexShader2, 1, &vertexSource2, NULL);
    glCompileShader(vertexShader2);
    checkShader(vertexShader2, "Vertex shader 2 error");
    
    // Create fragment shader from string
    unsigned int fragmentShader2 = glCreateShader(GL_FRAGMENT_SHADER);
    if (!fragmentShader2) { printf("Error in fragment shader 2 creation\n"); exit(1); }
    
    glShaderSource(fragmentShader2, 1, &fragmentSource2, NULL);
    glCompileShader(fragmentShader2);
    checkShader(fragmentShader2, "Fragment shader 2 error");
    
    // Attach shaders to a single program
    shaderProgram2 = glCreateProgram();
    if (!shaderProgram2) { printf("Error in shader program 2 creation\n"); exit(1); }
    
    glAttachShader(shaderProgram2, vertexShader2);
    glAttachShader(shaderProgram2, fragmentShader2);
    
    // Connect Attrib Arrays to input variables of the vertex shader
    glBindAttribLocation(shaderProgram2, 0, "vertexPosition");	// vertexPosition gets values from Attrib Array 0
    glBindAttribLocation(shaderProgram2, 1, "vertexTexCoord");  // vertexTexCoord gets values from Attrib Array 1
    glBindAttribLocation(shaderProgram2, 2, "vertexNormal");  // vertexNormalCoord gets values from Attrib Array 2
    
    // Connect the fragmentColor to the frame buffer memory
    glBindFragDataLocation(shaderProgram2, 0, "fragmentColor");	// fragmentColor goes to the frame buffer memory
    
    // Program packaging
    glLinkProgram(shaderProgram2);
    checkLinking(shaderProgram2);
    ///
    
    //envShaderProg
    
    // Create vertex shader from string
    unsigned int vertexShader3 = glCreateShader(GL_VERTEX_SHADER);
    if (!vertexShader3) { printf("Error in vertex shader 3 creation\n"); exit(1); }
    
    glShaderSource(vertexShader3, 1, &vertexSource3, NULL);
    glCompileShader(vertexShader3);
    checkShader(vertexShader3, "Vertex shader 3 error");
    
    // Create fragment shader from string
    unsigned int fragmentShader3 = glCreateShader(GL_FRAGMENT_SHADER);
    if (!fragmentShader3) { printf("Error in fragment shader 3 creation\n"); exit(1); }
    
    glShaderSource(fragmentShader3, 1, &fragmentSource3, NULL);
    glCompileShader(fragmentShader3);
    checkShader(fragmentShader3, "Fragment shader 3 error");
    
    // Attach shaders to a single program
    shaderProgram3 = glCreateProgram();
    if (!shaderProgram3) { printf("Error in shader program 3 creation\n"); exit(1); }
    
    glAttachShader(shaderProgram3, vertexShader3);
    glAttachShader(shaderProgram3, fragmentShader3);
    
    // Connect Attrib Arrays to input variables of the vertex shader
    glBindAttribLocation(shaderProgram3, 0, "vertexPosition");	// vertexPosition gets values from Attrib Array 0
    glBindAttribLocation(shaderProgram3, 1, "vertexTexCoord");  // vertexTexCoord gets values from Attrib Array 1
    glBindAttribLocation(shaderProgram3, 2, "vertexNormal");  // vertexNormalCoord gets values from Attrib Array 2
    
    // Connect the fragmentColor to the frame buffer memory
    glBindFragDataLocation(shaderProgram3, 0, "fragmentColor");	// fragmentColor goes to the frame buffer memory
    
    // Program packaging
    glLinkProgram(shaderProgram3);
    checkLinking(shaderProgram3);
    //////////////////
    
    scene.Initialize();

	for(int i = 0; i < 256; i++) keyPressed[i] = false;
}


void onKeyboard(unsigned char key, int x, int y)
{
	keyPressed[key] = true;
	glutPostRedisplay();
    
    scene.Control();
    
    
}

void onKeyboardUp(unsigned char key, int x, int y)
{
	keyPressed[key] = false;
	glutPostRedisplay();
    
    scene.Control();
}

void onExit() {
	glDeleteProgram(shaderProgram0);
	printf("exit");
}

void onDisplay() {
	
	glClearColor(0, 0.3, 0.5, 0);	 
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

	scene.Draw();
    scene.DrawShadow();

	glutSwapBuffers();									
}

void onReshape(int winWidth0, int winHeight0) 
{
	glViewport(0, 0, winWidth0, winHeight0);

	windowWidth = winWidth0, windowHeight = winHeight0;

	camera.SetAspectRatio((float)windowWidth / windowHeight);
}


void onIdle( ) {
    double t = glutGet(GLUT_ELAPSED_TIME) * 0.001;
    static double lastTime = 0.0;
    double dt = t - lastTime;
    lastTime = t;

	scene.Control();
	scene.Move(dt);
	scene.Draw();
    scene.DrawShadow();


    glutPostRedisplay();
}


int main(int argc, char * argv[]) {
	glutInit(&argc, argv);
#if !defined(__APPLE__)
	glutInitContextVersion(majorVersion, minorVersion);
#endif
	glutInitWindowSize(windowWidth, windowHeight);				
	glutInitWindowPosition(100, 100);							
#if defined(__APPLE__)
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_3_2_CORE_PROFILE);  
#else
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
#endif
	glutCreateWindow("3D Mesh");

#if !defined(__APPLE__)
	glewExperimental = true;
	glewInit();
#endif
	
	printf("GL Vendor    : %s\n", glGetString(GL_VENDOR));
	printf("GL Renderer  : %s\n", glGetString(GL_RENDERER));
	printf("GL Version (string)  : %s\n", glGetString(GL_VERSION));
	glGetIntegerv(GL_MAJOR_VERSION, &majorVersion);
	glGetIntegerv(GL_MINOR_VERSION, &minorVersion);
	printf("GL Version (integer) : %d.%d\n", majorVersion, minorVersion);
	printf("GLSL Version : %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));
	
	onInitialization();

	glutDisplayFunc(onDisplay);                
	glutReshapeFunc(onReshape);
	glutKeyboardFunc(onKeyboard);
	glutKeyboardUpFunc(onKeyboardUp);
	glutIdleFunc(onIdle);

	glutMainLoop();
	onExit();
	return 1;
}

