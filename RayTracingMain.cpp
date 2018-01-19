//
//  main.cpp
//  RayTracingHW
//  May 22nd
//  Graphics, AIT Spring 2017
//
//  Created by Samana Shrestha on 5/5/17.
//  Copyright © 2017 Samana Shrestha. All rights reserved.
//


#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vec2.h"
#include "vec3.h"
#include "vec4.h"
#include "mat4x4.h"

#include <algorithm>
#include "vec3.h"

using namespace std;

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

#include <vector>

//const unsigned int windowWidth = 384, windowHeight = 384;

const unsigned int windowWidth = 980, windowHeight = 1024;


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

const char *vertexSource0 = "\n\
#version 410 \n\
precision highp float; \n\
\n\
in vec2 vertexPosition;	\n\
in vec2 vertexTexCoord; \n\
out vec2 texCoord; \n\
\n\
void main() \n\
{ \n\
texCoord = vertexTexCoord; \n\
gl_Position = vec4(vertexPosition.x, vertexPosition.y, 0, 1); \n\
} \n\
";

const char *fragmentSource0 = " \n\
#version 410 \n\
precision highp float; \n\
\n\
uniform sampler2D samplerUnit; \n\
in vec2 texCoord;  \n\
out vec4 fragmentColor; \n\
\n\
void main() { \n\
fragmentColor = texture(samplerUnit, texCoord);  \n\
} \n\
";

unsigned int shaderProgram0;


// image to be computed by ray tracing
vec3 image[windowWidth * windowHeight];


// Ray structure.
class Ray
{
public:
    vec3 origin;
    vec3 dir;
    Ray(vec3 o, vec3 d)
    {
        origin = o;
        dir = d;
    }
};


class LightSource {
public:
    virtual vec3 getPowerDensityAt(vec3 x) = 0;
    virtual vec3 getLightDirAt(vec3 x) = 0;
    virtual float  getDistanceFrom(vec3 x) = 0;
};

class DirectionalLight : public LightSource {
    vec3 powerDensity, lightDir;
public:
    DirectionalLight(vec3 powerD, vec3 dir) : powerDensity(powerD), lightDir(dir) {}
    
    vec3 getPowerDensityAt(vec3 x) {
        return powerDensity;
    }
    
    vec3 getLightDirAt(vec3 x) {
        return lightDir.normalize();

    }
    
    float getDistanceFrom(vec3) {
        return 10000;
    }
};

class PointLight : public LightSource {
    vec3 powerDensity, lightPos;
public:
    PointLight(vec3 powerD, vec3 pos) : powerDensity(powerD), lightPos(pos) {}
    vec3 getPowerDensityAt(vec3 x) {
        float dist = (lightPos - x).norm();
        return powerDensity*(1 / dist);
    }
    vec3 getLightDirAt(vec3 x) {
        return (lightPos - x).normalize();
    }
    
    float getDistanceFrom(vec3 x) {
        return (lightPos - x).norm();
    }
};


// simple material class, with object color, and headlight shading
class Material
{
    vec3 innerColor, outerColor;
    
public:
    Material(vec3 outerColor, vec3 innerColor):
    innerColor(innerColor),
    outerColor(outerColor)
    {
    }

    
public:
    virtual vec3 getOuterColor()
    {
        return outerColor;
    }
    
    virtual vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray) {
        vec3 M = light->getPowerDensityAt(position);
        float dotVal = max(0.0f,normal.dot(light->getLightDirAt(position)));
        return M * outerColor * dotVal;
    }
    
    virtual bool reflective(){
        return false;
    }
    
    
    
};

class DiffuseMaterial : public Material
{
public:
    DiffuseMaterial():
    Material(vec3(1, 1, 1), vec3 (1, 1, 1))
    {
    }
    
    virtual vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray) {
        vec3 M = light->getPowerDensityAt(position);
        float dotVal = max(0.0f,normal.dot(light->getLightDirAt(position)));
        vec3 outerColor = getOuterColor();
        return M * outerColor * dotVal;
    }
    
};



float snoise(vec3 r) {
    unsigned int x = 0x0625DF73;
    unsigned int y = 0xD1B84B45;
    unsigned int z = 0x152AD8D0;
    float f = 0;
    for(int i=0; i<32; i++) {
        vec3 s(	x/(float)0xffffffff,
               y/(float)0xffffffff,
               z/(float)0xffffffff);
        f += sin(s.dot(r));
        x = x << 1 | x >> 31;
        y = y << 1 | y >> 31;
        z = z << 1 | z >> 31;
    }
    return f / 64.0 + 0.5;
}

class Ocean : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Ocean():
    Material(vec3(1, 1, 1), vec3 (1, 1, 1))
    {
        scale = 2;
        turbulence = 5;
        period = 2;
        sharpness = 1;
    }
    virtual vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray)
    {
        float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence;
        w = pow(sin(w)*0.5+0.5, 4);
        vec3 oColor = (vec3(4, 4, 4) * w + vec3(1.5, 1.90, 1.99) * (1-w)) * normal.dot(light->getLightDirAt(position));

//        vec3 oColor = (vec3(1, 1, 1) * w + vec3(0.5, 0.90, 0.99) * (1-w)) * normal.dot(light->getLightDirAt(position));

        vec3 h = (light->getLightDirAt(position) + (-ray.dir)).normalize();
        vec3 M = light->getPowerDensityAt(position);
        float spec = pow(max(normal.dot(h), 0.0f), 16.0);
        float dotVal = max(0.0f,normal.dot(light->getLightDirAt(position)));
        vec3 diffuse = M * oColor * dotVal;
        return (M * oColor * spec) + diffuse;
    }
};


class PhongMaterial : public Material
{
public:
    PhongMaterial():
    Material(vec3(1, 1, 1), vec3 (1, 1, 1))
    {
    }
    
    vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray) {
        
        float theta = (acos(normal.z/normal.norm()));
        theta = (theta / M_PI) * 100;
        float phi = (atan(normal.y/normal.x));
        phi = (phi / M_PI) * 180;
        vec3 finalColor;
        
        if (phi >-90.0 && phi < -45.0 && theta > 0.0 && theta < 180.0){
            finalColor = vec3(1,0,0.2);
        }
        else if (phi >-45.0 && phi < 0.0 && theta > 0.0 && theta < 180.0){
            finalColor = vec3(0,0.5,1);
        }
        else if (phi >0.0 && phi < 45.0 && theta > 0.0 && theta < 180.0){
            finalColor = vec3(0,1,0.4);
        }
        else{
            finalColor = vec3(1,1,0);
        }
        vec3 oColor = finalColor * fmax(0, normal.dot(light->getLightDirAt(position)));
        
        
        vec3 h = (light->getLightDirAt(position) + (-ray.dir)).normalize();
        vec3 M = light->getPowerDensityAt(position);
        float spec = pow(max(normal.dot(h), 0.0f), 32.0);
        float dotVal = max(0.0f,normal.dot(light->getLightDirAt(position)));
        vec3 diffuse = M * oColor * dotVal;
        
        return (M * oColor * spec) + diffuse;
    }
    
};


class Flare : public Material
{
public:
    Flare():
    Material(vec3(1, 1, 1), vec3 (1, 1, 1))
    {
    }
    
    vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray) {
        
        vec3 oColor = vec3(1, 1, 1) * fmax(0, normal.dot(light->getLightDirAt(position)));
        vec3 h = (light->getLightDirAt(position) + (-ray.dir)).normalize();
        vec3 M = light->getPowerDensityAt(position);
        float spec = pow(max(normal.dot(h), 0.0f), 2.0);
        float dotVal = max(0.0f,normal.dot(light->getLightDirAt(position)));
        vec3 diffuse = M * oColor * dotVal;
        
        return (M * oColor * spec);
    }
    
};



class Wood : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Wood():
    Material(vec3(1, 1, 1), vec3 (1, 1, 1))
    {
        scale = 16;
        turbulence = 500;
        period = 8;
        sharpness = 10;
    }

    
    virtual vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray)
    {
        float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence + 10000.0;
        w -= int(w);
        return (vec3(0.42, 0.31, 0.22) * w + vec3(0, 0, 0.05) * (1-w)) * normal.dot(light->getLightDirAt(position));
    }
};


class Marble : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Marble():
    Material(vec3(1, 1, 1), vec3(1, 1, 1))
    {
        scale = 32;
        turbulence = 50;
        period = 32;
        sharpness = 1;
    }
    virtual vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray)
    {
        float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence;
        w = pow(sin(w)*0.5+0.5, 4);
        return (vec3(0, 0, 1) * w + vec3(1, 1, 1) * (1-w)) * normal.dot(light->getLightDirAt(position));
    }
};

class Stone : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Stone():
    Material(vec3(1, 1, 1), vec3(1, 1, 1))
    {
        scale = 32;
        turbulence = 50;
        period = 32;
        sharpness = 1;
    }
    virtual vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray)
    {
        float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence;
        w = pow(sin(w)*0.5+0.5, 4);
        vec3 oColor = (vec3(0.1, 0.2, 0.3) * w + vec3(0.47, 0.46, 0.45) * (1-w)) * normal.dot(light->getLightDirAt(position));
        
        vec3 h = (light->getLightDirAt(position) + (-ray.dir)).normalize();
        vec3 M = light->getPowerDensityAt(position);
        float spec = pow(max(normal.dot(h), 0.0f), 8.0);
        float dotVal = max(0.0f,normal.dot(light->getLightDirAt(position)));
        vec3 diffuse = M * oColor * dotVal;
        
        return (M * oColor * spec) + diffuse;

    }
};


class Sand : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Sand():
    Material(vec3(1.5, 0.8, 0), vec3(1, 1, 1))
    {
        scale = 300;
        turbulence = 50;
        period = 1000;
        sharpness = 1;
    }
    virtual vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray)
    {
        float w = position.x * period + pow(snoise(position * scale), sharpness)*turbulence;

        w = pow(sin(w)*0.5+0.5, 4);
        return (vec3(0.38, 0.38, 0.34) * w + vec3(0.5, 0.38, 0) * (1-w)) * normal.dot(light->getLightDirAt(position));
    }
};


class ParaStripes : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    ParaStripes():
    Material(vec3(1, 1, 1), vec3(1, 1, 1))
    {
    }
    virtual vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray)
    {
        float theta = (acos(normal.z/normal.norm()));
        theta = (theta / M_PI) * 100;
        
        float phi = (atan(normal.y/normal.x));
        phi = (phi / M_PI) * 180;
        
        
        vec3 finalColor;
        
        if (phi >-90.0 && phi < -45.0 && theta > 0.0 && theta < 180.0){
            finalColor = vec3(0.8, 0.5, 1.05);
        }
        else if (phi >-45.0 && phi < 0.0 && theta > 0.0 && theta < 180.0){
            finalColor = vec3(1.42, 1.31, 1.22);
        }
        else if (phi >0.0 && phi < 45.0 && theta > 0.0 && theta < 180.0){
            finalColor = vec3(0.8, 0.5, 1.05);
        }
        else{
            finalColor = vec3(1.42, 1.31, 1.22);
        }
        
        return finalColor * fmax(0, normal.dot(light->getLightDirAt(position)));
    }
};


class Stripes : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Stripes():
    Material(vec3(0.8, 0.7, 0), vec3 (1, 1, 1))
    {
        scale = 16;
        turbulence = 0;
        period = 8;
        sharpness = 90;
    }
    virtual vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray)
    {
        float w = position.x * period + pow(snoise(position * scale), sharpness) + 10000.0;
        w -= int(w);
        return (vec3(2, 2, 0) * w + vec3(1, 0, 0.3) * (1-w)) * normal.dot(light->getLightDirAt(position));
    }
};

class Wavy1 : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Wavy1():
    Material(vec3(0.8, 0.7, 0), vec3 (1, 1, 1))
    {
        scale = 32;
        turbulence = 50;
        period = 2;
        sharpness = 1;
    }
    virtual vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray)
    {
        float w = position.x * period + pow(snoise(position * scale), sharpness) + 10000.0;
        w -= int(w);
        vec3 oColor = (vec3(0.6, 7, 0) * w + vec3(0, 0.9, 0.4) * (1-w)) * normal.dot(light->getLightDirAt(position));
        
        
        vec3 h = (light->getLightDirAt(position) + (-ray.dir)).normalize();
        vec3 M = light->getPowerDensityAt(position);
        float spec = pow(max(normal.dot(h), 0.0f), 16.0);
        float dotVal = max(0.0f,normal.dot(light->getLightDirAt(position)));
        vec3 diffuse = M * oColor * dotVal;
        
        return (M * oColor * spec) + diffuse;
        
    }
};


class Flotsam : public Material
{
    float scale;
    float turbulence;
    float period;
    float sharpness;
public:
    Flotsam():
    Material(vec3(0.8, 0.7, 0), vec3 (1, 1, 1))
    {
        scale = 32;
        turbulence = 50;
        period = 200;
        sharpness = 1;
    }
    virtual vec3 shade(vec3 position, vec3 normal, LightSource* light, Ray ray)
    {
        float w = position.x * period + pow(snoise(position * scale), sharpness) + 10000.0;
        w -= int(w);
        vec3 oColor = (vec3(0.40, 0.30, 0.16) * w + vec3(0.8, 0.8, 0.8) * (1-w)) * normal.dot(light->getLightDirAt(position));
        
        
        vec3 h = (light->getLightDirAt(position) + (-ray.dir)).normalize();
        vec3 M = light->getPowerDensityAt(position);
        float spec = pow(max(normal.dot(h), 0.0f), 16.0);
        float dotVal = max(0.0f,normal.dot(light->getLightDirAt(position)));
        vec3 diffuse = M * oColor * dotVal;
        
        return diffuse;
        
    }
};



// Camera class.
class Camera
{
    vec3 eye;		//< world space camera position
    vec3 lookAt;	//< center of window in world space
    vec3 right;		//< vector from window center to window right-mid (in world space)
    vec3 up;		//< vector from window center to window top-mid (in world space)
    
public:
    Camera()
    {
        eye = vec3(0.0, 0.0, 10.0);
        lookAt = vec3(0.0, 0.0, 2.0);
        right = vec3(1, 0, 0);
        up = vec3(0.0, 1.0, 0.0);
    }
    vec3 getEye()
    {
        return eye;
    }
    // compute ray through pixel at normalized device coordinates
    vec3 rayDirFromNdc(float x, float y) {
        return (lookAt - eye
                + right * x
                + up    * y
                ).normalize();
    }
};


// Hit record structure. Contains all data that describes a ray-object intersection point.
class Hit
{
public:
    Hit()
    {
        t = -1;
    }
    float t;				//< Ray paramter at intersection. Negative means no valid intersection.
    vec3 position;			//< Intersection coordinates.
    vec3 normal;			//< Surface normal at intersection.
    Material* material;		//< Material of intersected surface.
};

// Abstract base class.
class Intersectable
{
protected:
    Material* material;
public:
    Intersectable(Material* material):material(material) {}
    virtual Hit intersect(const Ray& ray)=0;
};

// Simple helper class to solve quadratic equations with the Quadratic Formula [-b +- sqrt(b^2-4ac)] / 2a, and store the results.
class QuadraticRoots
{
public:
    float t1;
    float t2;
    // Solves the quadratic a*t*t + b*t + c = 0 using the Quadratic Formula [-b +- sqrt(b^2-4ac)] / 2a, and sets members t1 and t2 to store the roots.
    QuadraticRoots(float a, float b, float c)
    {
        float discr = b * b - 4.0 * a * c;
        if ( discr < 0 ) // no roots
        {
            t1 = -1;
            t2 = -1;
            return;
        }
        float sqrt_discr = sqrt( discr );
        t1 = (-b + sqrt_discr)/2.0/a;
        t2 = (-b - sqrt_discr)/2.0/a;
    }
    // Returns the lesser of the positive solutions, or a negative value if there was no positive solution.
    float getLesserPositive()
    {
        return (0 < t1 && (t2 < 0 || t1 < t2)) ? t1 : t2;
    }
};

// Object realization.
class Sphere : public Intersectable
{
    vec3 center;
    float radius;
public:
    Sphere(const vec3& center, float radius, Material* material):
    Intersectable(material),
    center(center),
    radius(radius)
    {
    }
    QuadraticRoots solveQuadratic(const Ray& ray)
    {
        vec3 diff = ray.origin - center;
        float a = ray.dir.dot(ray.dir);
        float b = diff.dot(ray.dir) * 2.0;
        float c = diff.dot(diff) - radius * radius;
        return QuadraticRoots(a, b, c);
        
    }
    vec3 getNormalAt(vec3 r)
    {
        return (r - center).normalize();
    }
    Hit intersect(const Ray& ray)
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal
        float t = solveQuadratic(ray).getLesserPositive();
        
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = getNormalAt(hit.position);
        
        return hit;
    }
    
    
    
    
};

// CLASS PLANE

class Plane : public Intersectable
{
    vec3 normal, planePoint;
public:
    Plane(const vec3& normal, const vec3& planePoint, Material* material):
    Intersectable(material),
    normal(normal),
    planePoint(planePoint)
    {
    }

    vec3 getNormalAt(vec3 r)
    {
        return (r - planePoint).normalize();
    }
    Hit intersect(const Ray& ray)
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal
        float t;
        
        Hit hit;
        
        vec3 e = ray.origin;
        vec3 d = ray.dir;
        
        t = ((planePoint - e).dot(normal))/(d.dot(normal));
        
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = normal;
        
        return hit;
    }
    
    
    
    
};


// CLASS QUADRIC

class Quadric : public Intersectable
{
    mat4x4 coeffs;
    
public:
    Quadric(Material* material):
    
    Intersectable(material)
    {
        coeffs = mat4x4(1, 0, 0, 0,
                        0, 1, 0, 0,
                        0, 0, 1, 0,
                        0, 0, 0, -1);
    }
    
    QuadraticRoots solveQuadratic(const Ray& ray)
    {
        vec4 d = vec4(ray.dir.x, ray.dir.y, ray.dir.z, 0);
        vec4 e = vec4(ray.origin.x, ray.origin.y, ray.origin.z, 1);
        float a = d.dot(coeffs * d);
        float b = (d.dot(coeffs * e)) + (e.dot(coeffs * d));
        float c = e.dot(coeffs * e);
        return QuadraticRoots(a, b, c);
        
    }
    
    vec3 getNormalAt(vec3 position)
    {
        vec4 pos = vec4(position.x, position.y, position.z, 1);
        vec4 norm = (coeffs * pos) + (pos * coeffs);
        
        return vec3(norm.x, norm.y, norm.z).normalize();
    }
    
    
    
    Hit intersect(const Ray& ray)
    
    {
        // This is a generic intersect that works for any shape with a quadratic equation. solveQuadratic should solve the proper equation (+ ray equation) for the shape, and getNormalAt should return the proper normal
        
        float t = solveQuadratic(ray).getLesserPositive();
        Hit hit;
        
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = getNormalAt(hit.position);
        
        return hit;
    }
    
    
    
    Quadric* transform(mat4x4 t){
        mat4x4 tInv = t.invert();
        coeffs = tInv * coeffs * tInv.transpose();

        return this;
    }
    
    
    Quadric* sphere(){
        coeffs = mat4x4(1, 0, 0, 0,
                        0, 1, 0, 0,
                        0, 0, 1, 0,
                        0, 0, 0, -1);
        return this;
    }

    
    Quadric* sphere(float r){
        coeffs._33 = - (r*r);
        return this;
    }
    
    
    Quadric* ellipsoid(){
        coeffs = mat4x4(0.2, 0, 0, 0,
                        0, 0.8, 0, 0,
                        0, 0, 0.3, 0,
                        0, 0, 0, -1);
        return this;
    }
    
    
    Quadric* cylinder(){
        coeffs._11 = 0;
        return this;
    }
    
    Quadric* cylinder(float r){
        coeffs._11 = 0;
        coeffs._33 = - (r*r);
        return this;
    }
    
    
    
    Quadric* hyperB(){
        coeffs = mat4x4(1, 0, 0, 0,
                        0, -1, 0, 0,
                        0, 0, 1, 0,
                        0, 0, 0, -1);
        return this;
    }
    
    
    
    Quadric* parab(){
        coeffs = mat4x4(1, 0, 0, 0,
                        0, 0, 0, -1,
                        0, 0, 1, 0,
                        0, 0, 0, 0);
        return this;
    }
    
    
    
    Quadric* planeY(){
        coeffs._00 = 0;
        coeffs._11 = 1;
        coeffs._22 = 0;
        coeffs._33 = -0.23;
        return this;
    }
    
    
    Quadric* planeZ(){
        coeffs._00 = 0;
        coeffs._11 = 0;
        coeffs._22 = 1;
        coeffs._33 = -0.23;
        return this;
    }
    
    
    
    Quadric* planeX(){
        coeffs._00 = 1;
        coeffs._11 = 0;
        coeffs._22 = 0;
        coeffs._33 = -0.23;
        return this;
    }
    
    
    Quadric* hyperB(float r){
        coeffs._11 = -1;
        coeffs._33 = - (r*r);
        return this;
    }
    
    
    bool contains(vec3 r){
        vec4 rHomo(r);
        float res = rHomo.dot(coeffs * rHomo);
        
        if(res < 0){
            return true;
        }else {
            return false;
        }
    }
    
    
    
    Quadric* parallelPlanes(){
        coeffs._00 = 0;
        coeffs._11 = 1;
        coeffs._22 = 0;
        coeffs._33 = -1;
        return this;
    }

    
    Quadric* parallelPlanes(float height){
        coeffs._00 = 0;
        coeffs._11 = 1;
        coeffs._22 = 0;
        coeffs._33 = -(height/2)*(height/2);
        return this;
    }
    
    
    Quadric* parallelPlanesX(float height){
        coeffs._00 = 1;
        coeffs._11 = 0;
        coeffs._22 = 0;
        coeffs._33 = -(height/2)*(height/2);
        return this;
    }

    
    Quadric* parallelPlanesZ(float height){
        coeffs._00 = 0;
        coeffs._11 = 0;
        coeffs._22 = 1;
        coeffs._33 = -(height/2)*(height/2);
        return this;
    }

    
    Quadric* parallelPlanes(vec3 moveTo){
        coeffs._00 = 0;
        coeffs._11 = 1;
        coeffs._22 = 0;
        coeffs._33 = -1;
        return this->transform(mat4x4::scaling(vec3(1, 1, 1))                                             * mat4x4::translation(moveTo));
    }

    
};

// CLASS CLIPPEDQUADRIC


class ClippedQuadric : public Intersectable
{
    Quadric shape;
    Quadric clipper;
    Quadric clipper2;

public:
    
    ClippedQuadric(Material* material):
    Intersectable(material), shape(material), clipper(material), clipper2(material)
    {
        
    }
    
    
    Hit intersect(const Ray& ray)
    {
        //generic intersect that works for any shape
        
        QuadraticRoots roots = shape.solveQuadratic(ray);
        
        vec3 p1 = ray.origin + ray.dir * roots.t1;
        vec3 p2 = ray.origin + ray.dir * roots.t2;
        
        if(!clipper.contains(p1)) roots.t1 = -1;
        if(!clipper.contains(p2)) roots.t2 = -1;
    
        if(!clipper2.contains(p1)) roots.t1 = -1;
        if(!clipper2.contains(p2)) roots.t2 = -1;
        
        float t = roots.getLesserPositive();
        Hit hit;
        hit.t = t;
        hit.material = material;
        hit.position = ray.origin + ray.dir * t;
        hit.normal = shape.getNormalAt(hit.position);
        return hit;
    }
    
    
    
    ClippedQuadric* transform(mat4x4 t){
        shape.transform(t);
        clipper.transform(t);
        clipper2.transform(t);
        return this;
    }
    
    
    
    ClippedQuadric* cylinder(float r, float h){
        shape.cylinder(r);
        clipper.parallelPlanes(h);
        clipper2.parallelPlanes(h);
        return this;
    }
    
    
    
    ClippedQuadric* slabZ(float h){
        shape.planeZ();
        clipper.parallelPlanes(h);
        clipper2.parallelPlanesX(h);
        return this;
    }
    
    
    
    ClippedQuadric* slabX(float h){
        shape.planeX();
        clipper.parallelPlanes(h);
        clipper2.parallelPlanesZ(h);
        return this;
    }
    
    ClippedQuadric* slabY(float h){
        shape.planeY();
        clipper.parallelPlanesZ(h);
        clipper2.parallelPlanesX(h);
        return this;
    }

    
    ClippedQuadric* hemisphere(float r, float h){
        shape.sphere(r);
        clipper.parallelPlanes(vec3(0, 1, 0));
        clipper2.parallelPlanes(vec3(0, 1, 0));
        return this;
    }
    
    ClippedQuadric* hemisphere2(float r, float h){
        shape.sphere(r);
        clipper.parallelPlanes(vec3(0, 1.4, 0));
        clipper2.parallelPlanes(vec3(0, 1.4, 0));
        return this;
    }
    
    
    
    ClippedQuadric* clipsphere(float r, float h){
        shape.sphere(r);
        clipper.parallelPlanes(h);
        clipper2.parallelPlanes(h);
        return this;
    }
    
    
    ClippedQuadric* hyperB(float r, float h){
        shape.hyperB(r);
        clipper.parallelPlanes(vec3(0, -1, 0));
        clipper2.parallelPlanes(vec3(0, -1, 0));
        return this;
    }
    
};




class Scene
{
    Camera camera;
    std::vector<Intersectable*> objects;
    	std::vector<Material*> materials;
    std::vector<LightSource*> lightSources;
public:
    Scene()
    {
        //Materials
        materials.push_back(new Material(vec3(0.2, 0.6, 1), vec3(4, 0, 3)));
        materials.push_back(new Material(vec3(1, 1, 1), vec3(0, 1, 9)));
        materials.push_back(new Flotsam());
        materials.push_back(new Flare());
        materials.push_back(new Wood());
        materials.push_back(new Marble());
        materials.push_back(new Stripes());
        materials.push_back(new Wavy1());
        materials.push_back(new ParaStripes());
        materials.push_back(new Sand());
        materials.push_back(new PhongMaterial());
        materials.push_back(new Ocean());
        materials.push_back(new Stone());
        
        //LightSources
        lightSources.push_back(new DirectionalLight(vec3(0.9, 0.9, 0.9), vec3(1, 1, 0)));
        lightSources.push_back(new DirectionalLight(vec3(0.9, 0.9, 0.9), vec3(0, 0, 1)));
//        lightSources.push_back(new DirectionalLight(vec3(0.9, 0.9, 0.9), vec3(1, 0, 1)));
//        lightSources.push_back(new DirectionalLight(vec3(0.9, 0.9, 0.9), vec3(0, 1, 0)));
        lightSources.push_back(new PointLight(vec3(1, 0, 0), vec3(1.04, 0.8, 0.7)));
        
        
        
        
        //Stones
        Quadric* stone = (new Quadric(materials[12]))->ellipsoid();
        objects.push_back(stone->transform(mat4x4::scaling(vec3(0.03, 0.02, 0.03)) *
                                          mat4x4::translation(vec3(0.89, -0.7, 1.88))));
        Quadric* stone1 = (new Quadric(materials[12]))->ellipsoid();
        objects.push_back(stone1->transform(mat4x4::scaling(vec3(0.03, 0.02, 0.035)) *
                                           mat4x4::translation(vec3(0.7, -0.72, 1.9))));
        Quadric* stone2 = (new Quadric(materials[12]))->ellipsoid();
        objects.push_back(stone2->transform(mat4x4::scaling(vec3(0.035, 0.02, 0.04)) *
                                           mat4x4::translation(vec3(0.8, -0.72, 1.9))));
        Quadric* stone3 = (new Quadric(materials[12]))->ellipsoid();
        objects.push_back(stone3->transform(mat4x4::scaling(vec3(0.045, 0.02, 0.035)) *
                                            mat4x4::rotation(vec3(0, 0, 1), 0.2) *
                                           mat4x4::translation(vec3(0.79, -0.75, 1.88))));
        Quadric* stone4 = (new Quadric(materials[12]))->ellipsoid();
        objects.push_back(stone4->transform(mat4x4::scaling(vec3(0.045, 0.02, 0.035)) *
                                            mat4x4::rotation(vec3(0, 0, 1), 0.2) *
                                            mat4x4::translation(vec3(0.92, -0.78, 1.8))));
        Quadric* stone5 = (new Quadric(materials[12]))->ellipsoid();
        objects.push_back(stone5->transform(mat4x4::scaling(vec3(0.045, 0.02, 0.035)) *
                                            mat4x4::rotation(vec3(0, 0, 1), 0.2) *
                                            mat4x4::translation(vec3(0.94, -0.75, 1.9))));
        
        
        //FlotsamBox
        ClippedQuadric* boxX = new ClippedQuadric(materials[2]);
        objects.push_back(boxX->slabX(1)->transform(mat4x4::scaling(vec3(0.2, 0.2, 0.2)) *
                                                    mat4x4::rotation(vec3(0, 1, 0), 0.8) *
                                                    mat4x4::rotation(vec3(1, 0, 0), 0.5) *
                                                    mat4x4::translation(vec3(-0.6, -0.6, 2))));
        ClippedQuadric* boxY = new ClippedQuadric(materials[2]);
        objects.push_back(boxY->slabY(1)->transform(mat4x4::scaling(vec3(0.2, 0.2, 0.2)) *
                                                    mat4x4::rotation(vec3(0, 1, 0), 0.8) *
                                                    mat4x4::rotation(vec3(1, 0, 0), 0.5) *
                                                    mat4x4::translation(vec3(-0.6, -0.6, 2))));
        
        
        

        
        //Flare
        ClippedQuadric* flare = (new ClippedQuadric(materials[3]))->clipsphere(2, 1.8);
        objects.push_back(flare->transform(mat4x4::scaling(vec3(0.02, 0.05, 0.05)) *
                                          mat4x4::translation(vec3(1.08, 0.8, 0.6))));
        
        ClippedQuadric* flareCylin = (new ClippedQuadric(materials[1]))->cylinder(1, 2);
        objects.push_back(flareCylin->transform(mat4x4::scaling(vec3(0.015, 0.9, 0.015)) *
                                           mat4x4::translation(vec3(1.2, 0, -0.5))));
        
        
        //Ocean
        objects.push_back(new Plane(vec3(0, 1, 0), vec3(0.5, -5, 0), materials[11]));


        
        //Parasol
        ClippedQuadric* paraCylin = new ClippedQuadric(materials[1]);
        objects.push_back(paraCylin->cylinder(1, 2)->transform(mat4x4::scaling(vec3(0.015, 0.4, 0.015))*
                                                               mat4x4::rotation(vec3(0, 0, 1), 0.18) *
                                                               mat4x4::translation(vec3(-0.4, -0.2, 1.22))));
        ClippedQuadric* paraTop = (new ClippedQuadric(materials[8]))->hemisphere2(1, 1);
        objects.push_back(paraTop->transform(mat4x4::scaling(vec3(0.54, 0.4, 0.2)) *
                                          mat4x4::rotation(vec3(1, 0, 0), -0.15) *
                                          mat4x4::rotation(vec3(0, 0, 1), 0.3) *
                                          mat4x4::translation(vec3(-0.4, 0, 1.4))));
        
        
        
        //SandCastle
        
        ClippedQuadric* casPillar1 = (new ClippedQuadric(materials[9]))->cylinder(1, 2);
        objects.push_back(casPillar1->transform(mat4x4::scaling(vec3(0.06, 0.16, 0.06))*
                                                mat4x4::translation(vec3(0.3, -0.4, 1.9))));
        ClippedQuadric* casPillar2 = (new ClippedQuadric(materials[9]))->cylinder(1, 2);
        objects.push_back(casPillar2->transform(mat4x4::scaling(vec3(0.06, 0.16, 0.06))*
                                                mat4x4::translation(vec3(-0.2, -0.4, 1.9))));
        ClippedQuadric* casPillar3 = (new ClippedQuadric(materials[9]))->cylinder(1, 2);
        objects.push_back(casPillar3->transform(mat4x4::scaling(vec3(0.08, 0.32, 0.08))*
                                                mat4x4::translation(vec3(0.2, -0.4, 1.4))));
        ClippedQuadric* casPillar4 = (new ClippedQuadric(materials[9]))->cylinder(1, 2);
        objects.push_back(casPillar4->transform(mat4x4::scaling(vec3(0.08, 0.32, 0.08))*
                                                mat4x4::translation(vec3(-0.1, -0.4, 1.4))));
        
        
        
        ClippedQuadric* bastion1 = (new ClippedQuadric(materials[9]))->hyperB(0, 1);
        objects.push_back(bastion1->transform(mat4x4::scaling(vec3(0.04, 0.08, 0.04))*
                                              mat4x4::translation(vec3(0.3, -0.12, 1.9))));
        ClippedQuadric* bastion2 = (new ClippedQuadric(materials[9]))->hyperB(0, 1);
        objects.push_back(bastion2->transform(mat4x4::scaling(vec3(0.04, 0.08, 0.04))*
                                              mat4x4::translation(vec3(-0.2, -0.12, 1.9))));
        ClippedQuadric* bastion3 = (new ClippedQuadric(materials[9]))->hyperB(0, 1);
        objects.push_back(bastion3->transform(mat4x4::scaling(vec3(0.06, 0.08, 0.06))*
                                              mat4x4::translation(vec3(0.2, 0.08, 1.4))));
        ClippedQuadric* bastion4 = (new ClippedQuadric(materials[9]))->hyperB(0, 1);
        objects.push_back(bastion4->transform(mat4x4::scaling(vec3(0.06, 0.08, 0.08))*
                                              mat4x4::translation(vec3(-0.1, 0.08, 1.4))));
        
        
        
        ClippedQuadric* body = (new ClippedQuadric(materials[9]))->clipsphere(2, 2);
        objects.push_back(body->transform(mat4x4::scaling(vec3(0.1, 0.1, 0.1)) *
                                          mat4x4::translation(vec3(0.05, -0.45, 1.9))));
        ClippedQuadric* body1 = (new ClippedQuadric(materials[9]))->clipsphere(2, 2);
        objects.push_back(body1->transform(mat4x4::scaling(vec3(0.08, 0.08, 0.08)) *
                                           mat4x4::translation(vec3(0.05, -0.3, 1.9))));
        ClippedQuadric* body2 = (new ClippedQuadric(materials[9]))->clipsphere(2, 2);
        objects.push_back(body2->transform(mat4x4::scaling(vec3(0.06, 0.06, 0.06)) *
                                           mat4x4::translation(vec3(0.05, -0.18, 1.9))));
        
        ClippedQuadric* flagPole = (new ClippedQuadric(materials[9]))->cylinder(1, 1);
        objects.push_back(flagPole->transform(mat4x4::scaling(vec3(0.01, 0.2, 0.01))*
                                                mat4x4::translation(vec3(0.05, -0.1, 1.9))));
        ClippedQuadric* flag = (new ClippedQuadric(materials[5]))->hyperB(0, 1);
        objects.push_back(flag->transform(mat4x4::scaling(vec3(0.04, 0.04, 0.01))*
                                              mat4x4::rotation(vec3(0, 0, 1), 0.8) *
                                              mat4x4::translation(vec3(0.045, 0.08, 1.9))));
        

        
        
        //SandDune
        ClippedQuadric* sand1 = new ClippedQuadric(materials[9]);
        objects.push_back(sand1->hemisphere(1, 1)->transform(mat4x4::scaling(vec3(1.4, 0.4, 0.5))                                             * mat4x4::translation(vec3(0.4, -0.8, 1.1))));
        
        

        //BeachBall
        Quadric* ball = (new Quadric(materials[10]))->sphere();
        objects.push_back(ball->transform(mat4x4::scaling(vec3(0.1, 0.1, 0.1)) *
                                          mat4x4::rotation(vec3(1, 0, 0), 0.5) *
                                          mat4x4::translation(vec3(0.6, -0.4, 1.8))));
        
        
        
        
        
        
        
        
        //TreeTrunk
        ClippedQuadric* treeTrunk0 = (new ClippedQuadric(materials[4]))->hyperB(0, 1);
        objects.push_back(treeTrunk0->transform(mat4x4::scaling(vec3(0.06, 0.2, 0.06)) * mat4x4::translation(vec3(0.9, -0.1, 1.2))));
        ClippedQuadric* treeTrunk1 = (new ClippedQuadric(materials[4]))->hyperB(0, 1);
        objects.push_back(treeTrunk1->transform(mat4x4::scaling(vec3(0.045, 0.16, 0.045))*
                                                mat4x4::rotation(vec3(0, 0, 1), 0.15) *
                                                mat4x4::translation(vec3(0.85, 0, 1.22))
                                                ));
        ClippedQuadric* treeTrunk2 = (new ClippedQuadric(materials[4]))->hyperB(0, 1);
        objects.push_back(treeTrunk2->transform(mat4x4::scaling(vec3(0.03, 0.15, 0.03)) *
                                                mat4x4::rotation(vec3(0, 0, 1), 0.2) *
                                                mat4x4::translation(vec3(0.82, 0.07, 1.2))));
        ClippedQuadric* treeTrunk3 = (new ClippedQuadric(materials[4]))->hyperB(0, 1);
        objects.push_back(treeTrunk3->transform(mat4x4::scaling(vec3(0.03, 0.15, 0.03)) *
                                                mat4x4::rotation(vec3(0, 0, 1), 0.25) *
                                                mat4x4::translation(vec3(0.8, 0.14, 1.2))));
        ClippedQuadric* treeTrunk4 = (new ClippedQuadric(materials[4]))->hyperB(0, 1);
        objects.push_back(treeTrunk4->transform(mat4x4::scaling(vec3(0.03, 0.15, 0.03)) *
                                                mat4x4::rotation(vec3(0, 0, 1), 0.3) *
                                                mat4x4::translation(vec3(0.76, 0.21, 1.2))));
        ClippedQuadric* treeTrunk5 = (new ClippedQuadric(materials[4]))->hyperB(0, 1);
        objects.push_back(treeTrunk5->transform(mat4x4::scaling(vec3(0.03, 0.15, 0.03))*
                                                mat4x4::rotation(vec3(0, 0, 1), 0.32) *
                                                mat4x4::translation(vec3(0.72, 0.3, 1.2))));
        ClippedQuadric* treeTrunk6 = (new ClippedQuadric(materials[4]))->hyperB(0, 1);
        objects.push_back(treeTrunk6->transform(mat4x4::scaling(vec3(0.03, 0.15, 0.03))*
                                                mat4x4::rotation(vec3(0, 0, 1), 0.36) *
                                                mat4x4::translation(vec3(0.68, 0.37, 1.2))));
        ClippedQuadric* treeTrunk7 = (new ClippedQuadric(materials[4]))->hyperB(0, 1);
        objects.push_back(treeTrunk7->transform(mat4x4::scaling(vec3(0.03, 0.15, 0.03))*
                                                mat4x4::rotation(vec3(0, 0, 1), 0.4) *
                                                mat4x4::translation(vec3(0.64, 0.44, 1.2))));
        ClippedQuadric* treeTrunk8 = (new ClippedQuadric(materials[4]))->hyperB(0, 1);
        objects.push_back(treeTrunk8->transform(mat4x4::scaling(vec3(0.03, 0.15, 0.03))*
                                                mat4x4::rotation(vec3(0, 0, 1), 0.44) *
                                                mat4x4::translation(vec3(0.6, 0.51, 1.2))));
        ClippedQuadric* treeTrunk9 = (new ClippedQuadric(materials[4]))->hyperB(0, 1);
        objects.push_back(treeTrunk9->transform(mat4x4::scaling(vec3(0.02, 0.14, 0.02))*
                                                mat4x4::rotation(vec3(0, 0, 1), 0.46) *
                                                mat4x4::translation(vec3(0.56, 0.58, 1.2))));
        ClippedQuadric* treeTrunk10 = (new ClippedQuadric(materials[4]))->hyperB(0, 1);
        objects.push_back(treeTrunk10->transform(mat4x4::scaling(vec3(0.02, 0.14, 0.02))*
                                                mat4x4::rotation(vec3(0, 0, 1), 0.48) *
                                                mat4x4::translation(vec3(0.52, 0.64, 1.2))));
        
        
        ///TreeLeaves
        ClippedQuadric* leaf1 = new ClippedQuadric(materials[7]);
        objects.push_back(leaf1->hemisphere(1, 1)->transform(mat4x4::scaling(vec3(0.24, 0.1, 0.01)) *
                                                             mat4x4::rotation(vec3(0, 0, 1), 0.6) *
                                                             mat4x4::translation(vec3(0.34, 0.46, 1.4))));
        ClippedQuadric* leaf2 = new ClippedQuadric(materials[7]);
        objects.push_back(leaf2->hemisphere(1, 1)->transform(mat4x4::scaling(vec3(0.24, 0.1, 0.01)) *
                                                             mat4x4::rotation(vec3(0, 0, 1), -0.25) *
                                                             mat4x4::translation(vec3(0.76, 0.54, 1.4))));
        ClippedQuadric* leaf3 = new ClippedQuadric(materials[7]);
        objects.push_back(leaf3->hemisphere(1, 1)->transform(mat4x4::scaling(vec3(0.28, 0.1, 0.01)) *
                                                             mat4x4::rotation(vec3(0, 0, 1), -0.5) *
                                                             mat4x4::translation(vec3(0.32, 0.7, 1.4))));
        ClippedQuadric* leaf4 = new ClippedQuadric(materials[7]);
        objects.push_back(leaf4->hemisphere(1, 1)->transform(mat4x4::scaling(vec3(0.24, 0.1, 0.01)) *
                                                             mat4x4::rotation(vec3(0, 0, 1), 0.45) *
                                                             mat4x4::translation(vec3(0.7, 0.72, 1.4))));
    }
    
    
    //TRACE FUNCTION
    vec3 trace(const Ray& ray){
        Hit hit = firstIntersect(ray);
        if(hit.t < 0)
            return vec3(0.62, 0.90, 0.99);
//            return vec3 (0.2, 0.7, 0.94);
//            return vec3(0.619, 0.8, 1);
        
        float depth = 5;
        vec3 finalColor = vec3(0, 0, 0);
        vec3 totalColor = vec3(0, 0, 0);
        vec3 shadowRayDir, shadowRayOrigin;
        vec3 reflecRayDir, reflecRayOrigin;
        for (int i = 0; i < lightSources.size(); i++) {
            shadowRayDir = lightSources[i]->getLightDirAt(hit.position);
            shadowRayOrigin = hit.position + shadowRayDir * 0.005;
            Hit s = firstIntersect(Ray(shadowRayOrigin, shadowRayDir));
            if (s.t < 0 || lightSources[i]->getDistanceFrom(hit.position) < s.t) {
                finalColor += hit.material->shade(hit.position, hit.normal, lightSources[i], ray);
            }
            if (depth == 0){
                return vec3(0.62, 0.90, 0.99);
            }
            //My attempt at reflections:
//            if (hit.material->reflective()){
//                vec3 viewDir = -ray.dir;
//                reflecRayDir = hit.normal*(viewDir.dot(hit.normal))*2 - viewDir;
//                reflecRayOrigin = hit.position + reflecRayDir * 0.001;
//                Ray reflecRay = Ray (reflecRayOrigin, reflecRayDir);
//                totalColor = totalColor + trace(reflecRay);
//                depth = depth - 1;
//            }
        }
        return finalColor + totalColor;
    }
    

    
    ~Scene()
    {
        for (std::vector<Material*>::iterator iMaterial = materials.begin(); iMaterial != materials.end(); ++iMaterial)
        	delete *iMaterial;
        for (std::vector<Intersectable*>::iterator iObject = objects.begin(); iObject != objects.end(); ++iObject)
        	delete *iObject;
    }
    
    

    
public:


Camera& getCamera()
{
        return camera;
    }
    
    
    Hit firstIntersect(const Ray& ray)
    {
        Hit firstHit, hit;
        float t = 1000000;
        
        for(int i = 0; i < objects.size(); i++)
        {
            Intersectable *obj = objects[i];
            hit = obj->intersect(ray);
            if ((hit.t < t) && (hit.t > 0)) {
                firstHit = hit;
                t = hit.t;
            }
        }
        
        
        
        return firstHit;
    }
};


Scene scene;



bool computeImage()
{
    static unsigned int iPart = 0;
    
    if(iPart >= 64)
        return false;
    for(int j = iPart; j < windowHeight; j+=64)
    {
        for(int i = 0; i < windowWidth; i++)
        {
            float ndcX = (2.0 * i - windowWidth) / windowWidth;
            float ndcY = (2.0 * j - windowHeight) / windowHeight;
            Camera& camera = scene.getCamera();
            Ray ray = Ray(camera.getEye(), camera.rayDirFromNdc(ndcX, ndcY));
            
            image[j*windowWidth + i] = scene.trace(ray);
        }
    }
    iPart++;
    return true;
}

class Texture {
    unsigned int textureId;
public:
    Texture() {
        glGenTextures(1, &textureId);
        glBindTexture(GL_TEXTURE_2D, textureId);
        
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, windowWidth, windowHeight, 0, GL_RGB, GL_FLOAT, image);
        
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    }
    
    void Bind(unsigned int shader)
    {
        int samplerUnit = 0;
        int location = glGetUniformLocation(shader, "samplerUnit");
        glUniform1i(location, samplerUnit);
        glActiveTexture(GL_TEXTURE0 + samplerUnit);
        glBindTexture(GL_TEXTURE_2D, textureId);
        
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, windowWidth, windowHeight, 0, GL_RGB, GL_FLOAT, image);
    }
};

class TexturedQuad {
    Texture texture;
    unsigned int vao;
    
public:
    TexturedQuad()
    {
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);
        
        unsigned int vbo[2];
        glGenBuffers(2, &vbo[0]);
        
        glBindBuffer(GL_ARRAY_BUFFER, vbo[0]);
        //static float vertexCoords[] = { -1, -1,		1, -1,		1, 1,		-1, 1 };
        static float vertexCoords[] = { -1, -1,		1, -1,		-1, 1,		1, 1 };
        
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertexCoords), vertexCoords, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, NULL);
        
        glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
        //static float vertexTextureCoords[] = { 0, 0,	1, 0,		1, 1,		0, 1 };
        static float vertexTextureCoords[] = { 0, 0,	1, 0,		0, 1,		1, 1 };
        
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertexTextureCoords), vertexTextureCoords, GL_STATIC_DRAW);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, NULL);
    }
    
    void Draw()
    {
        glUseProgram(shaderProgram0);
        
        texture.Bind(shaderProgram0);
        
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glBindVertexArray(vao);
        //glDrawArrays(GL_QUADS, 0, 4);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glDisable(GL_BLEND);
    }
};

TexturedQuad *screen = 0;

void onDisplay( ) {
    glClearColor(0.1f, 0.2f, 0.3f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    if(computeImage())
        glutPostRedisplay();
    
    //glDrawPixels(windowWidth, windowHeight, GL_RGB, GL_FLOAT, image);
    
    screen->Draw();
    
    glutSwapBuffers();
}

void onInitialization()
{
    glViewport(0, 0, windowWidth, windowHeight);
    
    // create vertex shader from string
    unsigned int vertexShader0 = glCreateShader(GL_VERTEX_SHADER);
    if (!vertexShader0) { printf("Error in vertex shader creation\n"); exit(1); }
    
    glShaderSource(vertexShader0, 1, &vertexSource0, NULL);
    glCompileShader(vertexShader0);
    checkShader(vertexShader0, "Vertex shader error");
    
    // create fragment shader from string
    unsigned int fragmentShader0 = glCreateShader(GL_FRAGMENT_SHADER);
    if (!fragmentShader0) { printf("Error in fragment shader creation\n"); exit(1); }
    
    glShaderSource(fragmentShader0, 1, &fragmentSource0, NULL);
    glCompileShader(fragmentShader0);
    checkShader(fragmentShader0, "Fragment shader error");
    
    // attach shaders to a single program
    shaderProgram0 = glCreateProgram();
    if (!shaderProgram0) { printf("Error in shader program creation\n"); exit(1); }
    
    glAttachShader(shaderProgram0, vertexShader0);
    glAttachShader(shaderProgram0, fragmentShader0);
    
    // connect Attrib Array to input variables of the vertex shader
    glBindAttribLocation(shaderProgram0, 0, "vertexPosition"); // vertexPosition gets values from Attrib Array 0
    //glBindAttribLocation(shaderProgram, 1, "vertexColor"); // vertexColor gets values from Attrib Array 1
    glBindAttribLocation(shaderProgram0, 1, "vertexTexCoord"); // vertexColor gets values from Attrib Array 1
    
    // connect the fragmentColor to the frame buffer memory
    glBindFragDataLocation(shaderProgram0, 0, "fragmentColor"); // fragmentColor goes to the frame buffer memory
    
    // program packaging
    glLinkProgram(shaderProgram0);
    checkLinking(shaderProgram0);
    
    for(int i = 0; i < windowWidth * windowHeight; i++) image[i] = vec3(0.0, 0.0, 0.0);
    screen = new TexturedQuad();
}

void onExit() 
{
    delete screen;
    screen = 0;
    glDeleteProgram(shaderProgram0);
    printf("exit");
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
    glutCreateWindow("Ray Casting");
    
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
    
    glViewport(0, 0, windowWidth, windowHeight);
    
    onInitialization();
    
    glutDisplayFunc(onDisplay);                
    
    glutMainLoop();
    
    onExit();
    
    return 1;
}
