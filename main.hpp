#pragma once
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <map>
#include <cassert>
#include <algorithm>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#define NMAX 15
#define NMIN 4
// NO REASON TO EVER CHANGE NMIN

struct Halfedge;
struct Face;
struct Mesh;
struct Vertex;

struct Vertex {
    Halfedge* edge;

    double x; double y;
    bool boundary; // chosen as boundary of the embedding.
    int degree; // not necessarily set.
};

struct Halfedge {
    Halfedge* flip;
    Halfedge* next;
    Vertex* vertex;
    Face* face;

    bool boundary;
};

struct Face {
    Halfedge* edge;

    int signature; // 3 degrees sorted written in base 13. 169*a+13*b+c where a,b,c are sorted degrees.
    int degrees[3]; // redundant data to the signature, but more human readable.
};

struct Mesh {
    vector<Vertex*>* vertices;
    vector<Face*>* faces;
    vector<Halfedge*>* edges;
};

template <class T>
int FindInVector(const vector<T*>* vec, T* obj);
Mesh* NewTetrahedronMesh();
Mesh* CopyMesh(const Mesh* start);
void DeleteMesh(Mesh* mesh);
bool CheckS2Validity(Mesh* mesh);
bool IsIsomorphicMesh(vector<Mesh*>* meshes, Mesh* mesh);
bool IsIsomorphicMesh(Mesh* meshA, Mesh* meshB);
bool HasDegree3Vertex(Mesh* mesh);
bool HasDegree4Vertex(Mesh* mesh);
void RecomputeDegrees(Mesh* mesh);
bool propagate(map<Vertex*, Vertex*>* mapAB, map<Halfedge*, Halfedge*>* mapABE, map<Face*, bool>* AFaces, map<Face*, bool>* BFaces,
    Halfedge* AE, Halfedge* BE, bool GoBack);
bool validateIsomorphism(map<Vertex*, Vertex*>* mapAB, map<Halfedge*, Halfedge*>* mapABE, Mesh* meshA);
template <class Y>
bool mapSetOrVerify(map<Y*, Y*>* map, Y* val1, Y* val2);
void ValidateMeshes(vector<vector<Mesh*>*>* meshes);
void VisualizeMeshes(vector<vector<Mesh*>*>* meshes);
void OutputMeshes(vector<vector<Mesh*>*>* meshes);
void OutputMesh(Mesh* mesh);
void SaveMeshes(vector<vector<Mesh*>*>* meshes);
void SaveMesh(Mesh* mesh);
bool HasValidDegrees(Mesh* mesh);
bool HasLessThanDegreeSix(Mesh* mesh);
