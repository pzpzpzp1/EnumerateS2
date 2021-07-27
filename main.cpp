#include "main.hpp"

int _BASE = NMAX + 1;
int main(int argc, char** argv)
{
    map<Halfedge*, bool>* reachedEdges = new map<Halfedge*, bool>();
    vector<vector<Mesh*>*>* meshes = new vector<vector<Mesh*>*>();
    for (int i = 0; i <= NMAX; i++) {
        meshes->push_back(new vector<Mesh*>());
    }
    Mesh* mesh = NewTetrahedronMesh();
    RecomputeDegrees(mesh);

    (*meshes)[NMIN]->push_back(mesh);

    cout << "Mesh list made.\n";

    // Start iteration over number of vertices.
    for (int v = NMIN; v < NMAX; v++)
    {
        int numMeshes = (*meshes)[v]->size();
        // Start iteration over meshes of v vertices.
        for (int i = 0; i < numMeshes; i++)
        {
            if (i % 10 == 0) { cout << "Mesh i:v " << i << ":" << v << "nummeshes: " << numMeshes << "\n"; }
            // the mesh to act on. do not modify mesh directly.
            mesh = (*(*meshes)[v])[i];
            RecomputeDegrees(mesh);

            // case 1. loop on faces. apply subdivision.
            for (int f = 0; f < mesh->faces->size(); f++)
            {
                //cout << "case 1 iteration f:" << f << "\n";
                Mesh* copiedMesh = CopyMesh(mesh);
                // create new features
                Vertex* newV = new Vertex();
                Halfedge* newE1 = new Halfedge();
                Halfedge* newE2 = new Halfedge();
                Halfedge* newE3 = new Halfedge();
                Halfedge* newFE1 = new Halfedge();
                Halfedge* newFE2 = new Halfedge();
                Halfedge* newFE3 = new Halfedge();
                Face* newF1 = new Face();
                Face* newF2 = new Face();

                copiedMesh->vertices->push_back(newV);
                copiedMesh->edges->push_back(newE1);
                copiedMesh->edges->push_back(newE2);
                copiedMesh->edges->push_back(newE3);
                copiedMesh->edges->push_back(newFE1);
                copiedMesh->edges->push_back(newFE2);
                copiedMesh->edges->push_back(newFE3);
                copiedMesh->faces->push_back(newF1);
                copiedMesh->faces->push_back(newF2);

                // label current features
                Face* cFace = copiedMesh->faces->at(f);
                Halfedge* e0 = cFace->edge;
                Halfedge* e1 = e0->next;
                Halfedge* e2 = e1->next;

                Vertex* v0 = e0->vertex;
                Vertex* v1 = e1->vertex;
                Vertex* v2 = e2->vertex;

                // link old and new features
                newF1->edge = e2;
                newF2->edge = e1;

                newV->edge = newFE1; // new vert points out.

                newE1->face = newF1;
                newE2->face = newF2;
                newE3->face = cFace;
                newFE1->face = cFace;
                newFE2->face = newF1;
                newFE3->face = newF2;
                e1->face = newF2;
                e2->face = newF1;

                e0->next = newE3;
                e1->next = newE2;
                e2->next = newE1;
                newE3->next = newFE1;
                newE2->next = newFE3;
                newE1->next = newFE2;
                newFE1->next = e0;
                newFE3->next = e1;
                newFE2->next = e2;

                newE1->flip = newFE1;
                newE3->flip = newFE3;
                newE2->flip = newFE2;
                newFE1->flip = newE1;
                newFE2->flip = newE2;
                newFE3->flip = newE3;

                newE1->vertex = newV;
                newE2->vertex = newV;
                newE3->vertex = newV;
                newFE1->vertex = v2;
                newFE2->vertex = v1;
                newFE3->vertex = v0;

                CheckS2Validity(copiedMesh);
                // add to list.
                RecomputeDegrees(mesh);
                (meshes->at(v + 1))->push_back(copiedMesh);
            }

            // case 2. find back to back triangles by loop on edges. apply subdivision.
            reachedEdges->clear();
            for (int e = 0; e < mesh->edges->size(); e++)
            {
                //cout << "case 2 iteration e:" << e << "\n";


                Halfedge* et = mesh->edges->at(e);
                // if half edge represents an edge we iterated on already, then skip it.
                if (reachedEdges->find(et) != reachedEdges->end()
                    && reachedEdges->find(et->flip) != reachedEdges->end()) {
                    continue;
                }
                Mesh* copiedMesh = CopyMesh(mesh);

                // add new features
                Vertex* v4 = new Vertex();
                Halfedge* e6 = new Halfedge();
                Halfedge* e7 = new Halfedge();
                Halfedge* e8 = new Halfedge();
                Halfedge* e9 = new Halfedge();
                Halfedge* e10 = new Halfedge();
                Halfedge* e11 = new Halfedge();
                Face* f2 = new Face();
                Face* f3 = new Face();

                copiedMesh->vertices->push_back(v4);
                copiedMesh->edges->push_back(e6);
                copiedMesh->edges->push_back(e7);
                copiedMesh->edges->push_back(e8);
                copiedMesh->edges->push_back(e9);
                copiedMesh->edges->push_back(e10);
                copiedMesh->edges->push_back(e11);
                copiedMesh->faces->push_back(f2);
                copiedMesh->faces->push_back(f3);

                //current features.
                Halfedge* e0 = copiedMesh->edges->at(e);
                Halfedge* e1 = e0->next;
                //cout << "chkpnt1 " << e1 << endl;
                //cout << "chkpnt1 " << mesh->edges->at(0)->next << endl;
                Halfedge* e2 = e1->next;
                Halfedge* e3 = e0->flip;
                Halfedge* e4 = e3->next;
                Halfedge* e5 = e4->next;
                Face* f0 = e0->face;
                Face* f1 = e3->face;
                Vertex* v0 = e0->vertex;
                Vertex* v1 = e1->vertex;
                Vertex* v2 = e2->vertex;
                Vertex* v3 = e4->vertex;

                //link current features.
                v4->edge = e9;
                v0->edge = e6; v1->edge = e2; v2->edge = e4; v3->edge = e5;

                f2->edge = e10;
                f3->edge = e7;
                f0->edge = e2;
                f1->edge = e4;


                e6->flip = e9;
                e10->flip = e11;
                e7->flip = e8;
                e9->flip = e6;
                e11->flip = e10;
                e8->flip = e7;

                // face, vert, next
                e5->face = f3; e1->face = f2;
                e11->face = f0; e8->face = f1;
                e10->face = f2;
                e7->face = f3;
                e6->face = f3;
                e9->face = f2;

                e0->vertex = v4;
                e6->vertex = v4;
                e7->vertex = v3;
                e8->vertex = v4;
                e9->vertex = v0;
                e10->vertex = v4;
                e11->vertex = v1;

                e0->next = e11;
                e1->next = e10;
                e4->next = e8;
                e5->next = e6;
                e6->next = e7;
                e7->next = e5;
                e8->next = e3;
                e9->next = e1;
                e10->next = e9;
                e11->next = e2;

                //cout << "bloop" << endl;
                CheckS2Validity(copiedMesh);
                //cout << "bloopend" << endl;

                // if there's a degree 3 vertex then this mesh could have been generated with case 1.
                // so it's guaranteed to be isomorphic to something else already.
                if (HasDegree3Vertex(copiedMesh)) {
                    DeleteMesh(copiedMesh);
                    continue;
                }

                // add to list.
                RecomputeDegrees(mesh);
                (meshes->at(v + 1))->push_back(copiedMesh);
            }

            // case 3. find edges with degree 6+ vertices and form pentagons with the right degrees. apply subdivision.
            for (int e = 0; e < mesh->edges->size(); e++)
            {
                Halfedge* e0 = mesh->edges->at(e);
                Vertex* v0 = e0->flip->vertex;
                if (v0->degree < 6) { continue; }
                if (e0->vertex->degree < 4) { continue; }
                if (e0->next->vertex->degree < 5) { continue; }
                if (e0->next->next->flip->next->vertex->degree < 5) { continue; }
                if (e0->next->next->flip->next->next->flip->vertex->degree < 4) { continue; }

                // get current features
                Mesh* copiedMesh = CopyMesh(mesh);
                e0 = copiedMesh->edges->at(e);
                Halfedge* e1 = e0->next;
                Halfedge* e2 = e1->next;
                Halfedge* e3 = e2->flip;
                Halfedge* e4 = e3->next;
                Halfedge* e5 = e4->next;
                Halfedge* e6 = e5->flip;
                Halfedge* e7 = e6->next;
                Halfedge* e8 = e7->next;
                Halfedge* e9 = e8->flip;
                Halfedge* e10 = e0->flip;
                v0 = e10->vertex;
                Vertex* v1 = e9->vertex;
                Vertex* v2 = e6->vertex;
                Vertex* v3 = e3->vertex;
                Vertex* v4 = e0->vertex;
                Face* f0 = e0->face;
                Face* f1 = e3->face;
                Face* f2 = e4->face;

                //Make new features.
                Vertex* v5 = new Vertex();
                Halfedge* e11 = new Halfedge();
                Halfedge* e12 = new Halfedge();
                Halfedge* e13 = new Halfedge();
                Halfedge* e14 = new Halfedge();
                Halfedge* e15 = new Halfedge();
                Halfedge* e16 = new Halfedge();
                Face* f3 = new Face();
                Face* f4 = new Face();

                copiedMesh->edges->push_back(e11);
                copiedMesh->edges->push_back(e12);
                copiedMesh->edges->push_back(e13);
                copiedMesh->edges->push_back(e14);
                copiedMesh->edges->push_back(e15);
                copiedMesh->edges->push_back(e16);
                copiedMesh->vertices->push_back(v5);
                copiedMesh->faces->push_back(f3);
                copiedMesh->faces->push_back(f4);


                // reconnect structure.
                v0->edge = e9;
                v1->edge = e8;
                v2->edge = e5;
                v3->edge = e2;
                v4->edge = e10;
                v5->edge = e0;

                f3->edge = e11;
                f4->edge = e15;

                e0->flip = e14;
                e2->vertex = v5;
                e5->vertex = v5;
                e8->vertex = v5;
                e8->flip = e11;
                e9->flip = e12;
                e10->flip = e16;

                e11->vertex = v1;
                e11->next = e12;
                e11->flip = e8;
                e11->face = f3;

                e12->vertex = v0;
                e12->next = e13;
                e12->flip = e9;
                e12->face = f3;

                e13->vertex = v5;
                e13->next = e11;
                e13->flip = e15;
                e13->face = f3;

                e14->vertex = v5;
                e14->next = e15;
                e14->flip = e0;
                e14->face = f4;

                e15->vertex = v0;
                e15->next = e16;
                e15->flip = e13;
                e15->face = f4;

                e16->vertex = v4;
                e16->next = e14;
                e16->flip = e10;
                e16->face = f4;

                if (HasDegree3Vertex(copiedMesh) || HasDegree4Vertex(copiedMesh)) {
                    DeleteMesh(copiedMesh);
                    continue;
                }

                // add to list.
                RecomputeDegrees(mesh);
                (meshes->at(v + 1))->push_back(copiedMesh);
            }
        }

        // test isomorphism.
        vector<Mesh*>* reducedMeshes = new vector<Mesh*>();
        for (int i = 0; i < meshes->at(v + 1)->size(); i++) {
            if (i % 100 == 0) { cout << "isometries i:n " << i << " " << meshes->at(v + 1)->size() << endl; }

            Mesh* mesh = meshes->at(v + 1)->at(i);
            if ((HasLessThanDegreeSix(mesh) || v + 1 >= 11) && !IsIsomorphicMesh(reducedMeshes, mesh))
            {
                reducedMeshes->push_back(mesh);
            }
            else {
                DeleteMesh(mesh);
            }
        }
        delete meshes->at(v + 1);
        meshes->at(v + 1) = reducedMeshes;

    }

    // remove invalid degreed meshes
    for (int v = NMIN; v < NMAX; v++)
    {
        vector<Mesh*>* reducedMeshes = new vector<Mesh*>();
        for (int i = 0; i < meshes->at(v + 1)->size(); i++) {
            //if(i%100==0){cout<<"isometries i:n " << i << " " << meshes->at(v+1)->size() << endl;}

            Mesh* mesh = meshes->at(v + 1)->at(i);
            if (HasValidDegrees(mesh))
            {
                reducedMeshes->push_back(mesh);
            }
            else {
                DeleteMesh(mesh);
            }
        }
        delete meshes->at(v + 1);
        meshes->at(v + 1) = reducedMeshes;
    }

    ValidateMeshes(meshes);
    for (int i = NMIN; i <= NMAX; i++) { cout << " v:M " << i << ":" << meshes->at(i)->size() << endl; }


    //VisualizeMeshes(meshes);
    //OutputMeshes(meshes);
    SaveMeshes(meshes);


    return 0;
}


bool HasLessThanDegreeSix(Mesh* mesh) {
    RecomputeDegrees(mesh);
    for (int i = 0; i < mesh->vertices->size(); i++) {
        int d = mesh->vertices->at(i)->degree;
        if (d > 5) { return false; }
    }
    return true;
}


bool HasValidDegrees(Mesh* mesh) {
    RecomputeDegrees(mesh);
    for (int i = 0; i < mesh->vertices->size(); i++) {
        Vertex* vert = mesh->vertices->at(i);
        int d = vert->degree;
        if ((!vert->boundary && d != 3 && d != 4 && d != 5) ||
            (vert->boundary && d != 2 && d != 3 && d != 4)
            ) {
            return false;
        }
    }
    return true;
}

// compute tutte
void VisualizeMeshes(vector<vector<Mesh*>*>* meshes) {
    map<Vertex*, int>* vertlabel = new map<Vertex*, int>();
    map<Halfedge*, int>* edgelabel = new map<Halfedge*, int>();
    map<Face*, int>* facelabel = new map<Face*, int>();
    for (int v = NMIN; v <= NMAX; v++)
    {
        for (int i = 0; i < meshes->at(v)->size(); i++) {
            Mesh* mesh = meshes->at(v)->at(i);
            // clear boundary flag
            for (int v2 = 0; v2 < mesh->vertices->size(); v2++) {
                mesh->vertices->at(v2)->boundary = false;
            }
            //label boundary face.
            Face* bFace = mesh->faces->at(0);
            bFace->edge->vertex->x = 1; bFace->edge->vertex->y = sqrt(3) / 2; bFace->edge->vertex->boundary = true;
            bFace->edge->next->vertex->x = -1; bFace->edge->next->vertex->y = sqrt(3) / 2; bFace->edge->next->vertex->boundary = true;
            bFace->edge->next->next->vertex->x = 0; bFace->edge->next->next->vertex->y = -sqrt(3) / 2; bFace->edge->next->next->vertex->boundary = true;

            // label remaining vertices
            vertlabel->clear();
            int labelind = 0;
            for (int v2 = 0; v2 < mesh->vertices->size(); v2++) {
                if (mesh->vertices->at(v2)->boundary) { continue; }
                (*vertlabel)[mesh->vertices->at(v2)] = labelind; labelind++;
            }
            assert(labelind == mesh->vertices->size() - 3);

            int V = mesh->vertices->size() - 3;
            // populate matrix for problem solving
            MatrixXd m = MatrixXd::Zero(V, V); VectorXd X = VectorXd::Zero(V); VectorXd Y = VectorXd::Zero(V);
            for (int v = 0; v < V + 3; v++) {
                Vertex* vert = mesh->vertices->at(v);
                if (vert->boundary) { continue; }

                m(vertlabel->at(vert), vertlabel->at(vert)) = -1;
                int kick = true;
                for (Halfedge* neighbor = vert->edge; neighbor != vert->edge || kick; neighbor = neighbor->flip->next) {
                    kick = false;
                    Vertex* nvert = neighbor->vertex;
                    if (nvert->boundary) {
                        X(vertlabel->at(vert)) -= nvert->x / vert->degree;
                        Y(vertlabel->at(vert)) -= nvert->y / vert->degree;
                    }
                    else {
                        m(vertlabel->at(vert), vertlabel->at(nvert)) = 1.0 / vert->degree;
                        //cout << "REACHED" << endl;
                    }
                }
            }
            VectorXd r1 = m.inverse() * X;
            VectorXd r2 = m.inverse() * Y;

            cout << "m" << m << endl;
            cout << "X" << X << endl;
            cout << "Y" << Y << endl;

            //populate solution
            for (int v = 0; v < V; v++) {
                Vertex* vert = mesh->vertices->at(v);
                if (vert->boundary) { continue; }
                vert->x = r1(vertlabel->at(vert));
                vert->y = r2(vertlabel->at(vert));
            }
        }
    }
    delete vertlabel, edgelabel, facelabel;
}
void OutputMeshes(vector<vector<Mesh*>*>* meshes) {
    map<Vertex*, int>* vertlabel = new map<Vertex*, int>();
    map<Halfedge*, int>* edgelabel = new map<Halfedge*, int>();
    map<Face*, int>* facelabel = new map<Face*, int>();
    for (int v = NMIN; v <= NMAX; v++)
    {
        for (int i = 0; i < meshes->at(v)->size(); i++) {
            Mesh* mesh = meshes->at(v)->at(i);
            vertlabel->clear(); edgelabel->clear(); facelabel->clear();
            cout << "Visualization of mesh v:i " << v << " " << i << endl << endl;
            for (int v2 = 0; v2 < mesh->vertices->size(); v2++) {
                (*vertlabel)[mesh->vertices->at(v2)] = v2;
            }
            for (int v2 = 0; v2 < mesh->edges->size(); v2++) {
                (*edgelabel)[mesh->edges->at(v2)] = v2;
            }
            for (int v2 = 0; v2 < mesh->faces->size(); v2++) {
                (*facelabel)[mesh->faces->at(v2)] = v2;
            }

            for (int v2 = 0; v2 < mesh->faces->size(); v2++) {
                Face* face = mesh->faces->at(v2);
                cout << "    face " << facelabel->at(face) << " edges " <<
                    edgelabel->at(face->edge) << " " <<
                    edgelabel->at(face->edge->next) << " " <<
                    edgelabel->at(face->edge->next->next) << " signature " << face->signature << endl;
            }
            cout << endl;
            for (int v2 = 0; v2 < mesh->edges->size(); v2++) {
                Halfedge* edge = mesh->edges->at(v2);
                cout << "    edge " << edgelabel->at(edge) << " flipedge " <<
                    edgelabel->at(edge->flip) << " " << endl;
            }
            cout << endl;
            for (int v2 = 0; v2 < mesh->vertices->size(); v2++) {
                Vertex* vert = mesh->vertices->at(v2);
                cout << "    vert " << vertlabel->at(vert) << " edge " <<
                    edgelabel->at(vert->edge) << " degree " << vert->degree
                    << " x y " << vert->x << " " << vert->y << endl;
            }
            cout << endl;

        }

    }
    delete vertlabel, edgelabel, facelabel;
}

void OutputMesh(Mesh* mesh) {
    map<Vertex*, int>* vertlabel = new map<Vertex*, int>();
    map<Halfedge*, int>* edgelabel = new map<Halfedge*, int>();
    map<Face*, int>* facelabel = new map<Face*, int>();

    vertlabel->clear(); edgelabel->clear(); facelabel->clear();
    for (int v2 = 0; v2 < mesh->vertices->size(); v2++) {
        (*vertlabel)[mesh->vertices->at(v2)] = v2;
    }
    for (int v2 = 0; v2 < mesh->edges->size(); v2++) {
        (*edgelabel)[mesh->edges->at(v2)] = v2;
    }
    for (int v2 = 0; v2 < mesh->faces->size(); v2++) {
        (*facelabel)[mesh->faces->at(v2)] = v2;
    }

    for (int v2 = 0; v2 < mesh->faces->size(); v2++) {
        Face* face = mesh->faces->at(v2);
        cout << "    face " << facelabel->at(face) << " edges " <<
            edgelabel->at(face->edge) << " " <<
            edgelabel->at(face->edge->next) << " " <<
            edgelabel->at(face->edge->next->next) << " signature " << face->signature << endl;
    }
    cout << endl;
    for (int v2 = 0; v2 < mesh->edges->size(); v2++) {
        Halfedge* edge = mesh->edges->at(v2);
        cout << "    edge " << edgelabel->at(edge) << " flipedge " <<
            edgelabel->at(edge->flip) << " " << edge << endl;
    }
    cout << endl;
    for (int v2 = 0; v2 < mesh->vertices->size(); v2++) {
        Vertex* vert = mesh->vertices->at(v2);
        cout << "    vert " << vertlabel->at(vert) << " edge " <<
            edgelabel->at(vert->edge) << " degree " << vert->degree
            << " x y " << vert->x << " " << vert->y << endl;
    }
    cout << endl;
    delete vertlabel, edgelabel, facelabel;
}

void SaveMesh(Mesh* mesh) {
    map<Vertex*, int>* vertlabel = new map<Vertex*, int>();
    map<Halfedge*, int>* edgelabel = new map<Halfedge*, int>();
    map<Face*, int>* facelabel = new map<Face*, int>();
    const size_t bsize = 100;
    char buffer[bsize];
    FILE* pfile;


    vertlabel->clear(); edgelabel->clear(); facelabel->clear();
    int v = mesh->vertices->size();
    sprintf_s(buffer, bsize, "meshONETEMPV%.2d.txt", v);
    fopen_s(&pfile, buffer, "w");

    for (int v2 = 0; v2 < mesh->vertices->size(); v2++) {
        Vertex* vert = mesh->vertices->at(v2);
        (*vertlabel)[vert] = v2;
    }
    for (int v2 = 0; v2 < mesh->faces->size(); v2++) {
        Face* face = mesh->faces->at(v2);
        (*facelabel)[face] = v2;
    }

    for (int v2 = 0; v2 < mesh->edges->size(); v2++) {
        Halfedge* edge = mesh->edges->at(v2);
        (*edgelabel)[edge] = v2;
    }

    fprintf(pfile, "%.2d %.2d %.3d vfe\n",
        mesh->vertices->size(),
        mesh->faces->size(),
        mesh->edges->size()
    );
    for (int v2 = 0; v2 < mesh->vertices->size(); v2++) {
        Vertex* vert = mesh->vertices->at(v2);
        fprintf(pfile, "v %.2d e %.3d degree %.2d position %f %f\n",
            v2,
            edgelabel->at(vert->edge),
            vert->degree,
            vert->x,
            vert->y);
    }
    for (int v2 = 0; v2 < mesh->faces->size(); v2++) {
        Face* face = mesh->faces->at(v2);
        fprintf(pfile, "f %.2d v %.2d %.2d %.2d e %.3d\n", v2,
            vertlabel->at(face->edge->vertex),
            vertlabel->at(face->edge->next->vertex),
            vertlabel->at(face->edge->next->next->vertex),
            edgelabel->at(face->edge));
    }

    for (int v2 = 0; v2 < mesh->edges->size(); v2++) {
        Halfedge* edge = mesh->edges->at(v2);
        fprintf(pfile, "he %.3d flip %.3d next %.3d face %.2d vertex %.2d \n",
            v2,
            edgelabel->at(edge->flip),
            edgelabel->at(edge->next),
            facelabel->at(edge->face),
            vertlabel->at(edge->vertex));
    }

    fclose(pfile);
}

void SaveMeshes(vector<vector<Mesh*>*>* meshes) {
    map<Vertex*, int>* vertlabel = new map<Vertex*, int>();
    map<Halfedge*, int>* edgelabel = new map<Halfedge*, int>();
    map<Face*, int>* facelabel = new map<Face*, int>();
    const size_t bsize = 100;
    char buffer[bsize];
    FILE* pfile;
    for (int v = NMIN; v <= NMAX; v++) {
        for (int i = 0; i < meshes->at(v)->size(); i++) {
            Mesh* mesh = meshes->at(v)->at(i);
            vertlabel->clear(); edgelabel->clear(); facelabel->clear();

            sprintf_s(buffer, bsize, "output/meshV%.2d_I%.2d.txt", v, i);
            fopen_s(&pfile, buffer, "w");

            for (int v2 = 0; v2 < mesh->vertices->size(); v2++) {
                Vertex* vert = mesh->vertices->at(v2);
                (*vertlabel)[vert] = v2;
            }
            for (int v2 = 0; v2 < mesh->faces->size(); v2++) {
                Face* face = mesh->faces->at(v2);
                (*facelabel)[face] = v2;
            }

            for (int v2 = 0; v2 < mesh->edges->size(); v2++) {
                Halfedge* edge = mesh->edges->at(v2);
                (*edgelabel)[edge] = v2;
            }

            fprintf(pfile, "%.2d %.2d %.3d vfe\n",
                mesh->vertices->size(),
                mesh->faces->size(),
                mesh->edges->size()
            );
            for (int v2 = 0; v2 < mesh->vertices->size(); v2++) {
                Vertex* vert = mesh->vertices->at(v2);
                fprintf(pfile, "v %.2d e %.3d degree %.2d position %f %f\n",
                    v2,
                    edgelabel->at(vert->edge),
                    vert->degree,
                    vert->x,
                    vert->y);
            }
            for (int v2 = 0; v2 < mesh->faces->size(); v2++) {
                Face* face = mesh->faces->at(v2);
                fprintf(pfile, "f %.2d v %.2d %.2d %.2d e %.3d\n", v2,
                    vertlabel->at(face->edge->vertex),
                    vertlabel->at(face->edge->next->vertex),
                    vertlabel->at(face->edge->next->next->vertex),
                    edgelabel->at(face->edge));
            }

            for (int v2 = 0; v2 < mesh->edges->size(); v2++) {
                Halfedge* edge = mesh->edges->at(v2);
                fprintf(pfile, "he %.3d flip %.3d next %.3d face %.2d vertex %.2d \n",
                    v2,
                    edgelabel->at(edge->flip),
                    edgelabel->at(edge->next),
                    facelabel->at(edge->face),
                    vertlabel->at(edge->vertex));
            }

            fclose(pfile);
        }
    }
}
void ValidateMeshes(vector<vector<Mesh*>*>* meshes) {
    for (int v = NMIN; v <= NMAX; v++)
    {
        for (int i = 0; i < meshes->at(v)->size(); i++) {
            Mesh* mesh = meshes->at(v)->at(i);
            assert(mesh->vertices->size() == v);
            CheckS2Validity(mesh);
        }

    }
}
Mesh* NewTetrahedronMesh()
{
    Mesh* mesh = new Mesh();
    mesh->vertices = new vector<Vertex*>(); mesh->vertices->clear();
    mesh->edges = new vector<Halfedge*>(); mesh->edges->clear();
    mesh->faces = new vector<Face*>(); mesh->faces->clear();

    Vertex* v0 = new Vertex();
    Vertex* v1 = new Vertex();
    Vertex* v2 = new Vertex();
    Vertex* v3 = new Vertex();
    mesh->vertices->push_back(v0);
    mesh->vertices->push_back(v1);
    mesh->vertices->push_back(v2);
    mesh->vertices->push_back(v3);

    Halfedge* he0 = new Halfedge();
    Halfedge* he1 = new Halfedge();
    Halfedge* he2 = new Halfedge();
    Halfedge* he3 = new Halfedge();
    Halfedge* he4 = new Halfedge();
    Halfedge* he5 = new Halfedge();
    Halfedge* he6 = new Halfedge();
    Halfedge* he7 = new Halfedge();
    Halfedge* he8 = new Halfedge();
    Halfedge* he9 = new Halfedge();
    Halfedge* he10 = new Halfedge();
    Halfedge* he11 = new Halfedge();
    mesh->edges->push_back(he0);
    mesh->edges->push_back(he1);
    mesh->edges->push_back(he2);
    mesh->edges->push_back(he3);
    mesh->edges->push_back(he4);
    mesh->edges->push_back(he5);
    mesh->edges->push_back(he6);
    mesh->edges->push_back(he7);
    mesh->edges->push_back(he8);
    mesh->edges->push_back(he9);
    mesh->edges->push_back(he10);
    mesh->edges->push_back(he11);


    Face* f0 = new Face();
    Face* f1 = new Face();
    Face* f2 = new Face();
    Face* f3 = new Face();
    mesh->faces->push_back(f0);
    mesh->faces->push_back(f1);
    mesh->faces->push_back(f2);
    mesh->faces->push_back(f3);

    f0->edge = he0;
    f1->edge = he3;
    f2->edge = he6;
    f3->edge = he9;

    v0->edge = he0;
    v1->edge = he1;
    v2->edge = he2;
    v3->edge = he7;

    he0->face = f0;
    he1->face = f0;
    he2->face = f0;
    he3->face = f1;
    he4->face = f1;
    he5->face = f1;
    he6->face = f2;
    he7->face = f2;
    he8->face = f2;
    he9->face = f3;
    he10->face = f3;
    he11->face = f3;

    he0->vertex = v1;
    he1->vertex = v2;
    he2->vertex = v0;
    he3->vertex = v2;
    he4->vertex = v3;
    he5->vertex = v0;
    he6->vertex = v3;
    he7->vertex = v1;
    he8->vertex = v0;
    he9->vertex = v1;
    he10->vertex = v3;
    he11->vertex = v2;

    he0->next = he1;
    he1->next = he2;
    he2->next = he0;
    he3->next = he4;
    he4->next = he5;
    he5->next = he3;
    he6->next = he7;
    he7->next = he8;
    he8->next = he6;
    he9->next = he10;
    he10->next = he11;
    he11->next = he9;

    he0->flip = he8;
    he1->flip = he9;
    he2->flip = he3;
    he3->flip = he2;
    he4->flip = he11;
    he5->flip = he6;
    he6->flip = he5;
    he7->flip = he10;
    he8->flip = he0;
    he9->flip = he1;
    he10->flip = he7;
    he11->flip = he4;

    CheckS2Validity(mesh);

    cout << "Tetrahedron made.\n";

    return mesh;
}
Mesh* CopyMesh(const Mesh* start) {
    //cout << "CopyMesh started. \n";
    Mesh* result = new Mesh();
    result->vertices = new vector<Vertex*>(); result->vertices->clear();
    result->edges = new vector<Halfedge*>(); result->edges->clear();
    result->faces = new vector<Face*>(); result->faces->clear();

    // helper allocations
    map<Vertex*, Vertex*>* vertmap = new map<Vertex*, Vertex*>();
    map<Halfedge*, Halfedge*>* edgemap = new map<Halfedge*, Halfedge*>();
    map<Face*, Face*>* facemap = new map<Face*, Face*>();

    //duplicate verts
    for (int i = 0; i < start->vertices->size(); i++)
    {
        Vertex* v = new Vertex();
        result->vertices->push_back(v);
        (*vertmap)[start->vertices->at(i)] = v;
    }
    for (int i = 0; i < start->edges->size(); i++)
    {
        //cout << "CopyMesh started. edgemap "<< (start->edges->at(i)) << " "<<i <<endl;
        result->edges->push_back(new Halfedge());
        (*edgemap)[(*start->edges)[i]] = (*result->edges)[i];
    }
    for (int i = 0; i < start->faces->size(); i++)
    {
        result->faces->push_back(new Face());
        (*facemap)[(*start->faces)[i]] = (*result->faces)[i];
    }

    // copy structure
    for (int i = 0; i < result->vertices->size(); i++)
    {
        //cout << "CopyMesh started verts. "<<i << start->vertices->at(i)->edge<<endl;
        Vertex* resultVert = result->vertices->at(i);
        resultVert->edge = edgemap->at(start->vertices->at(i)->edge);
        resultVert->degree = start->vertices->at(i)->degree;
    }
    for (int i = 0; i < result->faces->size(); i++)
    {
        //cout << "CopyMesh started faces. "<<i<<start->faces->at(i)->edge<<endl;
        Face* resultFace = result->faces->at(i);
        resultFace->edge = edgemap->at(start->faces->at(i)->edge);
        resultFace->signature = start->faces->at(i)->signature;
        resultFace->degrees[0] = start->faces->at(i)->degrees[0];
        resultFace->degrees[1] = start->faces->at(i)->degrees[1];
        resultFace->degrees[2] = start->faces->at(i)->degrees[2];
    }
    for (int i = 0; i < result->edges->size(); i++)
    {

        Halfedge* resultEdge = result->edges->at(i);//cout << "CopyMesh started faces. 1 "<<i<<start->edges->at(i)->flip<<endl;
        resultEdge->flip = edgemap->at(start->edges->at(i)->flip);//cout << "CopyMesh started faces. 2 "<<i<<start->edges->at(i)->next<<endl;
        resultEdge->next = edgemap->at(start->edges->at(i)->next);//cout << "CopyMesh started faces. 3 "<<i<<start->edges->at(i)->face<<endl;
        resultEdge->face = facemap->at(start->edges->at(i)->face);//cout << "CopyMesh started faces. 4 "<<i<<start->edges->at(i)->vertex<<endl;
        resultEdge->vertex = vertmap->at(start->edges->at(i)->vertex);
        resultEdge->boundary = start->edges->at(i)->boundary;
    }

    // free allocations
    delete vertmap;
    delete edgemap;
    delete facemap;
    CheckS2Validity(result);
    //cout << "CopyMesh ended. \n";

    return result;
}
void DeleteMesh(Mesh* mesh) {
    //cout << "Delete mesh started.\n";
    for (int i = 0; i < mesh->vertices->size(); i++)
    {
        delete mesh->vertices->at(i);
    }
    for (int i = 0; i < mesh->edges->size(); i++)
    {
        delete mesh->edges->at(i);
    }
    for (int i = 0; i < mesh->faces->size(); i++)
    {
        delete mesh->faces->at(i);
    }

    delete mesh->vertices;
    delete mesh->edges;
    delete mesh->faces;
    delete mesh;

    //cout << "Delete mesh ended.\n";
}
bool CheckS2Validity(Mesh* mesh) {
    return true;
    assert(mesh != NULL);
    assert(mesh->vertices != NULL);
    assert(mesh->edges != NULL);
    assert(mesh->faces != NULL);
    int V = mesh->vertices->size();
    int HE = mesh->edges->size();
    int E = HE / 2;
    int F = mesh->faces->size();
    assert(3 * F == 2 * E);
    assert(3 * F == HE);
    assert(V - E + F == 2);

    assert(V - E / 3.0 == 2);
    assert(V - F / 2.0 == 2);



    RecomputeDegrees(mesh);
    //OutputMesh(mesh);
    //RecomputeDegrees(mesh);
    //SaveMesh(mesh);

    double counter = 0;
    double counter2 = 0;
    for (int i = 0; i < mesh->vertices->size(); i++)
    {
        //counter += 2.0/3.0*(1-mesh->vertices->at(i)->degree/4.0);
        counter += (1 - mesh->vertices->at(i)->degree / 6.0);
        counter2 += (1.0 / 3.0 + 2.0 / 3.0 * (1 - mesh->vertices->at(i)->degree / 4.0));
        //cout<<"ctrA "<< counter <<" " << mesh->vertices->at(i)->degree << endl;
        assert(-1 != FindInVector<Halfedge>(mesh->edges, mesh->vertices->at(i)->edge));
    }
    //counter += ((1.0/3.0) * (mesh->vertices->size()));
    //cout<<"ctrA "<< counter << endl;
    //cout<<"ctr "<< counter <<" " << mesh->vertices->size()<< endl;
    assert(counter > 1.999 && counter < 2.0001);
    assert(counter2 > 1.999 && counter2 < 2.0001);

    for (int i = 0; i < mesh->faces->size(); i++)
    {
        Face* face = mesh->faces->at(i);
        assert(-1 != FindInVector<Halfedge>(mesh->edges, face->edge));
        assert(face->edge->face == face);
        assert(face->edge->next->face == face);
        assert(face->edge->next->next->face == face);
    }
    for (int i = 0; i < mesh->edges->size(); i++)
    {
        assert(-1 != FindInVector<Halfedge>(mesh->edges, mesh->edges->at(i)->flip));
        assert(-1 != FindInVector<Halfedge>(mesh->edges, mesh->edges->at(i)->next));
        assert(-1 != FindInVector<Vertex>(mesh->vertices, mesh->edges->at(i)->vertex));
        assert(-1 != FindInVector<Face>(mesh->faces, mesh->edges->at(i)->face));
        assert(mesh->edges->at(i)->next->next->next == mesh->edges->at(i));
    }
}
template <class T>
int FindInVector(const vector<T*>* vec, T* obj) {
    assert(vec != NULL);
    for (int i = 0; i < vec->size(); i++) {
        assert(vec->at(i) != NULL);
        if (vec->at(i) == obj) { return i; }
    }
    return -1;
}
bool IsIsomorphicMesh(vector<Mesh*>* meshes, Mesh* mesh) {

    // computing degrees here at a stage where the mesh will not change again. This also results in meshes having degrees computed.
    RecomputeDegrees(mesh);

    // test if mesh is isomorphic to any in meshes.
    for (int i = 0; i < meshes->size(); i++)
    {
        Mesh* meshA = meshes->at(i);
        if (IsIsomorphicMesh(meshA, mesh)) {
            return true;
        }
    }

    return false;
}
void RecomputeDegrees(Mesh* mesh) {
    for (int v = 0; v < mesh->vertices->size(); v++) {
        Halfedge* end = mesh->vertices->at(v)->edge;
        Halfedge* e = end->flip->next;

        //cout << "end " << end << endl;
        //cout << "e " << e << endl;

        assert(e != end); // degree 2 is not allowed.
        int deg = 2;
        while (e != end) {

            e = e->flip->next;
            //cout << "e " << e << endl;
            deg++;
        }
        //cout << "recompute degrees "<< deg-1 << "vert "<< v <<endl;
        mesh->vertices->at(v)->degree = deg - 1; // returning to end results in an extra degree
    }

    for (int f = 0; f < mesh->faces->size(); f++) {
        Face* face = mesh->faces->at(f);
        int a = face->edge->vertex->degree;
        int b = face->edge->next->vertex->degree;
        int c = face->edge->next->next->vertex->degree;

        int t;

        if (b > a) { t = a; a = b; b = t; }
        if (c > b) { t = b; b = c; c = t; }
        if (b > a) { t = a; a = b; b = t; }

        face->signature = _BASE * _BASE * a + _BASE * b + c;
        face->degrees[0] = a; face->degrees[1] = b; face->degrees[2] = c;
    }
}
bool HasDegree3Vertex(Mesh* mesh) {
    //return false;
    //cout << "testing degree 3 vert" << endl;

    for (int v = 0; v < mesh->vertices->size(); v++) {
        Vertex* vert = mesh->vertices->at(v);
        if (vert->edge->flip->next->flip->next->flip->next == vert->edge) {
            //cout << "deg 3 found" << endl;
            return true;
        }
    }
    //cout << "no deg 3" << endl;
    return false;
}
bool HasDegree4Vertex(Mesh* mesh) {

    for (int v = 0; v < mesh->vertices->size(); v++) {
        Vertex* vert = mesh->vertices->at(v);
        if (vert->edge->flip->next->flip->next->flip->next->flip->next == vert->edge) {
            return true;
        }
    }
    return false;
}
bool OrderFaces(Face* a, Face* b) { return a->signature > b->signature; }
bool OrderVertices(Vertex* a, Vertex* b) { return a->degree > b->degree; }
bool OrderEdges(Halfedge* a, Halfedge* b) { return a->face->signature > b->face->signature; }
bool IsIsomorphicMesh(Mesh* meshA, Mesh* meshB) {
    //return false;

    int VA = meshA->vertices->size();
    int HEA = meshA->edges->size();
    int EA = meshA->edges->size() / 2;
    int FA = meshA->faces->size();
    int VB = meshB->vertices->size();
    int HEB = meshB->edges->size();
    int EB = meshB->edges->size() / 2;
    int FB = meshB->faces->size();

    // if the number of elements differ, then they aren't isomorphic
    if (VA != VB || HEA != HEB || FA != FB) {
        return false;
    }

    // map verts of A to verts of B
    map<Vertex*, Vertex*>* mapAB = new map<Vertex*, Vertex*>();
    map<Halfedge*, Halfedge*>* mapABE = new map<Halfedge*, Halfedge*>();
    map<Face*, bool>* AFaces = new map<Face*, bool>();
    map<Face*, bool>* BFaces = new map<Face*, bool>();

    // get biggest face and map to biggest faces
    sort(meshA->faces->begin(), meshA->faces->end(), OrderFaces);
    sort(meshB->faces->begin(), meshB->faces->end(), OrderFaces);
    /*sort(meshA->vertices->begin(), meshA->vertices->end(), OrderVertices);
    sort(meshB->vertices->begin(), meshB->vertices->end(), OrderVertices);
    sort(meshA->edges->begin(), meshA->edges->end(), OrderEdges);
    sort(meshB->edges->begin(), meshB->edges->end(), OrderEdges);*/

    Face* face = meshB->faces->at(0);
    int signature = meshA->faces->at(0)->signature;
    for (int f = 0; f < meshA->faces->size()
        /*comment this line out to get more brute force.*/
        && signature == meshA->faces->at(f)->signature
        ; f++) {
        Face* faceA = meshA->faces->at(f);

        // if iso between faceA and faceB x3, return true.
        bool kick = true;
        for (Halfedge* E1 = faceA->edge; E1 != faceA->edge || kick; E1 = E1->next) {
            kick = false;
            mapAB->clear(); mapABE->clear();
            AFaces->clear(); BFaces->clear();

            if (propagate(mapAB, mapABE, AFaces, BFaces, E1, face->edge, true)
                && validateIsomorphism(mapAB, mapABE, meshA)) {
                delete mapAB; delete mapABE; delete AFaces; delete BFaces;
                return true;
            };

            /*mapAB->clear(); mapABE->clear();
            AFaces->clear(); BFaces->clear();

            // this doesn't seem to do anything.
            if(propagateFlipped(mapAB,mapABE,AFaces,BFaces,E1,face->edge)
                    && validateIsomorphism(mapAB,mapABE,meshA)){
                delete mapAB; delete mapABE; delete AFaces; delete BFaces;
                return true;
            };*/
        }
    }

    delete mapAB; delete mapABE; delete AFaces; delete BFaces;
    return false;
}

// returns if isomorphism is correct or not.
bool validateIsomorphism(map<Vertex*, Vertex*>* mapAB, map<Halfedge*, Halfedge*>* mapABE, Mesh* meshA) {
    return true;
}

bool propagate(map<Vertex*, Vertex*>* mapAB, map<Halfedge*, Halfedge*>* mapABE,
    map<Face*, bool>* AFaces, map<Face*, bool>* BFaces,
    Halfedge* AE, Halfedge* BE, bool GoBack) {


    bool AFaceFound = AFaces->find(AE->face) != AFaces->end();
    bool BFaceFound = BFaces->find(BE->face) != BFaces->end();

    if (AFaceFound ^ BFaceFound) { return false; } // isomorphism ended in contradiction. give up.
    if (AFaceFound && BFaceFound) { return true; } // continue.
    (*AFaces)[AE->face] = true; (*BFaces)[BE->face] = true;

    if (!mapSetOrVerify<Halfedge>(mapABE, AE, BE)) { return false; };
    if (!mapSetOrVerify<Halfedge>(mapABE, AE->next, BE->next)) { return false; };
    if (!mapSetOrVerify<Halfedge>(mapABE, AE->next->next, BE->next->next)) { return false; };

    if (!mapSetOrVerify<Vertex>(mapAB, AE->vertex, BE->vertex)) { return false; };
    if (!mapSetOrVerify<Vertex>(mapAB, AE->next->vertex, BE->next->vertex)) { return false; };
    if (!mapSetOrVerify<Vertex>(mapAB, AE->next->next->vertex, BE->next->next->vertex)) { return false; };

    if (GoBack)
    {
        if (!propagate(mapAB, mapABE, AFaces, BFaces, AE->flip, BE->flip, false)) { return false; };
    }


    if (!propagate(mapAB, mapABE, AFaces, BFaces, AE->next->flip, BE->next->flip, false)) { return false; };
    if (!propagate(mapAB, mapABE, AFaces, BFaces, AE->next->next->flip, BE->next->next->flip, false)) { return false; };

    return true;
}

// If map contains the matching already, verify that this attempt to set values doesn't contradict.
// otherwise just set the matching. return if there's a contradiction.
template <class Y>
bool mapSetOrVerify(map<Y*, Y*>* map, Y* val1, Y* val2) {
    if (map->find(val1) == map->end())
    {
        (*map)[val1] = val2;
    }
    else {
        return (map->at(val1) == val2);
    }
    return true;
}
