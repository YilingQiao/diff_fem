/*
  Copyright ©2013 The Regents of the University of California
  (Regents). All Rights Reserved. Permission to use, copy, modify, and
  distribute this software and its documentation for educational,
  research, and not-for-profit purposes, without fee and without a
  signed licensing agreement, is hereby granted, provided that the
  above copyright notice, this paragraph and the following two
  paragraphs appear in all copies, modifications, and
  distributions. Contact The Office of Technology Licensing, UC
  Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA 94720-1620,
  (510) 643-7201, for commercial licensing opportunities.

  IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT,
  INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
  LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS
  DOCUMENTATION, EVEN IF REGENTS HAS BEEN ADVISED OF THE POSSIBILITY
  OF SUCH DAMAGE.

  REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING
  DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS
  IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/

#include "collisionutil.hpp"
#include "mesh.hpp"
#include "bvh.hpp"
#include <algorithm>
#include <fstream>
#include <set>
// #include <omp.h>
using namespace std;


struct Impact {
    enum Type {VF, EE} type;
    double t;
    Node *nodes[4];
    double w[4], d;
    Vec3 n, contactPoint;
    Impact () {}
    Impact (Type type, Node *n0, Node *n1, Node *n2,
            Node *n3): type(type) {
        nodes[0] = (Node*)n0;
        nodes[1] = (Node*)n1;
        nodes[2] = (Node*)n2;
        nodes[3] = (Node*)n3;
    }
    const static double constexpr thickness = 1e-4;
};



pair<vector<Impact>,vector<Impact>> find_impacts (vector<AccelStruct*> &acc,
                             vector<AccelStruct*> &obs_accs);

pair<vector<Impact>, vector<Impact>> collision_detection (ArcsimMesh *mesh, ArcsimMesh *obs_mesh);


// Impacts


void store_face_impacts (Face *face0, Face *face1, vector<pair<Face*,Face*>> *faces);

bool vf_collision_test (Node *vert, Face *face, Impact &impact, double mult);
// bool ee_collision_test (const Edge *edge0, const Edge *edge1, Impact &impact);

bool collision_test (Impact::Type type, Node *node0, Node *node1,
                     Node *node2, Node *node3, Impact &impact, double mult);

int solve_cubic (double a3, double a2, double a1, double a0, double t[3]);

Vec3 pos (Node *node, double t);

// Solving cubic equations

double newtons_method (double a, double b, double c, double d, double x0,
                       int init_dir);

// solves a x^3 + b x^2 + c x + d == 0
int solve_cubic (double a, double b, double c, double d, double x[3]);

