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

#include "collision.hpp"

pair<vector<Impact>, vector<Impact>> collision_detection (ArcsimMesh *mesh, ArcsimMesh *obs_mesh) {
    vector<ArcsimMesh*> meshes = {mesh};
    vector<ArcsimMesh*> obs_meshes = {obs_mesh};
    // std::cout << "mesh size " << mesh->nodes.size() << " " << mesh->faces.size() << std::endl;
    // std::cout << "obs_mesh size " << obs_mesh->nodes.size() << " " << obs_mesh->faces.size() << std::endl;
    vector<AccelStruct*> accs = create_accel_structs(meshes, true),
                         obs_accs = create_accel_structs(obs_meshes, true);
    return find_impacts(accs, obs_accs);
}


// Impacts

void store_face_impacts (Face *face0, Face *face1, vector<pair<Face*,Face*>> *faces) {
    int t = 0;//omp_get_thread_num();
    faces[t].push_back(make_pair(face0, face1));
}

void do_push_back(vector<Impact> &impacts, set<vector<int>> &impset, const Impact &impact) {
    vector<int> id = {};
    for (int i = 0; i < 4; ++i)
        id.push_back(impact.nodes[i]->index);
    if (impset.find(id) != impset.end())
        return;
    impset.insert(id);
    impacts.push_back(impact);
}

void find_face_impacts (Face *face0, Face *face1, vector<Impact> &impacts, set<vector<int>> &impset) {
    Impact impact;
    for (int v = 0; v < 3; v++)
        if (vf_collision_test(face0->n[v], face1, impact, 2))
            do_push_back(impacts, impset, impact);
    for (int v = 0; v < 3; v++)
        if (vf_collision_test(face1->n[v], face0, impact, 2))
            do_push_back(impacts, impset, impact);
    // for (int e0 = 0; e0 < 3; e0++)
    //     for (int e1 = 0; e1 < 3; e1++)
    //         if (ee_collision_test(face0->adje[e0], face1->adje[e1], impact))
    //             impacts.push_back(impact);
}

void find_face_impacts1 (Face *face0, Face *face1, vector<Impact> &impacts, set<vector<int>> &impset) {
    Impact impact;
    for (int v = 0; v < 3; v++)
        if (vf_collision_test(face0->n[v], face1, impact, 1))
            do_push_back(impacts, impset, impact);
    //original collision code not handling obs v and soft f. omitting here as well
    // for (int v = 0; v < 3; v++)
    //     if (vf_collision_test(face1->n[v], face0, impact)) {
    //         Impact imp = impact;
    //         for (int l = 0; l < 3; ++l)
    //             imp.nodes[l+1] = face0->n[l];
    //         imp.nodes[0] = face1->n[v];
    //         imp.n = -impact.n;
    //     //vertex0 should always be soft body's.when it's obs's v and soft body's f, we need to 3x it
    //         for (int k = 0; k < 3; ++k) {
    //             swap(imp.nodes[0], imp.nodes[k+1]);
    //             do_push_back(impacts, impset, imp);
    //             swap(imp.nodes[0], imp.nodes[k+1]);
    //         }
    //     }
    // for (int e0 = 0; e0 < 3; e0++)
    //     for (int e1 = 0; e1 < 3; e1++)
    //         if (ee_collision_test(face0->adje[e0], face1->adje[e1], impact))
    //             impacts.push_back(impact);
}

pair<vector<Impact>,vector<Impact>> find_impacts (vector<AccelStruct*> &accs,
                             vector<AccelStruct*> &obs_accs) {
    static int nthreads = 0;
    vector<pair<Face*,Face*>> *faces = NULL, *obs_faces = NULL;
    if (!faces) {
        nthreads = 1;//omp_get_max_threads();
        faces = new vector<pair<Face*,Face*>>[nthreads];
        obs_faces = new vector<pair<Face*,Face*>>[nthreads];
    }
    for (int t = 0; t < nthreads; t++) {
        faces[t].clear();
        obs_faces[t].clear();
    }
    for_overlapping_faces(accs, obs_accs, Impact::thickness*2, faces, obs_faces, store_face_impacts);
    vector<Impact> impacts, obs_impacts;
    set<vector<int>> impset, obsimpset;
    impacts.clear();
    obs_impacts.clear();
    for (int t = 0; t < nthreads; ++t) {
        for (auto it : faces[t])
            find_face_impacts(it.first, it.second, impacts, impset);
        for (auto it : obs_faces[t])
            find_face_impacts1(it.first, it.second, obs_impacts, obsimpset);
    }
    return make_pair(impacts, obs_impacts);
}

Vec3 pos (Node *node, double t) {
    return node->x0 + t*(node->x - node->x0);}

bool vf_collision_test (Node *node, Face *face, Impact &impact, double mult) {
    if (node == face->n[0]
     || node == face->n[1]
     || node == face->n[2])
        return false;
    if (!overlap(node_box(node, true), face_box(face, true), Impact::thickness*mult))
        return false;
    // if (node->index == 29)
    //     std::cout << "???" << std::endl;
    return collision_test(Impact::VF, node, face->n[0], face->n[1],
                          face->n[2], impact, mult);
}

// bool ee_collision_test (const Edge *edge0, const Edge *edge1, Impact &impact) {
//     if (edge0->n[0] == edge1->n[0] || edge0->n[0] == edge1->n[1]
//         || edge0->n[1] == edge1->n[0] || edge0->n[1] == edge1->n[1])
//         return false;
//     if (!overlap(edge_box(edge0, true), edge_box(edge1, true), Impact::thickness))
//         return false;
//     return collision_test(Impact::EE, edge0->n[0], edge0->n[1],
//                           edge1->n[0], edge1->n[1], impact);
// }


template <typename T> T sgn (const T &x) {return x<0 ? -1 : 1;}
const double infinity = numeric_limits<double>::infinity();

int solve_quadratic (double a, double b, double c, double x[2]) {
    // http://en.wikipedia.org/wiki/Quadratic_formula#Floating_point_implementation
    double d = b*b - 4*a*c;
    if (d < 0) {
        x[0] = -b/(2*a);
        return 0;
    }
    double q = -(b + sgn(b)*sqrt(d))/2;
    int i = 0;
    if (abs(a) > 1e-12*abs(q))
        x[i++] = q/a;
    if (abs(q) > 1e-12*abs(c))
        x[i++] = c/q;
    if (i==2 && x[0] > x[1])
        swap(x[0], x[1]);
    return i;
}

double signed_vf_distance (const Vec3 &x,
                           const Vec3 &y0, const Vec3 &y1, const Vec3 &y2,
                           Vec3 *n, double *w) {
    Vec3 _n; if (!n) n = &_n;
    double _w[4]; if (!w) w = _w;
    *n = cross(normalize(y1-y0), normalize(y2-y0));
    if (norm2(*n) < 1e-6)
        return infinity;
    *n = normalize(*n);
    double h = dot(x-y0, *n);
    double b0 = stp(y1-x, y2-x, *n),
           b1 = stp(y2-x, y0-x, *n),
           b2 = stp(y0-x, y1-x, *n);
    w[0] = 1;
    w[1] = -b0/(b0 + b1 + b2);
    w[2] = -b1/(b0 + b1 + b2);
    w[3] = -b2/(b0 + b1 + b2);
    return h;
}

double signed_ee_distance (const Vec3 &x0, const Vec3 &x1,
                           const Vec3 &y0, const Vec3 &y1,
                           Vec3 *n, double *w) {
    Vec3 _n; if (!n) n = &_n;
    double _w[4]; if (!w) w = _w;
    *n = cross(normalize(x1-x0), normalize(y1-y0));
    if (norm2(*n) < 1e-6)
        return infinity;
    *n = normalize(*n);
    double h = dot(x0-y0, *n);
    double a0 = stp(y1-x1, y0-x1, *n), a1 = stp(y0-x0, y1-x0, *n),
           b0 = stp(x0-y1, x1-y1, *n), b1 = stp(x1-y0, x0-y0, *n);
    w[0] = a0/(a0 + a1);
    w[1] = a1/(a0 + a1);
    w[2] = -b0/(b0 + b1);
    w[3] = -b1/(b0 + b1);
    return h;
}

bool collision_test (Impact::Type type, Node *node0, Node *node1,
                     Node *node2, Node *node3, Impact &impact, double mult) {
    impact.type = type;
    impact.nodes[0] = (Node*)node0;
    impact.nodes[1] = (Node*)node1;
    impact.nodes[2] = (Node*)node2;
    impact.nodes[3] = (Node*)node3;
    const Vec3 &x0 = node0->x0, v0 = node0->x - x0;
    Vec3 x1 = node1->x0 - x0, x2 = node2->x0 - x0, x3 = node3->x0 - x0;
    Vec3 v1 = (node1->x - node1->x0) - v0, v2 = (node2->x - node2->x0) - v0,
         v3 = (node3->x - node3->x0) - v0;
    double a0 = stp(x1, x2, x3),
           a1 = stp(v1, x2, x3) + stp(x1, v2, x3) + stp(x1, x2, v3),
           a2 = stp(x1, v2, v3) + stp(v1, x2, v3) + stp(v1, v2, x3),
           a3 = stp(v1, v2, v3);
    double t[10];
    int nsol = solve_cubic(a3, a2, a1, a0, t);
    t[nsol++] = 1; // also check at end of timestep
    for (int i = 0; i < nsol; i++) {
        if (t[i] < 0 || t[i] > 1)
            continue;
        impact.t = t[i];
        Vec3 x0 = pos(node0,t[i]), x1 = pos(node1,t[i]),
             x2 = pos(node2,t[i]), x3 = pos(node3,t[i]);
        Vec3 &n = impact.n;
        double *w = impact.w;
        double d;
        bool inside;
        if (type == Impact::VF) {
            d = signed_vf_distance(x0, x1, x2, x3, &n, w);
            inside = (min(min(-w[1], -w[2]), -w[3]) >= -1e-6);
            impact.contactPoint = x0;
            impact.d = abs(d);
        } else {// Impact::EE
            // d = signed_ee_distance(x0, x1, x2, x3, &n, w);
            // inside = (min(w[0], w[1], -w[2], -w[3]) >= -1e-6);
        }
        if (t[i] == 1) {
            //if thickness too large, need to examine the ball as well
            inside = inside || norm(x0 - x1) < Impact::thickness * mult;
            inside = inside || norm(x0 - x2) < Impact::thickness * mult;
            inside = inside || norm(x0 - x3) < Impact::thickness * mult;
            inside = inside || norm(x0-x1-dot(x0-x1,x2-x1)*(x2-x1)/norm2(x2-x1)) < Impact::thickness * mult;
            inside = inside || norm(x0-x2-dot(x0-x2,x3-x2)*(x3-x2)/norm2(x3-x2)) < Impact::thickness * mult;
            inside = inside || norm(x0-x3-dot(x0-x3,x1-x3)*(x1-x3)/norm2(x1-x3)) < Impact::thickness * mult;
        }
        // if (node0->index == 29 && t[i] == 1) {
        //     for (int i = 0; i < 4; ++i)
        //         std::cout << impact.nodes[i]->index << " ";
        //     std::cout << n << " ";
        //     for (int i = 0; i < 4; ++i)
        //         std::cout << w[i] << " ";
        //     std::cout << x0 << " ";
        //     std::cout << x1 << " ";
        //     std::cout << x3 << " ";
        //     std::cout << x0-x3-dot(x0-x3,normalize(x1-x3))*normalize(x1-x3) << " ";
        //     std::cout << d<<" "<<inside << std::endl;
        // }
        if (dot(n, w[1]*v1 + w[2]*v2 + w[3]*v3) > 0)
            // n = -n;
            n = -cross(normalize(node2->x0-node1->x0), normalize(node3->x0-node1->x0));
        else
            n = cross(normalize(node2->x0-node1->x0), normalize(node3->x0-node1->x0));
        if (abs(d) < Impact::thickness * mult && inside)
            return true;
    }
    return false;
}

// Solving cubic equations

// solves a x^3 + b x^2 + c x + d == 0
int solve_cubic (double a, double b, double c, double d, double x[3]) {
    double xc[2];
    int ncrit = solve_quadratic(3*a, 2*b, c, xc);
    if (ncrit == 0) {
        x[0] = newtons_method(a, b, c, d, xc[0], 0);
        return 1;
    } else if (ncrit == 1) {// cubic is actually quadratic
        return solve_quadratic(b, c, d, x);
    } else {
        double yc[2] = {d + xc[0]*(c + xc[0]*(b + xc[0]*a)),
                        d + xc[1]*(c + xc[1]*(b + xc[1]*a))};
        int i = 0;
        if (yc[0]*a >= 0)
            x[i++] = newtons_method(a, b, c, d, xc[0], -1);
        if (yc[0]*yc[1] <= 0) {
            int closer = abs(yc[0])<abs(yc[1]) ? 0 : 1;
            x[i++] = newtons_method(a, b, c, d, xc[closer], closer==0?1:-1);
        }
        if (yc[1]*a <= 0)
            x[i++] = newtons_method(a, b, c, d, xc[1], 1);
        return i;
    }
}

double newtons_method (double a, double b, double c, double d, double x0,
                       int init_dir) {
    if (init_dir != 0) {
        // quadratic approximation around x0, assuming y' = 0
        double y0 = d + x0*(c + x0*(b + x0*a)),
               ddy0 = 2*b + x0*(6*a);
        x0 += init_dir*sqrt(abs(2*y0/ddy0));
    }
    for (int iter = 0; iter < 100; iter++) {
        double y = d + x0*(c + x0*(b + x0*a));
        double dy = c + x0*(2*b + x0*3*a);
        if (dy == 0)
            return x0;
        double x1 = x0 - y/dy;
        if (abs(x0 - x1) < 1e-6)
            return x0;
        x0 = x1;
    }
    return x0;
}
