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

#ifndef MESH_HPP
#define MESH_HPP
#include <utility>
#include <Eigen/Core>
#include "vectors.hpp"
#include <vector>


struct Node {
    int label;
    Vec3 x, x0; // position, old (collision-free) position,
    int index; // position in mesh.nodes
    Node () {}
    explicit Node (const Vec3 &x):
        x(x), x0(x) {}
};

struct Face {
    Node* n[3]; // verts
    int index; // position in mesh.faces
    Face () {}
    explicit Face (Node *vert0, Node *vert1, Node *vert2){
        n[0] = vert0;
        n[1] = vert1;
        n[2] = vert2;
    }
};

struct ArcsimMesh {
    std::vector<Node*> nodes;
    std::vector<Face*> faces;

    void Initialize(const std::vector<Eigen::Vector3d> &pos, 
        const std::vector<Eigen::Vector3d> &pos0, 
        const std::vector<Eigen::Vector3i> &fcs) {
        for (int i = 0; i < pos.size(); ++i) {
            Vec3 tmp(pos[i][0],pos[i][1],pos[i][2]);
            nodes.push_back(new Node(tmp));
            nodes.back()->x0 = Vec3(pos0[i][0],pos0[i][1],pos0[i][2]);
            nodes.back()->index = i;
        }
        for (int i = 0; i < fcs.size(); ++i) {
            const Eigen::Vector3i &face = fcs[i];
            faces.push_back(new Face(nodes[face[0]], nodes[face[1]], nodes[face[2]]));
            faces.back()->index = i;
        }
    }
};


#endif
