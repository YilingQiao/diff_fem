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

#ifndef COLLISIONUTIL_HPP
#define COLLISIONUTIL_HPP

#include "bvh.hpp"

typedef DeformBVHNode BVHNode;
typedef DeformBVHTree BVHTree;

struct AccelStruct {
    BVHTree tree;
    BVHNode *root;
    std::vector<BVHNode*> leaves;
    AccelStruct (const ArcsimMesh &mesh, bool ccd);
};

void update_accel_struct (AccelStruct &acc);

void mark_all_inactive (AccelStruct &acc);
void mark_active (AccelStruct &acc, const Face *face);

// callback must be safe to parallelize via OpenMP
typedef void (*BVHCallback) (Face *face0, Face *face1, std::vector<std::pair<Face*,Face*>> *faces);

// void for_overlapping_faces (BVHNode *node, float thickness,
//                             BVHCallback callback);
// void for_overlapping_faces (BVHNode *node0, BVHNode *node1, float thickness,
//                             BVHCallback callback);
void for_overlapping_faces (const std::vector<AccelStruct*> &accs,
                            const std::vector<AccelStruct*> &obs_accs,
                            double thickness,
                            std::vector<std::pair<Face*,Face*>> *faces,
                            std::vector<std::pair<Face*,Face*>> *obs_faces,
                            BVHCallback callback,
                            bool parallel=true);
// void for_faces_overlapping_obstacles (const std::vector<AccelStruct*> &accs,
//                                       const std::vector<AccelStruct*> &obs_accs,
//                                       double thickness, BVHCallback callback,
//                                       bool parallel=true);

std::vector<AccelStruct*> create_accel_structs
    (const std::vector<ArcsimMesh*> &meshes, bool ccd);
void destroy_accel_structs (std::vector<AccelStruct*> &accs);

// find index of mesh containing specified element
template <typename Prim>
int find_mesh (const Prim *p, const std::vector<ArcsimMesh*> &meshes);

#endif
