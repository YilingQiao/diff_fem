#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__


template <typename TinyScalar, typename TinyConstants> 
class Transformation {
public:
    Eigen::Matrix<TinyScalar, 3, 3> R;
    Eigen::Matrix<TinyScalar, 3, 1> d;
    Transformation()
        :R(Eigen::Matrix<TinyScalar, 3, 3>::Identity()),d(Eigen::Matrix<TinyScalar, 3, 1>::Zero()){}
    Transformation(const Eigen::Matrix<TinyScalar, 3, 3> &R_, const Eigen::Matrix<TinyScalar, 3, 1> &d_)
        :R(R_),d(d_){}
    Transformation(const Eigen::Matrix<TinyScalar, 3, 1> &d_)
        :R(Eigen::Matrix<TinyScalar, 3, 3>::Identity()),d(d_){}
    Transformation operator*(const Transformation &T) {
        return Transformation(R*T.R, R*T.d+d);
    }
    Transformation linearize() {
        Eigen::JacobiSVD<Eigen::Matrix<TinyScalar, 3, 3>> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
        return Transformation(svd.matrixU()*svd.matrixV().transpose(), d);
    }
    Eigen::Matrix<TinyScalar, 3, 1> operator*(const Eigen::Matrix<TinyScalar, 3, 1> &v) {
        return R*v+d;
    }
};

template <typename TinyScalar, typename TinyConstants> 
class Link {
public:
    std::vector<Link*> sons;
    Eigen::Matrix<TinyScalar, 3, 1> joint;
    Link *parent;
    Transformation<TinyScalar, TinyConstants> prefix, T, dT, Tr, suffix, prevTr;
    enum JointType {PRISMATIC, ROTATIONAL, ALL, ERR};
    Eigen::Matrix<TinyScalar, 3, 1> axis, d0;
    JointType jointType;
    int bodyIdx, rigidoffset, dofIdx, storedofIdx;
    std::vector<int> indices;
    std::vector<Eigen::Matrix<TinyScalar, 3, 1> > restpos;
    TinyScalar magnitude, prevmag;
    Eigen::Matrix<TinyScalar, 3, 3> nx, nnt;
    Link() {
        sons.clear();
        joint = Eigen::Matrix<TinyScalar, 3, 1>::Zero();
        parent = NULL;
        axis = Eigen::Matrix<TinyScalar, 3, 1>::Zero();
        jointType = ERR;
        bodyIdx = rigidoffset = 0;
        prefix = T = dT = Tr = Transformation<TinyScalar, TinyConstants>();
        prevTr = Tr;
        indices.clear();
        magnitude = 0;
        prevmag = magnitude;
    }
    Link(Link *pa, Eigen::Matrix<TinyScalar, 3, 1> j, JointType jt, Eigen::Matrix<TinyScalar, 3, 1> ax, int bid) {
        sons.clear();
        joint = j;
        parent = pa;
        axis = ax;
        jointType = jt;
        bodyIdx = bid;
        prefix = T = dT = Tr = Transformation<TinyScalar, TinyConstants>();
        prevTr = Tr;
        if (pa != NULL)
            pa->sons.push_back(this);
        indices.clear();
        magnitude = 0;
        prevmag = magnitude;
        if (jointType == ROTATIONAL) {
            nx<<
                0,-axis[2],axis[1],
                axis[2],0,-axis[0],
                -axis[1],axis[0],0;
            nnt = axis * axis.transpose();
        }
        d0 = joint;
        if (parent != NULL)
            d0 = d0 - parent->joint;
    }

    void guessNext() {
        if (jointType == ALL) {
            T.d = 2 * Tr.d - prevTr.d;
            T.R = Tr.R * prevTr.R.fullPivLu().solve(Tr.R);
            prevTr = Tr;
        } else if (jointType == ROTATIONAL || jointType == PRISMATIC) {
            TinyScalar tmp = magnitude;
            magnitude = 2 * magnitude - prevmag;
            prevmag = tmp;
        } else
            assert(false);
    }

    void ComputeV0NTr(Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &V0, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &minitX) {
        //compute Tr
        if (jointType == ALL)
            Tr = T.linearize();
        else if (jointType == ROTATIONAL) {
            TinyScalar C = cos(magnitude);
            TinyScalar S = sin(magnitude);
            Tr.R = C * Eigen::Matrix<TinyScalar, 3, 3>::Identity() + S * nx + (1-C) * nnt;
            Tr.d = d0;
        } else if (jointType == PRISMATIC)
            Tr.d = magnitude * axis + d0;
        else
            assert(false);
        //compute prefix
        if (parent != NULL)
            prefix = parent->prefix * parent->Tr;
        //compute V0
        for (int k = 0; k < indices.size(); ++k) {
            int rowIdx = 3 * (rigidoffset + k);
            V0.template segment<3>(rowIdx) = prefix * Tr * (restpos[k] - joint);
            // V0.template segment<3>(rowIdx) = prefix * Tr * (minitX.template segment<3>(indices[k]*3) - joint);
        }
        //recursive
        for (auto son : sons)
            son->ComputeV0NTr(V0, minitX);
    }

    void ComputeTildaV(std::vector<Eigen::Triplet<TinyScalar>> &trips, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &V0, const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &minitX) {
        //compute dT
        if (jointType == ROTATIONAL) {
            TinyScalar C = cos(magnitude);
            TinyScalar S = sin(magnitude);
            dT.R = -S * Eigen::Matrix<TinyScalar, 3, 3>::Identity() + C * nx + S * nnt;
        } else if (jointType == PRISMATIC) {
            dT.R = Eigen::Matrix<TinyScalar, 3, 3>::Zero();
            dT.d = axis;
        }
        //compute suffix
        suffix = Transformation<TinyScalar, TinyConstants>();
        Link *cur = this;
        while (cur->parent != NULL) {
            cur->parent->suffix = cur->Tr * cur->suffix;
            cur = cur->parent;
        }
        //compute tildaV
        for (int k = 0; k < indices.size(); ++k) {
            int rowIdx = 3 * (rigidoffset + k);
            Link *cur = this;
            // Eigen::Matrix<TinyScalar, 3, 1> vinit = 
            //     minitX.template segment<3>(indices[k]*3) - joint;
            Eigen::Matrix<TinyScalar, 3, 1> vinit = restpos[k] - joint;
            while (cur->parent != NULL) {
                int colIdx = cur->dofIdx;
                Eigen::Matrix<TinyScalar, 3, 1> tmp;
                if (cur->jointType == ROTATIONAL)
                    tmp = cur->prefix.R * (cur->dT.R * (cur->suffix * vinit));
                else if (cur->jointType == PRISMATIC)
                    tmp = cur->prefix.R * cur->dT.d;
                else
                    assert(false);
                for (int i = 0; i < 3; ++i)
                    trips.push_back(Eigen::Triplet<TinyScalar>(rowIdx+i, colIdx, tmp[i]));
                cur = cur->parent;
            }
            //cur == root
            int colIdx = cur->dofIdx;
            Eigen::Matrix<TinyScalar, 3, 1> v0 = V0.template segment<3>(rowIdx);
            trips.push_back(Eigen::Triplet<TinyScalar>(rowIdx+0,colIdx+1,v0[2]));
            trips.push_back(Eigen::Triplet<TinyScalar>(rowIdx+0,colIdx+2,-v0[1]));
            trips.push_back(Eigen::Triplet<TinyScalar>(rowIdx+1,colIdx+0,-v0[2]));
            trips.push_back(Eigen::Triplet<TinyScalar>(rowIdx+1,colIdx+2,v0[0]));
            trips.push_back(Eigen::Triplet<TinyScalar>(rowIdx+2,colIdx+0,v0[1]));
            trips.push_back(Eigen::Triplet<TinyScalar>(rowIdx+2,colIdx+1,-v0[0]));
            trips.push_back(Eigen::Triplet<TinyScalar>(rowIdx+0,colIdx+3,1));
            trips.push_back(Eigen::Triplet<TinyScalar>(rowIdx+1,colIdx+4,1));
            trips.push_back(Eigen::Triplet<TinyScalar>(rowIdx+2,colIdx+5,1));
        }
        //recursive
        for (auto son : sons)
            son->ComputeTildaV(trips, V0, minitX);
    }

    void UpdateT(const Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &s) {
        if (jointType == PRISMATIC) {
            magnitude += s[dofIdx];
        }
        else if (jointType == ROTATIONAL) {
            TinyScalar C = cos(magnitude);
            TinyScalar S = sin(magnitude);
            magnitude = atan2(S+C*s[dofIdx], C-S*s[dofIdx]);
            // magnitude = magnitude + s[dofIdx];
            // if (bodyIdx == 21) {
            //  std::cout << "s=" << s[dofIdx] << std::endl;
            //  std::cout << "pa s=" << s[parent->dofIdx] << std::endl;
            //  std::cout << "papa s=" << s[parent->parent->dofIdx] << std::endl;
            //  std::cout << "papapa s=" << s[parent->parent->parent->dofIdx] << std::endl;
            //  std::cout << "magnitude=" << magnitude << std::endl;
            //  std::cout << std::endl;
            // }
        }
        else if (jointType == ALL) {
            Eigen::Matrix<TinyScalar, 3, 1> w = s.template segment<3>(dofIdx);
            Eigen::Matrix<TinyScalar, 3, 1> l = s.template segment<3>(dofIdx+3);
            Eigen::Matrix<TinyScalar, 3, 3> tmp;
            tmp <<
                1,-w[2],w[1],
                w[2],1,-w[0],
                -w[1],w[0],1;
            T.R = tmp * Tr.R;
            T.d = tmp * Tr.d + l;
        }
        else
            assert(false);
    }

    void WriteT(Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &q, Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &qd) { // qd is last frame q, not velocity!
        if (jointType == PRISMATIC || jointType == ROTATIONAL) {
            q[storedofIdx] = magnitude;
            qd[storedofIdx] = prevmag;
        }
        else if (jointType == ALL) {
            for (int i = 0; i < 9; i++) 
                q[storedofIdx+i] = T.R(i/3,i%3);
            for (int i = 0; i < 3; i++) 
                q[storedofIdx+i+9] = T.d[i];
            for (int i = 0; i < 9; i++) 
                qd[storedofIdx+i] = prevTr.R(i/3,i%3);
            for (int i = 0; i < 3; i++) 
                qd[storedofIdx+i+9] = prevTr.d[i];
        }
        else
            assert(false);
    }

    void ReadT(Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &q, Eigen::Matrix<TinyScalar, Eigen::Dynamic, 1> &qd) {
        if (jointType == PRISMATIC || jointType == ROTATIONAL) {
            magnitude = q[storedofIdx];
            prevmag = qd[storedofIdx];
        }
        else if (jointType == ALL) {
            for (int i = 0; i < 9; i++) 
                T.R(i/3,i%3) = q[storedofIdx+i];
            for (int i = 0; i < 3; i++) 
                T.d[i] = q[storedofIdx+i+9];
            for (int i = 0; i < 9; i++) 
                prevTr.R(i/3,i%3) = qd[storedofIdx+i];
            for (int i = 0; i < 3; i++) 
                prevTr.d[i] = qd[storedofIdx+i+9];
        }
        else
            assert(false);
    }
};

#endif