#pragma once

#include <array>
#include <vector.h>
// #include <cstring>
// #include <algorithm>


class Mat3 {
public:
    Mat3(){
        p.fill(0);
    }
    Mat3(const Mat3 &m){p = m.p;}

    double& operator()(int i, int j) {
        return p[i*3+j];
    }

    Mat3& operator=(const Mat3& m) {
        p = m.p;
        return *this;
    }

    Mat3 operator+(const Mat3& m) const {
        Mat3 ret;
        for (int i=0; i < 9; ++i) ret.p[i] = p[i] + m.p[i];
        return ret;
    }

    Mat3 operator-(const Mat3& m) const {
        Mat3 ret;
        for (int i=0; i < 9; ++i) ret.p[i] = p[i] - m.p[i];
        return ret;
    }

    Mat3 operator*(const Mat3& m) const {
        Mat3 ret;
        for (int i = 0; i < 3; i++){ // cols
            for (int j = 0; j < 3; j++){ // rows
                ret(i, j) = p[j*3+0]*m.p[i*3+j] + p[j*3+1]*m.p[i*3+j] + p[j*3+2]*m.p[i*3+j];
            }
        }
        return ret;
    }

    Vec3 operator*(const Vec3& v) const {
        return Vec3(
            p[0]*v.x + p[1]*v.y + p[2]*v.z,
            p[3]*v.x + p[4]*v.y + p[5]*v.z,
            p[6]*v.x + p[7]*v.y + p[8]*v.z
        );
    }

    Mat3 operator*(const double scalar) const {
        Mat3 ret;
        for (int i=0; i < 9; ++i) ret.p[i] = scalar * p[i];
        return ret;
    }

    Mat3 operator/(const double scalar) const {
        Mat3 ret;
        for (int i=0; i < 9; ++i) ret.p[i] = p[i]/scalar;
        return ret;
    }

    double trace() const {
        return p[0] + p[4] + p[8];
    }

    double det() const {
        double d = 0.0;
        //
        return d;
    }

    Mat3 inverse() const {
        Mat3 ret;
        //
        return ret;
    }

    Mat3 T() const {
        Mat3 ret;
        // for (int i=0; i<3; ++i)
        //     for (int j=0; j<3; ++j)
        //         // ret(j,i) = p[i*3+j];
        //         ret.p[j*3+i] = p[i*3+j];
        for (int i=0; i<8; i++) ret.p[(3*i)%8] = p[i];
        ret.p[8] = p[8];
        return ret;
    }

    std::array<double,9> p;
};

std::ostream& operator<< (std::ostream &out, Mat3 const& m) {
    out << '|' << m.p[0] << ' ' << m.p[1] << ' ' << m.p[2] << '|' << std::endl;
    out << '|' << m.p[3] << ' ' << m.p[4] << ' ' << m.p[5] << '|' << std::endl;
    out << '|' << m.p[6] << ' ' << m.p[7] << ' ' << m.p[8] << '|' << std::endl;
    return out;
}

Mat3 Eye(){
    Mat3 ret;
    ret.p[0] = 1.0;
    ret.p[4] = 1.0;
    ret.p[8] = 1.0;
    return ret;
}
