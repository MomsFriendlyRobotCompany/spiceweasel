// MIT License
//
// Copyright (c) 2022 Mom's Friendly Robot Company
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
#pragma once

#if defined(__linux__) || defined(__APPLE__)
  #include <Eigen/Dense>
  #include <functional>
  using cemd = const Eigen::MatrixXd;
  using cevd = const Eigen::VectorXd;
  using emd  = Eigen::MatrixXd;
  using evd  = Eigen::VectorXd;

  #include "ekf.hpp"
  #include "kf.hpp"
#endif

#include "kf_1d.hpp"
#include "median_filter.hpp"

// #include <cstring>
// #include <cstddef>
// #include <stdexcept>
// #include <ostream>
// #include <iostream>
// #include <algorithm>
// #include <stdexcept>

// using std::cout;
// using std::endl;
// using std::size_t;
// using std::memcpy;

// #define CHECK_ARRAY 1

// class Array {
//   public:
//   Array(size_t r, size_t c): r(r), c(c), size(r*c), buffer(new float[r*c]) {}
//   Array(Array&& a) noexcept : r(0), c(0), size(0), buffer(nullptr) {
//     *this = std::move(a);
//     cout << "move constructor "<< r << " " << c << endl;
//   }
//   Array(const Array& a): r(a.r), c(a.c), size(a.size), buffer(new float[r*c])
//   {
//     cout << "copy constructor"<<endl;
//     memcpy(buffer, a.data(), size*sizeof(buffer[0]));
//   }

//   ~Array() { delete[] buffer; }

//   friend std::ostream &operator<<(std::ostream &s, const Array &a);

//   Array& operator=(Array&& a) {
//     cout << "move constructor ="<<endl;
//     // check(r,c,a.r,a.c,__LINE__);
//     if (this != &a) {
//       delete[] buffer;
//       // r = a.r;
//       // c = a.c;
//       // size = a.size;
//       // buffer = a.buffer;
//       // a.buffer = nullptr;
//       // a.r = 0;
//       // a.c = 0;
//       // a.size = 0;
//       *this = std::move(a);
//     }
//     return *this;
//   }
//   Array& operator=(const Array& a) {
//     cout << "copy constructor ="<<endl;
//     // check(r,c,a.r,a.c,__LINE__);
//     memcpy(buffer, a.data(), size*sizeof(buffer[0]));
//     return *this;
//   }

//   float& operator()(size_t rr, size_t cc) { return buffer[rr*c+cc]; }
//   float operator()(size_t rr, size_t cc) const { return buffer[rr*c+cc]; }
//   float& operator[](size_t i) { return buffer[i]; }
//   float operator[](size_t i) const { return buffer[i]; }

//   float* data() const { return buffer; }

//   Array operator+(const Array& a) const {
//     check(r,c,a.r,a.c,__PRETTY_FUNCTION__,__LINE__);
//     Array ret = a;
//     for (int i=0; i < size; ++i) ret[i] += buffer[i];
//     return ret;
//   }
//   Array operator-(const Array& a) const {
//     // check(r,c,a.r,a.c,__LINE__);
//     Array ret = a;
//     for (int i=0; i < size; ++i) ret[i] -= buffer[i];
//     return ret;
//   }
//   // Array operator*(const Array& a) const {
//   //   check(r,c,a.r,a.c,__LINE__);
//   //   Array ret;
//   //   for (int i = 0; i < c; i++){ // cols
//   //     for (int j = 0; j < r; j++){ // rows
//   //         ret(i, j) = buffer[j*3+0]*m.p[i*3+j];
//   //     }
//   //   }
//   //   return ret;
//   // }

//   const size_t r,c,size;

//   // protected:
//   float *buffer=nullptr;

//   protected:
//   inline void check(size_t ar, size_t ac, size_t br, size_t bc,
//                     std::string fun, int line) const {
//     #if CHECK_ARRAY
//     if (ar != br || ac != bc) throw std::runtime_error(
//     std::string("\033[0;31m")
//                           + "\n\ndimensions don't match: \n"
//                           + std::string(__FILE__) + ":" +
//                           std::to_string(line) + "\n"
//                           + fun + "\033[0m\n");
//     #endif
//   }
// };

// std::ostream &operator<<(std::ostream &s, const Array &a) {
//   s << "(" << a.r << "," << a.c << ") [";
//   for (int i=0; i < a.r*a.c; ++i) s << a[i] << ",";
//   s << "]";
//   return s;
// }

// namespace np {

// template<typename T, size_t R, size_t C>
// class Array {
//   public:
//   Array(): size(R*C), r(R), c(C) { memset(&buffer, 0, size); }
//   Array(const Array& a): size(a.r*a.c), r(a.r), c(a.c)  {
//     memcpy(&buffer, &a.buffer, a.size);
//   }

//   friend std::ostream &operator<<(std::ostream &s, const Array &a);
//   //   s << "(" << a.r << "," << a.c << ") [";
//   //   for (int i=0; i < a.size; ++i) s << a[i] << ",";
//   //   s << "]";
//   //   return s;
//   // }
//   // std::ostream &operator<<(std::ostream &s) {
//   //   s << "(" << r << "," << c << ") [";
//   //   for (int i=0; i < size; ++i) s << ",";
//   //   s << "]";
//   //   return s;
//   // }

//   void operator=(const Array& a) {
//     check(a.r, a.c);
//     memcpy(&buffer, &a.buffer, a.size);
//   }
//   T operator()(size_t r, size_t c) { return buffer[r*c+c]; }
//   T operator[](size_t i) { return buffer[i]; }
//   const size_t size;
//   const size_t r;
//   const size_t c;

//   protected:
//   T buffer[R*C];
//   inline void check(size_t row, size_t col) {
//     if (row != r) throw std::runtime_error("rows don't match");
//     if (col != c) throw std::runtime_error("cols don't match");
//   }
// };

// template<typename T, size_t R, size_t C>
// Array<T,R, C> zeros() { return Array<T,R,C>(); }

// template<typename T, size_t R, size_t C>
// std::ostream &operator<<(std::ostream &s, const Array<T,R,C> &a) {
//   s << "(" << a.r << "," << a.c << ") [";
//   for (int i=0; i < a.size; ++i) s << a[i] << ",";
//   s << "]";
//   return s;
// }

// };
