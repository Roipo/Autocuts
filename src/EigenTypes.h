#pragma once

#ifndef EIGENTYPES_H
#define EIGENTYPES_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

// Dense matrices
typedef Eigen::MatrixXd Mat;
typedef Eigen::Matrix4d Mat4;
typedef Eigen::Matrix3d Mat3;
typedef Eigen::Matrix2d Mat2;
typedef Eigen::MatrixX3d MatX3;
typedef Eigen::MatrixX2d MatX2;
typedef Eigen::Matrix3Xd Mat3X;
typedef Eigen::Matrix2Xd Mat2X;
typedef Eigen::MatrixX4d MatX4;
typedef Eigen::Matrix<double, 8, 2> Mat82;
typedef Eigen::Matrix<double, 6, 6> Mat6;

// special matrix used in Position.h/cpp
typedef Eigen::Matrix<double, 3, 2> Mat32;

// special matrix used in SolverPlugin.h/.cpp
typedef Eigen::Matrix<double, 4, 3> Mat43;

typedef Eigen::MatrixXi Mati;
typedef Eigen::Matrix4i Mat4i;
typedef Eigen::Matrix3i Mat3i;
typedef Eigen::Matrix2i Mat2i;
typedef Eigen::MatrixX3i MatX3i;
typedef Eigen::MatrixX2i MatX2i;
typedef Eigen::Matrix3Xi Mat3Xi;
typedef Eigen::Matrix2Xi Mat2Xi;
typedef Eigen::MatrixX4i MatX4i;
typedef Eigen::Matrix<int, 4, 3> Mat43i;

// Sparse matrices
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::SparseMatrix<int> SpMati;

// Column Vectors
typedef Eigen::VectorXd Vec;
typedef Eigen::Vector4d Vec4;
typedef Eigen::Vector3d Vec3;
typedef Eigen::Vector2d Vec2;
typedef Eigen::Matrix<double, 6, 1> Vec6;
typedef Eigen::Matrix<int, 6, 1> Vec6i;

typedef Eigen::VectorXi Veci;
typedef Eigen::Vector4i Vec4i;
typedef Eigen::Vector3i Vec3i;
typedef Eigen::Vector2i Vec2i;

// Row Vectors
typedef Eigen::RowVectorXd RVec;
typedef Eigen::RowVector4d RVec4;
typedef Eigen::RowVector3d RVec3;
typedef Eigen::RowVector2d RVec2;

typedef Eigen::RowVectorXi RVeci;
typedef Eigen::RowVector4i RVec4i;
typedef Eigen::RowVector3i RVec3i;
typedef Eigen::RowVector2i RVec2i;

// Some in floats
typedef Eigen::Vector2f Vec2f;

// Triplet vectors
typedef Eigen::Triplet<double> Tripletd;
typedef Eigen::Triplet<int> Tripleti;
typedef std::vector<Tripletd> Triplets;

// boolean stuff
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> Vecb;
typedef Eigen::Matrix<bool, 1, 2> RVec2b;

typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> Matb;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 2> MatX2b;

// unsigned char
typedef Eigen::Matrix<unsigned char, Eigen::Dynamic, 1> VecC;
typedef Eigen::Matrix<unsigned char, 1, Eigen::Dynamic> RVecC;
typedef Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> MatC;
typedef Eigen::Matrix<unsigned char, 2, 2> Mat2C;
typedef Eigen::Matrix<unsigned char, 3, 1> Vec3C;

// float stuff
typedef Eigen::Matrix<float, 3, 1> Vec3f;

#endif
