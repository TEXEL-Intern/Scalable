#ifndef SLIM_INTERN_H
#define SLIM_INTERN_H


#include "helpers.h"
#include <ultimaille/all.h>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <filesystem>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <thread>
#include <chrono>
#include <cmath>
#include <eigen3/Eigen/Sparse>
#include <fstream>
#include <random>
#include <igl/grad.h>
#include <igl/local_basis.h>
#include <igl/read_triangle_mesh.h>
#include <igl/polar_svd.h>
#include <igl/boundary_loop.h>
#include <param_parser/param_parser.h>



using namespace UM;



class TrianglesMapping {
public:
    TrianglesMapping(const int acount, char** avariable);

    const char* getOutput() const;
    void LocalGlobalParametrization(const char* map);

    Eigen::VectorXd xk;
    Eigen::VectorXd xk_1;
    Eigen::VectorXd dk;
    Eigen::MatrixXd uv;
    Eigen::MatrixXd distance;

private:
    Triangles mOri;
    Triangles mTut;
    Triangles mLocGlo;
    std::set<int> blade;
    std::set<int> bound;
    std::vector<int> bound_sorted;
    std::unordered_map<int, double> fOriMap;
    std::vector<double> distortion_energy;
    Eigen::MatrixXd EigenMap;
    char output_name_geo[250];
    char output_name_obj[250];
    char times_txt[150];
    int num_vertices;
	int num_triangles;
    Eigen::MatrixXd V; // card(V) by 3, list of mesh vertex positions
    Eigen::MatrixXi F; // card(F) by 3/3, list of mesh faces (triangles/tetrahedra)
    Eigen::MatrixXd V_1;
    Eigen::MatrixXi F_1;
    Eigen::MatrixXd Ji;
    std::vector<Eigen::Matrix2d> Rot, Jac, Wei;
    Eigen::SparseMatrix<double> D_x, D_y, D_z;
    Eigen::VectorXd flattened_weight_matrix;
    Eigen::VectorXd mass;
    double weight_option = 1.0;
    double exponential_factor_1 = 1e-3;
    double exponential_factor_2 = 1e-1;
    Eigen::VectorXd rhs;
    double alpha;
    Eigen::VectorXd M;
    double mesh_area;
    double energumene;
    double E_previous = -1;
    double lambda_polyconvex = 1;
    double epsilon = 1e-1;
    double minimum_determinant;
    double maximum_divide_singular;
    int number_inverted;
    int number_inverted_previous = 1e3;
    bool sanity_check = false;
    bool texture_coordinates = false;
    int dimension = 2;
    const char* energy;
    int max_iterations;
    long long totalTime;

    double unsignedArea(const vec3 &A, const vec3 &B, const vec3 &C);
    double calculateCotan(const vec3& v0, const vec3& v1, const vec3& v2, const vec3& v3);
    void swapVertices(Eigen::MatrixXd& V_1, Triangles& map, int point1, int point2);
    void map_vertices_to_circle_area_normalized(Triangles& map, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd, Eigen::MatrixXd& UV);
    void Tutte1963(const char* name, int weights);
    void jacobian_rotation_area(Triangles& map, bool lineSearch);
    void least_squares();
    double step_singularities(const Eigen::MatrixXi& F, const Eigen::MatrixXd& uv, const Eigen::MatrixXd& d);
    double smallest_position_quadratic_zero(double a, double b, double c);
    double add_energies_jacobians(Eigen::MatrixXd& V_new, bool flips_linesearch);
    double lineSearch(Eigen::MatrixXd& xk_current, Eigen::MatrixXd& dk);
    void updateEpsilon(double instant_E);
    void nextStep(Triangles& map);
};


#endif // #ifndef SLIM_INTERN_H