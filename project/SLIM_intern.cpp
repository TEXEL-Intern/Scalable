/**
 * Scalable Locally Injective Mappings
*/

#include "SLIM_intern.h"



TrianglesMapping::TrianglesMapping(const int acount, char** avariable) {
    const char* name = nullptr;
    int weights;
    energy = nullptr;

    // Initialize the param-parser
    Parameters params;
    params.help = "This program processes a mesh file into a 2D map with specific energy and iteration settings.";
    params.add("input", "name", "project/mesh_test/hemisphere.obj").description("Input mesh file.");
    params.add("int", "weights", "1").description("Weights parameter.");
    params.add("int", "max_iterations", "50").description("Maximum number of iterations for the algorithm Generalized Reweighted local/global.");
    params.add("string", "energy", "SYMMETRIC-DIRICHLET").description("Energy type: ARAP, SYMMETRIC-DIRICHLET, EXPONENTIAL-SYMMETRIC-DIRICHLET, HENCKY-STRAIN, AMIPS, CONFORMAL-AMIPS-2D, UNTANGLE-2D.");
    params.add("float", "epsilon", "0.1").description("Epsilon value for the computation of UNTANGLE-2D.");
    params.init_from_args(acount, avariable);

    weights = std::stoi(params["weights"]);
    max_iterations = std::stoi(params["max_iterations"]);
    epsilon = std::stod(params["epsilon"]);

    if (acount > 1) {
        for (int i = 1; i < acount; ++i) {
            if (strlen(avariable[i]) > 1) {
                if (strncmp(avariable[i], "name=", 5) == 0) {
                    name = avariable[i] + 5;
                }
                if (strncmp(avariable[i], "energy=", 7) == 0) {
                    energy = avariable[i] + 7;
                }
            }
        }
    }

    if (name == nullptr) {
        #ifdef _WIN32
            name = "mesh_test/hemisphere.obj";
        #endif
        #ifdef linux
            name = "project/mesh_test/hemisphere.obj";
        #endif
    }

    if (energy == nullptr) {
        energy = "SYMMETRIC-DIRICHLET";
    }

    Tutte1963(name, weights);
}

const char* TrianglesMapping::getOutput() const {
    return output_name_obj;
}

double TrianglesMapping::unsignedArea(const vec3 &A, const vec3 &B, const vec3 &C) {
    return 0.5*cross(B-A, C-A).norm();
}

double TrianglesMapping::calculateCotan(const vec3& v0, const vec3& v1, const vec3& v2, const vec3& v3) {
    vec3 v = v0 - v1;
    vec3 w = v2 - v3;
    double cotan = v * w / cross(v, w).norm();
    return cotan;
}

void compute_surface_gradient_matrix(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &F1,
                                    const Eigen::MatrixXd &F2, Eigen::SparseMatrix<double> &D1, Eigen::SparseMatrix<double> &D2) {
    Eigen::SparseMatrix<double> G;
    igl::grad(V, F, G);
    Eigen::SparseMatrix<double> Dx = G.block(0, 0, F.rows(), V.rows());
    Eigen::SparseMatrix<double> Dy = G.block(F.rows(), 0, F.rows(), V.rows());
    Eigen::SparseMatrix<double> Dz = G.block(2 * F.rows(), 0, F.rows(), V.rows());

    D1 = F1.col(0).asDiagonal() * Dx + F1.col(1).asDiagonal() * Dy + F1.col(2).asDiagonal() * Dz;
    D2 = F2.col(0).asDiagonal() * Dx + F2.col(1).asDiagonal() * Dy + F2.col(2).asDiagonal() * Dz;
}

void TrianglesMapping::jacobian_rotation_area(Triangles& map, bool lineSearch) {
    if (!lineSearch) {xk_1 = Eigen::VectorXd::Zero(2 * num_vertices);}

    for (auto f : map.iter_facets()) {
        for (int j = 0; j < 3; ++j) {
            int v_ind = int(f.vertex(j));
            if (!lineSearch) {
                xk_1(v_ind) = f.vertex(j).pos()[0];
                xk_1(v_ind + num_vertices) = f.vertex(j).pos()[1];
            }
        }
    }

    Jac.clear();
    Rot.clear();
    if (!lineSearch) {Wei.clear();}

    if (strcmp(energy, "UNTANGLE-2D") == 0) {
        if (E_previous != -1) {
            double energy_sum = 0;

            Ji = Eigen::MatrixXd::Zero(num_triangles, dimension * dimension);
            Ji.col(0) = D_x * V_1.col(0);
            Ji.col(1) = D_y * V_1.col(0);
            Ji.col(2) = D_x * V_1.col(1);
            Ji.col(3) = D_y * V_1.col(1);

            Eigen::Matrix<double, 2, 2> ji, ri, ti, ui, vi;
            Eigen::Matrix<double, 2, 1> sing;
            ji(0, 0) = Ji(0, 0);
            ji(0, 1) = Ji(0, 1);
            ji(1, 0) = Ji(0, 2);
            ji(1, 1) = Ji(0, 3);
            igl::polar_svd(ji, ri, ti, ui, sing, vi);
            double detJ = ji.determinant();

            for (int i = 0; i < num_triangles; i++) {
                Eigen::Matrix<double, 2, 2> ji, ri, ti, ui, vi;
                Eigen::Matrix<double, 2, 1> sing;

                ji(0, 0) = Ji(i, 0);
                ji(0, 1) = Ji(i, 1);
                ji(1, 0) = Ji(i, 2);
                ji(1, 1) = Ji(i, 3);

                igl::polar_svd(ji, ri, ti, ui, sing, vi);

                double traceJTJ = (ji.transpose() * ji).trace();
                energy_sum += M(i) * 2 * ((traceJTJ / (detJ + sqrt(pow(epsilon, 2) + pow(detJ, 2)))) +
                                lambda_polyconvex * (pow(detJ, 2) + 1) / (detJ + sqrt(pow(epsilon, 2) + pow(detJ, 2))));
            }

            updateEpsilon(energy_sum);
        }
    }

    // Ji=[D1*u, D2*u, D1*v, D2*v];
    Ji = Eigen::MatrixXd::Zero(num_triangles, dimension * dimension);
    Ji.col(0) = D_x * V_1.col(0);
    Ji.col(1) = D_y * V_1.col(0);
    Ji.col(2) = D_x * V_1.col(1);
    Ji.col(3) = D_y * V_1.col(1);

    for (int i = 0; i < Ji.rows(); ++i) {
        Eigen::Matrix<double, 2, 2> ji, ri, ti, ui, vi;
        Eigen::Matrix<double, 2, 1> sing;
        Eigen::Matrix<double, 2, 1> closest_sing_vec;
        Eigen::Matrix<double, 2, 2> mat_W;
        Eigen::Matrix<double, 2, 1> m_sing_new;
        double s1, s2;

        ji(0, 0) = Ji(i, 0);
        ji(0, 1) = Ji(i, 1);
        ji(1, 0) = Ji(i, 2);
        ji(1, 1) = Ji(i, 3);

        igl::polar_svd(ji, ri, ti, ui, sing, vi);

        s1 = sing(0);
        s2 = sing(1);

        if (!lineSearch) {
            if (strcmp(energy, "ARAP") == 0) {
                m_sing_new << 1, 1;
            }
            else if (strcmp(energy, "SYMMETRIC-DIRICHLET") == 0) {
                double s1_g = 2 * (s1 - pow(s1, -3));
                double s2_g = 2 * (s2 - pow(s2, -3));

                m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
            }
            else if (strcmp(energy, "EXPONENTIAL-SYMMETRIC-DIRICHLET") == 0) {
                double s1_g = 2 * (s1 - pow(s1, -3));
                double s2_g = 2 * (s2 - pow(s2, -3));
                double inside_exponential = exponential_factor_1 * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2));
                double exponential_term = exp(inside_exponential);

                s1_g *= exponential_term * exponential_factor_1;
                s2_g *= exponential_term * exponential_factor_1;

                m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
            }
            else if (strcmp(energy, "HENCKY-STRAIN") == 0) {
                double s1_g = 2 * (log(s1) / s1);
                double s2_g = 2 * (log(s2) / s2);

                m_sing_new << sqrt(s1_g / (2 * (s1 - 1))), sqrt(s2_g / (2 * (s2 - 1)));
            }
            else if (strcmp(energy, "AMIPS") == 0) {
                double s1_g = 1;
                double s2_g = 1;
                double s1_lambda = sqrt((2 * pow(s2, 2) + 1) / (pow(s2, 2) + 2));
                double s2_lambda = sqrt((2 * pow(s1, 2) + 1) / (pow(s1, 2) + 2));
                double inside_exponential_1 = exponential_factor_2 * (0.5 * ((s1 / s2) + (s2 / s1)) + 0.25 * ((s1 * s2) + (1. / (s1 * s2))));
                double exponential_term_1 = exp(inside_exponential_1);
                double inside_exponential_2 = exponential_factor_2 * (0.5 * ((s1 / s2) + (s2 / s1)) + 0.25 * ((s1 * s2) + (1. / (s1 * s2))));
                double exponential_term_2 = exp(inside_exponential_2);

                s1_g *= exponential_term_1 * exponential_factor_2 * (0.5 * ((1. / s2) - (s2 / pow(s1, 2))) + 0.25 * (s2 - (1. / (s2 * pow(s1, 2)))));
                s2_g *= exponential_term_2 * exponential_factor_2 * (0.5 * ((1. / s1) - (s1 / pow(s2, 2))) + 0.25 * (s1 - (1. / (s1 * pow(s2, 2)))));

                m_sing_new << sqrt(s1_g / (2 * (s1 - s1_lambda))), sqrt(s2_g / (2 * (s2 - s2_lambda)));

                // Replacing the closest rotation R with another matrix Λ, which depends on the energy
                closest_sing_vec << s1_lambda, s2_lambda;
                ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
            }
            else if (strcmp(energy, "CONFORMAL-AMIPS-2D") == 0) {
                double s1_g = 1 / s2 - s2 / pow(s1, 2);
                double s2_g = 1 / s1 - s1 / pow(s2, 2);
                double geometric_average = sqrt(s1 * s2);
                double s1_lambda = geometric_average;
                double s2_lambda = geometric_average;

                m_sing_new << sqrt(s1_g / (2 * (s1 - s1_lambda))), sqrt(s2_g / (2 * (s2 - s2_lambda)));

                // Replacing the closest rotation R with another matrix Λ, which depends on the energy
                closest_sing_vec << s1_lambda, s2_lambda;
                ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();
            }
            else if (strcmp(energy, "UNTANGLE-2D") == 0) {
                double s1_g = 2 * (2 * s1 / (s1 * s2 + sqrt(pow(epsilon, 2) + pow(s1, 2) * pow(s2, 2))) -
                                ((pow(s1, 2) + pow(s2, 2)) * (s2 + (s1 * pow(s2, 2)) / sqrt(pow(epsilon, 2) + pow(s1, 2) * pow(s2, 2)))) /
                                pow((s1 * s2 + sqrt(pow(epsilon, 2) + pow(s1, 2) * pow(s2, 2))), 2) +
                                lambda_polyconvex * (2 * s1 * pow(s2, 2) / (s1 * s2 + sqrt(pow(epsilon, 2) + pow(s1, 2) * pow(s2, 2))) -
                                (pow(s1, 2) * pow(s2, 2) + 1) * (s2 + (s1 * pow(s2, 2)) / sqrt(pow(epsilon, 2) + pow(s1, 2) * pow(s2, 2)))) /
                                pow((s1 * s2 + sqrt(pow(epsilon, 2) + pow(s1, 2) * pow(s2, 2))), 2));

                double s2_g = 2 * (2 * s2 / (s1 * s2 + sqrt(pow(epsilon, 2) + pow(s1, 2) * pow(s2, 2))) -
                                ((pow(s2, 2) + pow(s1, 2)) * (s1 + (s2 * pow(s1, 2)) / sqrt(pow(epsilon, 2) + pow(s2, 2) * pow(s1, 2)))) /
                                pow((s1 * s2 + sqrt(pow(epsilon, 2) + pow(s1, 2) * pow(s2, 2))), 2) +
                                lambda_polyconvex * (2 * s2 * pow(s1, 2) / (s1 * s2 + sqrt(pow(epsilon, 2) + pow(s1, 2) * pow(s2, 2))) -
                                (pow(s2, 2) * pow(s1, 2) + 1) * (s1 + (s2 * pow(s1, 2)) / sqrt(pow(epsilon, 2) + pow(s2, 2) * pow(s1, 2)))) /
                                pow((s1 * s2 + sqrt(pow(epsilon, 2) + pow(s1, 2) * pow(s2, 2))), 2));

                double solution1 = sqrt((2 * (1 + lambda_polyconvex * pow(s2, 2)) * sqrt(pow(pow(s2, 2) + lambda_polyconvex, 2) * pow(s2, 4) - 
                                    (1 + lambda_polyconvex * pow(s2, 2)) * (pow(s2, 2) + lambda_polyconvex) * pow(s2, 2) * pow(epsilon, 2) + 
                                    pow(1 + lambda_polyconvex * pow(s2, 2), 2) * pow(epsilon, 4)) - 
                                    (1 + lambda_polyconvex * pow(s2, 2)) * (2 * pow(epsilon, 2) * (1 + lambda_polyconvex * pow(s2, 2)) - 
                                    pow(s2, 2) * (pow(s2, 2) + lambda_polyconvex))) / 
                                    (3 * pow(s2, 2) * pow(1 + lambda_polyconvex * pow(s2, 2), 2)));
                
                double solution2 = sqrt((-2 * (1 + lambda_polyconvex * pow(s2, 2)) * sqrt(pow(pow(s2, 2) + lambda_polyconvex, 2) * pow(s2, 4) - 
                                    (1 + lambda_polyconvex * pow(s2, 2)) * (pow(s2, 2) + lambda_polyconvex) * pow(s2, 2) * pow(epsilon, 2) + 
                                    pow(1 + lambda_polyconvex * pow(s2, 2), 2) * pow(epsilon, 4)) - 
                                    (1 + lambda_polyconvex * pow(s2, 2)) * (2 * pow(epsilon, 2) * (1 + lambda_polyconvex * pow(s2, 2)) - 
                                    pow(s2, 2) * (pow(s2, 2) + lambda_polyconvex))) / 
                                    (3 * pow(s2, 2) * pow(1 + lambda_polyconvex * pow(s2, 2), 2)));
                
                double solution3 = -sqrt((2 * (1 + lambda_polyconvex * pow(s2, 2)) * sqrt(pow(pow(s2, 2) + lambda_polyconvex, 2) * pow(s2, 4) - 
                                    (1 + lambda_polyconvex * pow(s2, 2)) * (pow(s2, 2) + lambda_polyconvex) * pow(s2, 2) * pow(epsilon, 2) + 
                                    pow(1 + lambda_polyconvex * pow(s2, 2), 2) * pow(epsilon, 4)) - 
                                    (1 + lambda_polyconvex * pow(s2, 2)) * (2 * pow(epsilon, 2) * (1 + lambda_polyconvex * pow(s2, 2)) - 
                                    pow(s2, 2) * (pow(s2, 2) + lambda_polyconvex))) / 
                                    (3 * pow(s2, 2) * pow(1 + lambda_polyconvex * pow(s2, 2), 2)));
                
                double solution4 = -sqrt((-2 * (1 + lambda_polyconvex * pow(s2, 2)) * sqrt(pow(pow(s2, 2) + lambda_polyconvex, 2) * pow(s2, 4) - 
                                    (1 + lambda_polyconvex * pow(s2, 2)) * (pow(s2, 2) + lambda_polyconvex) * pow(s2, 2) * pow(epsilon, 2) + 
                                    pow(1 + lambda_polyconvex * pow(s2, 2), 2) * pow(epsilon, 4)) - 
                                    (1 + lambda_polyconvex * pow(s2, 2)) * (2 * pow(epsilon, 2) * (1 + lambda_polyconvex * pow(s2, 2)) - 
                                    pow(s2, 2) * (pow(s2, 2) + lambda_polyconvex))) / 
                                    (3 * pow(s2, 2) * pow(1 + lambda_polyconvex * pow(s2, 2), 2)));
                
                double s1_lambda;
                s1_lambda = solution1;
                
                solution1 = sqrt((2 * (1 + lambda_polyconvex * pow(s1, 2)) * sqrt(pow(pow(s1, 2) + lambda_polyconvex, 2) * pow(s1, 4) - 
                                    (1 + lambda_polyconvex * pow(s1, 2)) * (pow(s1, 2) + lambda_polyconvex) * pow(s1, 2) * pow(epsilon, 2) + 
                                    pow(1 + lambda_polyconvex * pow(s1, 2), 2) * pow(epsilon, 4)) - 
                                    (1 + lambda_polyconvex * pow(s1, 2)) * (2 * pow(epsilon, 2) * (1 + lambda_polyconvex * pow(s1, 2)) - 
                                    pow(s1, 2) * (pow(s1, 2) + lambda_polyconvex))) / 
                                    (3 * pow(s1, 2) * pow(1 + lambda_polyconvex * pow(s1, 2), 2)));
                
                solution2 = sqrt((-2 * (1 + lambda_polyconvex * pow(s1, 2)) * sqrt(pow(pow(s1, 2) + lambda_polyconvex, 2) * pow(s1, 4) - 
                                    (1 + lambda_polyconvex * pow(s1, 2)) * (pow(s1, 2) + lambda_polyconvex) * pow(s1, 2) * pow(epsilon, 2) + 
                                    pow(1 + lambda_polyconvex * pow(s1, 2), 2) * pow(epsilon, 4)) - 
                                    (1 + lambda_polyconvex * pow(s1, 2)) * (2 * pow(epsilon, 2) * (1 + lambda_polyconvex * pow(s1, 2)) - 
                                    pow(s1, 2) * (pow(s1, 2) + lambda_polyconvex))) / 
                                    (3 * pow(s1, 2) * pow(1 + lambda_polyconvex * pow(s1, 2), 2)));
                
                solution3 = -sqrt((2 * (1 + lambda_polyconvex * pow(s1, 2)) * sqrt(pow(pow(s1, 2) + lambda_polyconvex, 2) * pow(s1, 4) - 
                                    (1 + lambda_polyconvex * pow(s1, 2)) * (pow(s1, 2) + lambda_polyconvex) * pow(s1, 2) * pow(epsilon, 2) + 
                                    pow(1 + lambda_polyconvex * pow(s1, 2), 2) * pow(epsilon, 4)) - 
                                    (1 + lambda_polyconvex * pow(s1, 2)) * (2 * pow(epsilon, 2) * (1 + lambda_polyconvex * pow(s1, 2)) - 
                                    pow(s1, 2) * (pow(s1, 2) + lambda_polyconvex))) / 
                                    (3 * pow(s1, 2) * pow(1 + lambda_polyconvex * pow(s1, 2), 2)));
                
                solution4 = -sqrt((-2 * (1 + lambda_polyconvex * pow(s1, 2)) * sqrt(pow(pow(s1, 2) + lambda_polyconvex, 2) * pow(s1, 4) - 
                                    (1 + lambda_polyconvex * pow(s1, 2)) * (pow(s1, 2) + lambda_polyconvex) * pow(s1, 2) * pow(epsilon, 2) + 
                                    pow(1 + lambda_polyconvex * pow(s1, 2), 2) * pow(epsilon, 4)) - 
                                    (1 + lambda_polyconvex * pow(s1, 2)) * (2 * pow(epsilon, 2) * (1 + lambda_polyconvex * pow(s1, 2)) - 
                                    pow(s1, 2) * (pow(s1, 2) + lambda_polyconvex))) / 
                                    (3 * pow(s1, 2) * pow(1 + lambda_polyconvex * pow(s1, 2), 2)));
                
                double s2_lambda;
                s2_lambda = solution1;

                m_sing_new << sqrt(s1_g / (2 * (s1 - s1_lambda))), sqrt(s2_g / (2 * (s2 - s2_lambda)));

                // Replacing the closest rotation R with another matrix Λ, which depends on the energy
                closest_sing_vec << s1_lambda, s2_lambda;
                ri = ui * closest_sing_vec.asDiagonal() * vi.transpose();

                if (!std::isnan(s1_g / (2 * (s1 - s1_lambda)))) {
                    m_sing_new(0) = 1;
                }
                if (!std::isnan(s2_g / (2 * (s2 - s2_lambda)))) {
                    m_sing_new(1) = 1;
                }
            }

            Jac.push_back(ji);
            Rot.push_back(ri);
            if (std::abs(s1 - 1) < 1e-8) m_sing_new(0) = 1;
            if (std::abs(s2 - 1) < 1e-8) m_sing_new(1) = 1;
            mat_W = ui * m_sing_new.asDiagonal() * ui.transpose();
            Wei.push_back(mat_W);
        }
    }
}

void TrianglesMapping::least_squares() {
    Eigen::SparseMatrix<double> A(dimension * dimension * num_triangles, dimension * num_vertices);
    Eigen::VectorXd b(dimension * dimension * num_triangles);
    b.setZero();
    
	// Create diagonal matrices W11, W12, W21, W22 as vectors
	Eigen::VectorXd W11_diag = Eigen::VectorXd::Zero(num_triangles);
	Eigen::VectorXd W12_diag = Eigen::VectorXd::Zero(num_triangles);
	Eigen::VectorXd W21_diag = Eigen::VectorXd::Zero(num_triangles);
	Eigen::VectorXd W22_diag = Eigen::VectorXd::Zero(num_triangles);

	// R vectors
	Eigen::VectorXd R11 = Eigen::VectorXd::Zero(num_triangles);
	Eigen::VectorXd R12 = Eigen::VectorXd::Zero(num_triangles);
	Eigen::VectorXd R21 = Eigen::VectorXd::Zero(num_triangles);
	Eigen::VectorXd R22 = Eigen::VectorXd::Zero(num_triangles);

	for (int i = 0; i < num_triangles; ++i) {
		W11_diag(i) = Wei[i](0, 0);
		W12_diag(i) = Wei[i](0, 1);
		W21_diag(i) = Wei[i](1, 0);
		W22_diag(i) = Wei[i](1, 1);

		R11(i) = Rot[i](0, 0);
		R12(i) = Rot[i](0, 1);
		R21(i) = Rot[i](1, 0);
		R22(i) = Rot[i](1, 1);
	}

    std::vector<Eigen::Triplet<double>> triplet;
    triplet.reserve(4 * (D_x.outerSize() + D_y.outerSize()));
    // Fill in the A matrix
    /*A = [W11*Dx, W12*Dx;
           W11*Dy, W12*Dy;
           W21*Dx, W22*Dx;
           W21*Dy, W22*Dy];*/
    for (int k = 0; k < D_x.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(D_x, k); it; ++it) {
        int dx_row = it.row();
        int dx_col = it.col();
        double val = it.value();

        triplet.push_back(Eigen::Triplet<double>(dx_row, dx_col, val * W11_diag(dx_row)));
        triplet.push_back(Eigen::Triplet<double>(dx_row, num_vertices + dx_col, val * W12_diag(dx_row)));

        triplet.push_back(Eigen::Triplet<double>(2 * num_triangles + dx_row, dx_col, val * W21_diag(dx_row)));
        triplet.push_back(Eigen::Triplet<double>(2 * num_triangles + dx_row, num_vertices + dx_col, val * W22_diag(dx_row)));
        }
    }
    
    for (int k = 0; k < D_y.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(D_y, k); it; ++it) {
        int dy_row = it.row();
        int dy_col = it.col();
        double val = it.value();

        triplet.push_back(Eigen::Triplet<double>(num_triangles + dy_row, dy_col, val * W11_diag(dy_row)));
        triplet.push_back(Eigen::Triplet<double>(num_triangles + dy_row, num_vertices + dy_col, val * W12_diag(dy_row)));

        triplet.push_back(Eigen::Triplet<double>(3 * num_triangles + dy_row, dy_col, val * W21_diag(dy_row)));
        triplet.push_back(Eigen::Triplet<double>(3 * num_triangles + dy_row, num_vertices + dy_col, val * W22_diag(dy_row)));
        }
    }

	A.setFromTriplets(triplet.begin(), triplet.end());
    Eigen::SparseMatrix<double> At = A.transpose();
    At.makeCompressed();

    Eigen::SparseMatrix<double> identity_dimA(At.rows(), At.rows());
    identity_dimA.setIdentity();

    // Add a proximal term
    flattened_weight_matrix.resize(dimension * dimension * num_triangles);
    mass.resize(num_triangles);
	mass.setConstant(weight_option); // All the weights are equal to 1
    for (int i = 0; i < dimension * dimension; i++)
        for (int j = 0; j < num_triangles; j++)
            flattened_weight_matrix(i * num_triangles + j) = mass(j);
    Eigen::SparseMatrix<double> L;
    double lambda = 1e-4;
    Eigen::SparseMatrix<double> uv_dimA(At.rows(), At.rows());
    L = At * flattened_weight_matrix.asDiagonal() * A + lambda * identity_dimA;
    L.makeCompressed();

	// Fill in the b vector
    /*b = [W11*R11 + W12*R21;
           W11*R12 + W12*R22;
           W21*R11 + W22*R21;
           W21*R12 + W22*R22];*/
    for (int i = 0; i < num_triangles; i++) {
        b(i + 0 * num_triangles) = W11_diag(i) * R11(i) + W12_diag(i) * R21(i);
        b(i + 1 * num_triangles) = W11_diag(i) * R12(i) + W12_diag(i) * R22(i);
        b(i + 2 * num_triangles) = W21_diag(i) * R11(i) + W22_diag(i) * R21(i);
        b(i + 3 * num_triangles) = W21_diag(i) * R12(i) + W22_diag(i) * R22(i);
    }

    Eigen::VectorXd uv_flat_1(dimension * num_vertices);
    for (int j = 0; j < num_vertices; j++) {
        uv_flat_1(0 * num_vertices + j) = V_1(j, 0);
        uv_flat_1(1 * num_vertices + j) = V_1(j, 1);
    }
    rhs.resize(dimension * num_vertices);
    rhs = (At * flattened_weight_matrix.asDiagonal() * b + lambda * uv_flat_1);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    xk = solver.compute(L).solve(rhs);
    if(solver.info() != Eigen::Success) {
		// Solving failed
		std::cerr << "Solving failed" << std::endl;
		return;
	}

    uv = V_1.block(0, 0, V_1.rows(), 2);

    for (int i = 0; i < dimension; i++) {
        uv.col(i) = xk.block(i * num_vertices, 0, num_vertices, 1);
    }

    distance = uv - V_1.block(0, 0, V_1.rows(), 2);
    dk = xk - xk_1;
}

double TrianglesMapping::step_singularities(const Eigen::MatrixXi& F, const Eigen::MatrixXd& uv, const Eigen::MatrixXd& d) {
    double maximum_step = INFINITY;
    for (int f = 0; f < F.rows(); f++) {
        // Get quadratic coefficients (ax^2 + bx + c)

        int v1 = F(f, 0); int v2 = F(f, 1); int v3 = F(f, 2);

        const double& U11 = uv(v1, 0);
        const double& U12 = uv(v1, 1);
        const double& U21 = uv(v2, 0);
        const double& U22 = uv(v2, 1);
        const double& U31 = uv(v3, 0);
        const double& U32 = uv(v3, 1);

        const double& V11 = d(v1, 0);
        const double& V12 = d(v1, 1);
        const double& V21 = d(v2, 0);
        const double& V22 = d(v2, 1);
        const double& V31 = d(v3, 0);
        const double& V32 = d(v3, 1);
        
        
        double a = V11*V22 - V12*V21 - V11*V32 + V12*V31 + V21*V32 - V22*V31;
        double b = U11*V22 - U12*V21 - U21*V12 + U22*V11 - U11*V32 + U12*V31 + U31*V12 - U32*V11 + U21*V32 - U22*V31 - U31*V22 + U32*V21;
        double c = U11*U22 - U12*U21 - U11*U32 + U12*U31 + U21*U32 - U22*U31;
        
        double minimum_positive_root = smallest_position_quadratic_zero(a, b, c);
        maximum_step = std::min(maximum_step, minimum_positive_root);
    }
    return maximum_step;
}

double TrianglesMapping::smallest_position_quadratic_zero(double a, double b, double c) {
    double x1, x2;
    if (a != 0) {
        double delta = pow(b, 2) - 4*a*c;
        if (delta < 0) {
            return INFINITY;
        }
        delta = sqrt(delta);
        x1 = (-b + delta)/ (2*a);
        x2 = (-b - delta)/ (2*a);
    } else {
        x1 = x2 = -b/c;
    }
    assert(std::isfinite(x1));
    assert(std::isfinite(x2));

    double temp = std::min(x1, x2);
    x1 = std::max(x1, x2); x2 = temp;
    if (x1 == x2) {
        return INFINITY; // Means the orientation flips twice, so does it flip?
    }
    // Return the smallest negative root if it exists, otherwise return infinity
    if (x1 > 0) {
        if (x2 > 0) {
            return x2;
        } else {
            return x1;
        }
    } else {
        return INFINITY;
    }
}

double TrianglesMapping::add_energies_jacobians(Eigen::MatrixXd& V_new, bool flips_linesearch) {
	double energy_sum = 0;
    distortion_energy.clear();
    number_inverted = 0;

    Ji = Eigen::MatrixXd::Zero(num_triangles, dimension * dimension);
    Ji.col(0) = D_x * V_new.col(0);
    Ji.col(1) = D_y * V_new.col(0);
    Ji.col(2) = D_x * V_new.col(1);
    Ji.col(3) = D_y * V_new.col(1);

    Eigen::Matrix<double, 2, 2> ji, ri, ti, ui, vi;
    Eigen::Matrix<double, 2, 1> sing;
    ji(0, 0) = Ji(0, 0);
    ji(0, 1) = Ji(0, 1);
    ji(1, 0) = Ji(0, 2);
    ji(1, 1) = Ji(0, 3);
    igl::polar_svd(ji, ri, ti, ui, sing, vi);
    double s1 = sing(0);
    double s2 = sing(1);
    minimum_determinant = s1 * s2;
    maximum_divide_singular = s1 / s2;

	for (int i = 0; i < num_triangles; i++) {
        Eigen::Matrix<double, 2, 2> ji, ri, ti, ui, vi;
        Eigen::Matrix<double, 2, 1> sing;
        double mini_energy = 0;

        ji(0, 0) = Ji(i, 0);
        ji(0, 1) = Ji(i, 1);
        ji(1, 0) = Ji(i, 2);
        ji(1, 1) = Ji(i, 3);

        igl::polar_svd(ji, ri, ti, ui, sing, vi);
        double s1 = sing(0);
        double s2 = sing(1);

        double detJ = ji.determinant();
        double divide = s1 / s2;
        if (detJ <= 0) {
            number_inverted++;
        }
        minimum_determinant = std::min(minimum_determinant, detJ);
        maximum_divide_singular = std::max(maximum_divide_singular, divide);

		if (flips_linesearch) {
            if (strcmp(energy, "ARAP") == 0) {
                mini_energy = pow(s1 - 1, 2) + pow(s2 - 1, 2);
                energy_sum += M(i) * mini_energy;
            }
            else if (strcmp(energy, "SYMMETRIC-DIRICHLET") == 0) {
                mini_energy = pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2);
                energy_sum += M(i) * mini_energy;
            }
            else if (strcmp(energy, "EXPONENTIAL-SYMMETRIC-DIRICHLET") == 0) {
                mini_energy = exp(exponential_factor_1 * (pow(s1, 2) + pow(s1, -2) + pow(s2, 2) + pow(s2, -2)));
                energy_sum += M(i) * mini_energy;
            }
            else if (strcmp(energy, "HENCKY-STRAIN") == 0) {
                mini_energy = pow(log(s1), 2) + pow(log(s2), 2);
                energy_sum += M(i) * mini_energy;
            }
            else if (strcmp(energy, "AMIPS") == 0) {
                mini_energy = exp(exponential_factor_2 * (0.5 * ((s1 / s2) + (s2 / s1)) + 0.25 * ((s1 * s2) + (1. / (s1 * s2)))));
                energy_sum += M(i) * mini_energy;
            }
            else if (strcmp(energy, "CONFORMAL-AMIPS-2D") == 0) {
                mini_energy = (pow(s1, 2) + pow(s2, 2)) / (s1 * s2);
                energy_sum += M(i) * mini_energy;
            }
            else if (strcmp(energy, "UNTANGLE-2D") == 0) {
                double traceJTJ = (ji.transpose() * ji).trace();
                mini_energy = 2 * ((traceJTJ / (detJ + sqrt(pow(epsilon, 2) + pow(detJ, 2)))) +
                                lambda_polyconvex * (pow(detJ, 2) + 1) / (detJ + sqrt(pow(epsilon, 2) + pow(detJ, 2))));
                energy_sum += M(i) * mini_energy;
            }
		} else {
			if (ui.determinant() * vi.determinant() > 0) {
			energy_sum += M(i) * (pow(s1-1,2) + pow(s2-1,2));
			} else {
			vi.col(1) *= -1;
			energy_sum += M(i) * (Jac[i]-ui*vi.transpose()).squaredNorm();
			}
		}
        distortion_energy.push_back(mini_energy);
	}
    return energy_sum;
}

double TrianglesMapping::lineSearch(Eigen::MatrixXd& xk_current, Eigen::MatrixXd& dk_current) {
    // Line search using the bisection method
    double alphaMax;
    if (strcmp(energy, "UNTANGLE-2D") != 0) {
        alphaMax = step_singularities(F, xk_current, dk_current);
    } else {
        alphaMax = 1.25;
    }
    double alphaBisectionMethod = std::min(1.0, 0.8 * alphaMax);

    Eigen::MatrixXd V_old = xk_current.block(0, 0, xk_current.rows(), 2);
    double ener, new_ener;
    double current_energy;
    ener = add_energies_jacobians(V_old, true);
    new_ener = ener;

    int max_iter = 12;
    int iter = 0;

    if (strcmp(energy, "UNTANGLE-2D") != 0) {
        while (new_ener >= ener && iter < max_iter) {
            Eigen::MatrixXd V_new = V_old + alphaBisectionMethod * dk_current;
            current_energy = add_energies_jacobians(V_new, true);

            if (current_energy >= ener) {
                alphaBisectionMethod /= 2.0;
            } else {
                V_old = V_new;
                new_ener = current_energy;
            }

            iter++;
        }
    } else {
        while (new_ener >= ener && iter < max_iter) {
            Eigen::MatrixXd V_new = V_old + alphaBisectionMethod * dk_current;
            current_energy = add_energies_jacobians(V_new, true);

            if (current_energy >= ener) {
                alphaBisectionMethod /= 2.0;
            } else {
                V_old = V_new;
                new_ener = current_energy;
            }

            iter++;
        }
    }
    
    energumene = new_ener;
    V_1 = V_old;
    return alphaBisectionMethod;
}

void TrianglesMapping::updateEpsilon(double instant_E) {
    double E = instant_E;

    double sigma = std::max(1. - E / E_previous, 1e-1);
    if (minimum_determinant >= 0) {
        if (minimum_determinant >= 0) {
            epsilon *= (1 - sigma);
        } else {
            epsilon *= 1 - (sigma * std::sqrt(minimum_determinant * minimum_determinant + epsilon * epsilon)) / (std::abs(minimum_determinant) + std::sqrt(minimum_determinant * minimum_determinant + epsilon * epsilon));
        }
    } else {
        double mu = (1 - sigma) * ((minimum_determinant + std::sqrt(epsilon * epsilon + minimum_determinant * minimum_determinant)) / 2);
        if (minimum_determinant < mu) {
            epsilon = std::max(1e-9, 2 * std::sqrt(mu * (mu - minimum_determinant)));
        } else {
            epsilon = 1e-9;
        }
    }
}

void TrianglesMapping::nextStep(Triangles& map) {
	// Perform line search to find step size alpha
	alpha = lineSearch(V_1, distance);
    std::cout << "Energy: " << energumene << " | Inverted triangles: " << number_inverted << " | Minimum determinant: " << minimum_determinant << " | Epsilon: " << epsilon << " | Alpha: " << alpha << std::endl;
	
    for (int i = 0; i < V_1.rows(); i++) {
        map.points[i][0] = V_1(i, 0);
        map.points[i][1] = V_1(i, 1);
    }

    E_previous = energumene;
    number_inverted_previous = number_inverted;
}

void TrianglesMapping::swapVertices(Eigen::MatrixXd& V_1, Triangles& map, int point1, int point2) {
    auto x = V_1(point1, 0);
    auto y = V_1(point1, 1);

    V_1(point1, 0) = V_1(point2, 0);
    V_1(point1, 1) = V_1(point2, 1);
    V_1(point2, 0) = x;
    V_1(point2, 1) = y;

    map.points[point1][0] = V_1(point1, 0);
    map.points[point1][1] = V_1(point1, 1);
    map.points[point2][0] = V_1(point2, 0);
    map.points[point2][1] = V_1(point2, 1);
}

void TrianglesMapping::map_vertices_to_circle_area_normalized(Triangles& map, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi& bnd, Eigen::MatrixXd& UV) {
    Eigen::VectorXd dblArea_orig;
    igl::doublearea(V, F, dblArea_orig);
    double area = dblArea_orig.sum() / 2;
    double radius = sqrt(area / (M_PI)); 

    // Get sorted list of boundary vertices
    std::vector<int> interior, map_ij;
    map_ij.resize(V.rows());
    interior.reserve(V.rows() - bnd.size());

    std::vector<bool> isOnBnd(V.rows(), false);
    for (int i = 0; i < bnd.size(); i++) {
        isOnBnd[bnd[i]] = true;
        map_ij[bnd[i]] = i;
    }

    for (int i = 0; i < (int)isOnBnd.size(); i++) {
        if (!isOnBnd[i]) {
                map_ij[i] = interior.size();
                interior.push_back(i);
            }
    }

    std::vector<double> len(bnd.size());
    len[0] = 0.;

    for (int i = 1; i < bnd.size(); i++) {
        len[i] = len[i-1] + (V.row(bnd[i-1]) - V.row(bnd[i])).norm();
    }
    double total_len = len[len.size()-1] + (V.row(bnd[0]) - V.row(bnd[bnd.size()-1])).norm();

    UV.resize(bnd.size(), 2);
    for (int i = 0; i < bnd.size(); i++) {
        double frac = len[i] * (2. * M_PI) / total_len;
        UV.row(map_ij[bnd[i]]) << radius * cos(frac), radius * sin(frac);
        map.points[bnd[i]][0] = radius * cos(frac);
        map.points[bnd[i]][1] = radius * sin(frac);
    }
}

// Initialization: Tutte embedding
void TrianglesMapping::Tutte1963(const char* name, int weights) {
    std::filesystem::path filepath = name;
    std::string filepath_str_ext = filepath.extension().string();
    std::string filepath_str_stem = filepath.stem().string();
    const char* stem = filepath_str_stem.c_str();

    char ext2[12] = ".geogram";
    char method[20] = "_barycentre";
    char weight1[20] = "_uniform";
    char weight2[20] = "_cotan";
    char weight3[20] = "_random";
    char weight4[20] = "_sanity_check";
    char weight5[25] = "_texture_coordinates";
    bool texture_coordinates = false;

    // Create directory if it does not exist
    std::filesystem::path stem_dir(stem);
    if (!std::filesystem::exists(stem_dir)) {
        std::filesystem::create_directory(stem_dir);
    }

    std::filesystem::path energy_dir = stem_dir / energy;
    if (!std::filesystem::exists(energy_dir)) {
        std::filesystem::create_directory(energy_dir);
    }

    strcpy(output_name_geo, stem);
    strcat(output_name_geo, "/");
    strcat(output_name_geo, energy);
    strcat(output_name_geo, "/");
    strcat(output_name_geo, stem);
    strcat(output_name_geo, method);
    if (weights == 1) {
        strcat(output_name_geo, weight1);
    } else if (weights == 2) {
        strcat(output_name_geo, weight2);
    } else if (weights == 3) {
        strcat(output_name_geo, weight3);
    } else if (weights == 4) {
        strcat(output_name_geo, weight4);
        sanity_check = true;
        weights = 1;
    } else if (weights == 5) {
        strcat(output_name_geo, weight5);
        texture_coordinates = true;
        weights = 1;
    }
    strcpy(output_name_obj, output_name_geo);
    strcat(output_name_geo, ext2);
    strcat(output_name_obj, ".obj");

    strcpy(times_txt, stem);
    strcat(times_txt, "/");
    strcat(times_txt, energy);
    strcat(times_txt, "/");
    strcat(times_txt, stem);
    strcat(times_txt, ".txt");
    std::ofstream timeFile(times_txt); // Open a file for writing times
    auto start = std::chrono::high_resolution_clock::now();
    totalTime = 0; // Initialize total time accumulator
    
    igl::read_triangle_mesh(name, V, F);
    SurfaceAttributes attr = read_by_extension(name, mTut);
    PointAttribute<vec2> tex_coord("tex_coord", attr, mTut);
    mTut.connect();
    read_by_extension(name, mOri);
    mOri.connect();
    igl::doublearea(V, F, M);
    M /= 2.;
    mesh_area = M.sum();
    V_1 = V;

    int nverts = mOri.nverts();
    int nfacets = mOri.nfacets();
    int ncorners = mOri.ncorners();
    double seaLevel = 0.;

    std::cout << "The number of vertices is: " << nverts << ", facets: " << nfacets << ", corners: " << ncorners << std::endl;

    std::unordered_map<int, std::vector<int>> neighbor;
    std::unordered_set<int> bound;
    std::set<int> bound_halfedges;

    for (auto he : mOri.iter_halfedges()) {
        if (!he.opposite().active()) {
            neighbor[he.from()].push_back(he.to());
            neighbor[he.to()].push_back(he.from());
            bound.insert(he.from());
            bound.insert(he.to());

            bound_halfedges.insert(he);
        }
    }
    
    int missing = 0;
    if (!bound.empty()) {
        std::unordered_set<int> visited;
        int start = *bound.begin();
        bound_sorted.push_back(start);
        visited.insert(start);

        while (bound_sorted.size() < bound.size()) {
            int last = bound_sorted.back();
            bool found = false;
            
            for (int next : neighbor[last]) {
                if (visited.find(next) == visited.end()) {
                    bound_sorted.push_back(next);
                    visited.insert(next);
                    found = true;
                    break;
                }
            }

            if (!found) {
                missing++;
                break;
            }
        }
    }

    int fixed = 0;
    Eigen::VectorXd x_B_ = Eigen::VectorXd::Zero(nverts);
    for (int i = 0; i < mOri.nverts(); i++) {
        Surface::Vertex vi = Surface::Vertex(mOri, i);
        if (bound.contains(vi)) {
            blade.insert(vi);
            x_B_(i) = fixed;
            fixed++;
        }
    }

    Eigen::VectorXd x_I_ = Eigen::VectorXd::Zero(nverts);
    int insider = 0;
    std::set<int> plane;
    for (int i = 0; i < mOri.nverts(); i++) {
        Surface::Vertex vi = Surface::Vertex(mOri, i);
        if (!bound.contains(vi)) {
            plane.insert(vi);
            x_I_(i) = insider;
            insider++;
        }
    }

    if (!texture_coordinates) {
        Eigen::SparseMatrix<double> A_II(nverts - fixed, nverts - fixed);
        Eigen::SparseMatrix<double> A_IB(nverts - fixed, fixed);

        Eigen::SparseMatrix<double> A_II_A_BB(nverts, nverts);
        for (int i = 0; i < fixed; ++i) {
            A_II_A_BB.insert(i, i) = 1;
        }

        Eigen::MatrixXd lhsF = Eigen::MatrixXd::Zero(nverts, 2);

        Eigen::VectorXi b;
        igl::boundary_loop(F, b);
        Eigen::MatrixXd bc;
        map_vertices_to_circle_area_normalized(mTut, V, F, b, bc);

        int vert = 0;
        int bb = 0;
        for (int i = 0; i < mTut.nverts(); i++) {
            Surface::Vertex vi = Surface::Vertex(mOri, i);
            if (bound.contains(vi)) {
                lhsF(bb, 0) = mTut.points[i][0];
                lhsF(bb, 1) = mTut.points[i][1];
                bb++;
            } else {
                Surface::Halfedge depart = Surface::Vertex(mOri, i).halfedge();
                Surface::Halfedge variable = depart;
                double count = 0;
                std::vector<int> neighbors;
                if (weights == 1) {
                    std::map<int, double> cotan;
                    if (depart.opposite().active())
                        variable = variable.opposite().next();
                    neighbors.push_back(depart.to());
                    cotan.insert(std::make_pair(neighbors.back(), 1));
                    count += 1;
                    while (depart != variable && variable.active()) {
                        neighbors.push_back(variable.to());
                        cotan.insert(std::make_pair(neighbors.back(), 1));
                        count += 1;
                        if (!variable.opposite().active())
                            break;
                        variable = variable.opposite().next();
                    }

                    int ree = x_I_(i);
                    A_II.insert(ree, ree) = -count;
                    A_II_A_BB.insert(ree + fixed, ree + fixed) = -count;

                    for (auto const& [key, val] : cotan) {
                        if (blade.contains(key)) {
                            int re_ne2 = x_B_(key);
                            A_IB.insert(ree, re_ne2) = val;
                            A_II_A_BB.insert(ree + fixed, re_ne2) = val;
                        } else {
                            int re_ne = x_I_(key);
                            A_II.insert(ree, re_ne) = val;
                            A_II_A_BB.insert(ree + fixed, re_ne + fixed) = val;
                        }
                    }
                } else if (weights == 2) {
                    std::map<int, double> cotan;
                    if (depart.opposite().active())
                        variable = variable.opposite().next();

                    double cotan_alpha = calculateCotan(depart.next().from().pos(), depart.next().to().pos(), depart.prev().to().pos(), depart.prev().from().pos());
                    double cotan_beta = calculateCotan(depart.opposite().prev().to().pos(), depart.opposite().prev().from().pos(), depart.opposite().next().from().pos(), depart.opposite().next().to().pos());
                    double cotan_gamma = calculateCotan(depart.next().to().pos(), depart.next().from().pos(), depart.from().pos(), depart.to().pos());

                    double w_ij = cotan_alpha + cotan_beta;
                    double voronoi = 0.125 * (cotan_gamma * (depart.from().pos() - depart.to().pos()).norm2() + cotan_alpha * (depart.prev().to().pos() - depart.prev().from().pos()).norm2());

                    neighbors.push_back(depart.to());
                    cotan.insert(std::make_pair(neighbors.back(), w_ij));
                    count += w_ij;

                    while (depart != variable && variable.active()) {
                        cotan_alpha = calculateCotan(variable.next().from().pos(), variable.next().to().pos(), variable.prev().to().pos(), variable.prev().from().pos());
                        cotan_beta = calculateCotan(variable.opposite().prev().to().pos(), variable.opposite().prev().from().pos(), variable.opposite().next().from().pos(), variable.opposite().next().to().pos());
                        cotan_gamma = calculateCotan(variable.next().to().pos(), variable.next().from().pos(), variable.from().pos(), variable.to().pos());

                        w_ij = cotan_alpha + cotan_beta;
                        voronoi += 0.125 * (cotan_gamma * (variable.from().pos() - variable.to().pos()).norm2() + cotan_alpha * (variable.prev().to().pos() - variable.prev().from().pos()).norm2());

                        neighbors.push_back(variable.to());
                        cotan.insert(std::make_pair(neighbors.back(), w_ij));
                        count += w_ij;
                        if (!variable.opposite().active())
                            break;
                        variable = variable.opposite().next();
                    }

                    int ree = x_I_(i);
                    A_II.insert(ree, ree) = -count / (2 * voronoi);
                    A_II_A_BB.insert(ree + fixed, ree + fixed) = -count / (2 * voronoi);

                    for (auto const& [key, val] : cotan) {
                        if (blade.contains(key)) {
                            int re_ne2 = x_B_(key);
                            A_IB.insert(ree, re_ne2) = val / (2 * voronoi);
                            A_II_A_BB.insert(ree + fixed, re_ne2) = val / (2 * voronoi);
                        } else {
                            int re_ne = x_I_(key);
                            A_II.insert(ree, re_ne) = val / (2 * voronoi);
                            A_II_A_BB.insert(ree + fixed, re_ne + fixed) = val / (2 * voronoi);
                        }
                    }
                } else if (weights == 3) {
                    std::map<int, double> cotan;
                    std::random_device rd; // Obtain a random number from hardware
                    std::mt19937 gen(rd()); // Seed the generator
                    std::uniform_real_distribution<double> distr(-0.05, -0.01);
                    std::uniform_real_distribution<double> decision(0, 1);
                    double random_number;
                    double decision_number = decision(gen);
                    if (decision_number < 0.8) {
                        random_number = 1;
                    } else {
                        random_number = distr(gen);
                    }
                    if (depart.opposite().active())
                        variable = variable.opposite().next();
                    neighbors.push_back(depart.to());
                    cotan.insert(std::make_pair(neighbors.back(), random_number));
                    count += random_number;
                    while (depart != variable && variable.active()) {
                        double random_number_G;
                        double decision_number_G = decision(gen);
                        if (decision_number_G < 0.8) {
                            random_number_G = 1;
                        } else {
                            random_number_G = distr(gen);
                        }
                        neighbors.push_back(variable.to());
                        cotan.insert(std::make_pair(neighbors.back(), random_number_G));
                        count += random_number_G;
                        if (!variable.opposite().active())
                            break;
                        variable = variable.opposite().next();
                    }

                    int ree = x_I_(i);
                    A_II.insert(ree, ree) = -count;
                    A_II_A_BB.insert(ree + fixed, ree + fixed) = -count;

                    for (auto const& [key, val] : cotan) {
                        if (blade.contains(key)) {
                            int re_ne2 = x_B_(key);
                            A_IB.insert(ree, re_ne2) = val;
                            A_II_A_BB.insert(ree + fixed, re_ne2) = val;
                        } else {
                            int re_ne = x_I_(key);
                            A_II.insert(ree, re_ne) = val;
                            A_II_A_BB.insert(ree + fixed, re_ne + fixed) = val;
                        }
                    }
                }
            }
            vert++;
        }

        A_II_A_BB.makeCompressed();
        Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
        solver.analyzePattern(A_II_A_BB);
        solver.factorize(A_II_A_BB);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Decomposition failed" << std::endl;
            return;
        }
        
        Eigen::MatrixXd x_I_full = solver.solve(lhsF);
        if (solver.info() != Eigen::Success) {
            std::cerr << "Solving failed" << std::endl;
            return;
        }

        EigenMap = x_I_full;

        for (int plan : plane) {
            int re = x_I_(plan);
            mTut.points[plan][0] = x_I_full(re + fixed, 0);
            mTut.points[plan][1] = x_I_full(re + fixed, 1);
            mTut.points[plan][2] = seaLevel;
            V_1(plan, 0) = x_I_full(re + fixed, 0);
            V_1(plan, 1) = x_I_full(re + fixed, 1);
            V_1(plan, 2) = seaLevel;
        }

        for (int blad : blade) {
            mTut.points[blad][2] = seaLevel;
            V_1(blad, 0) = mTut.points[blad][0];
            V_1(blad, 1) = mTut.points[blad][1];
            V_1(blad, 2) = seaLevel;
        }

        if (sanity_check) {
            int n = 1; // Set n to the desired number of pairs to swap (n=1 for 2 elements, n=2 for 4 elements, etc.)
            std::set<int>::iterator it = plane.begin();
            for (int i = 0; i < n && it != plane.end(); ++i) {
                int point1 = *it;
                it++;
                if (it == plane.end()) break;
                int point2 = *it;
                it++;

                swapVertices(V_1, mTut, point1, point2);
            }
        }
    } else {
        for (int i = 0; i < mTut.nverts(); i++) {
            for (int d : range(2)) {
                mTut.points[i][d] = tex_coord[i][d];
                V_1(i, d) = mTut.points[i][d];
            }
            mTut.points[i][2] = seaLevel;
            V_1(i, 2) = seaLevel;
        }
    }

    Eigen::MatrixXd F1, F2, F3;
    igl::local_basis(V, F, F1, F2, F3);
    compute_surface_gradient_matrix(V, F, F1, F2, D_x, D_y);
    D_x.makeCompressed();
    D_y.makeCompressed();
    D_z.makeCompressed();
    
    num_vertices = nverts;
    num_triangles = nfacets;

    Eigen::MatrixXd uv_old = V_1.block(0, 0, V_1.rows(), 2);
    energumene = add_energies_jacobians(uv_old, true);
    std::cout << "Energy: " << energumene << " | Inverted triangles: " << number_inverted << " | Minimum determinant: " << minimum_determinant << std::endl;

    FacetAttribute<double> fa2(mOri);
    for (auto f : mOri.iter_facets()) {
        double area = unsignedArea(f.vertex(0).pos(), f.vertex(1).pos(), f.vertex(2).pos());
        fa2[f] = area;
        fOriMap[int(f)] = area;
    }

    CornerAttribute<double> he(mTut);
    for (auto f : mTut.iter_halfedges()) {
        if (blade.contains(f.from()) || blade.contains(f.to())) {
            he[f] = 404;
        } else {
            he[f] = 0;
        }
    }

    Surface::Facet f(mTut, 0);
    double minArea = unsignedArea(f.vertex(0).pos(), f.vertex(1).pos(), f.vertex(2).pos()) / fOriMap[0];
    double maxArea = minArea;
    FacetAttribute<double> fa_a(mTut);
    for (auto f : mTut.iter_facets()) {
        fa_a[f] = unsignedArea(f.vertex(0).pos(), f.vertex(1).pos(), f.vertex(2).pos()) / fa2[f];
        if (fa_a[f] < minArea) {
            minArea = fa_a[f];
        } else if (fa_a[f] > maxArea) {
            maxArea = fa_a[f];
        }
    }

    double minEnergy = distortion_energy[0];
    double maxEnergy = distortion_energy[0];
    FacetAttribute<double> fa(mTut);
    for (auto f : mTut.iter_facets()) {
        fa[f] = distortion_energy[int(f)];
        if (fa[f] < minEnergy) {
            minEnergy = fa[f];
        } else if (fa[f] > maxEnergy) {
            maxEnergy = fa[f];
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    totalTime += duration; // Accumulate total time
    if (timeFile.is_open()) {
        timeFile << "0" << "|" << duration << "|"; // Write iteration number and duration to file
    }

    if (timeFile.is_open()) {
        timeFile << totalTime << "|"; // Log total time
        timeFile << minArea << "|" << maxArea << "|" << minEnergy << "|" << maxEnergy << "|";
        timeFile << mLocGlo.nverts() << "|" << mLocGlo.nfacets() << "|" << mLocGlo.ncorners() << "|" << alpha << "|" << energumene << "|" << minimum_determinant << "|" << maximum_divide_singular;
        if (strcmp(energy, "UNTANGLE-2D") == 0) {
            timeFile << "|" << epsilon << "|" << lambda_polyconvex << "|" << number_inverted << "\n";
        }
        else {
            timeFile << "\n";
        }
    }

    if (timeFile.is_open()) {
        timeFile.close();
    }

    write_by_extension(output_name_geo, mTut, { {}, {{"Energy", fa.ptr}, {"AreaRatio", fa_a.ptr}}, {{"Halfedge", he.ptr}} });
    std::string output_name_str(output_name_obj);
    size_t pos = output_name_str.rfind(".obj");
    if (pos != std::string::npos) {
        output_name_str = output_name_str.substr(0, pos);
    }
    std::string modified_output_name_obj = output_name_str + "_with_MAP.obj";
    for (int i = 0; i < mTut.nverts(); i++) {
        for (int d : range(2)) {
            tex_coord[i][d] = mTut.points[i][d];
        }
    }
    write_by_extension(modified_output_name_obj.c_str(), mOri, attr);
    write_by_extension(output_name_obj, mTut);
}

void TrianglesMapping::LocalGlobalParametrization(const char* map) {
    read_by_extension(map, mLocGlo);
    mLocGlo.connect();

    std::filesystem::path filepath = map;
    std::string filepath_str_ext = filepath.extension().string();
    std::string filepath_str_stem = filepath.stem().string();
    const char* stem = filepath_str_stem.c_str();
    const char* first_space = strchr(stem, '_'); // Find the first space in stem
    size_t first_word_length = first_space ? (size_t)(first_space - stem) : strlen(stem); // Calculate length of the first word

    std::ofstream timeFile(times_txt, std::ios::app); // Append mode

    char ext2[12] = ".geogram";
    char method[20] = "_local_global_";
    char numStr[20];
    for (int i = 1; i <= max_iterations; ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        jacobian_rotation_area(mLocGlo, false);
        least_squares();
        nextStep(mLocGlo);
        std::cout << "FIN ITERATION " << i << std::endl;
        
        output_name_geo[0] = '\0'; // Clear output_name
        strncpy(output_name_geo, stem, first_word_length);
        output_name_geo[first_word_length] = '\0'; // Ensure null-termination
        strcat(output_name_geo, "/");
        strcat(output_name_geo, energy);
        strcat(output_name_geo, "/");
        strncat(output_name_geo, stem, first_word_length);

        strcat(output_name_geo, method);
        strcat(output_name_geo, energy);
        strcat(output_name_geo, "_");
        sprintf(numStr, "%d", i);
        strcat(output_name_geo, numStr);
        strcpy(output_name_obj, output_name_geo);
        strcat(output_name_geo, ext2);
        strcat(output_name_obj, ".obj");

        CornerAttribute<double> he(mLocGlo);
        for (auto f : mLocGlo.iter_halfedges()) {
            if (blade.contains(f.from()) || blade.contains(f.to())) {
                he[f] = 404;
            } else {
                he[f] = 0;
            }
        }

        Surface::Facet f(mLocGlo, 0);
        double minArea = unsignedArea(f.vertex(0).pos(), f.vertex(1).pos(), f.vertex(2).pos()) / fOriMap[0];
        double maxArea = minArea;
        FacetAttribute<double> fa_a(mLocGlo);
        for (auto f : mLocGlo.iter_facets()) {
            fa_a[f] = unsignedArea(f.vertex(0).pos(), f.vertex(1).pos(), f.vertex(2).pos()) / fOriMap[int(f)];
            if (fa_a[f] < minArea) {
                minArea = fa_a[f];
            } else if (fa_a[f] > maxArea) {
                maxArea = fa_a[f];
            }
        }

        double minEnergy = distortion_energy[0];
        double maxEnergy = distortion_energy[0];
        FacetAttribute<double> fa(mLocGlo);
        for (auto f : mLocGlo.iter_facets()) {
                fa[f] = distortion_energy[int(f)];
                if (fa[f] < minEnergy) {
                    minEnergy = fa[f];
                } else if (fa[f] > maxEnergy) {
                    maxEnergy = fa[f];
                }
        }
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        totalTime += duration;
        if (timeFile.is_open()) {
            timeFile << i << "|" << duration << "|";
        }

        if (timeFile.is_open()) {
            timeFile << totalTime << "|";
            timeFile << minArea << "|" << maxArea << "|" << minEnergy << "|" << maxEnergy << "|";
            timeFile << mLocGlo.nverts() << "|" << mLocGlo.nfacets() << "|" << mLocGlo.ncorners() << "|" << alpha << "|" << energumene << "|" << minimum_determinant << "|" << maximum_divide_singular;
            if (strcmp(energy, "UNTANGLE-2D") == 0) {
                timeFile << "|" << epsilon << "|" << lambda_polyconvex << "|" << number_inverted << "\n";
            }
            else {
                timeFile << "\n";
            }
        }

        write_by_extension(output_name_geo, mLocGlo, { {}, {{"Energy", fa.ptr}, {"AreaRatio", fa_a.ptr}}, {{"Halfedge", he.ptr}} });
        write_by_extension(output_name_obj, mLocGlo);
	}

    if (timeFile.is_open()) {
        timeFile.close();
    }
}

int main(int argc, char** argv) {

    std::cout << "Initialization..." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    TrianglesMapping Init(argc, argv);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Time taken: " << duration << " milliseconds" << std::endl;
    
	Init.LocalGlobalParametrization(Init.getOutput());

    return 0;
}