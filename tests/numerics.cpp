#include "../src/numerics.hpp"

namespace {
constexpr int GRID_SIZE = 8;
constexpr double TOLERANCE = 1e-12;
const Complex LOWER(-0.3, 0.2);
const Complex DIAGONAL(3.0, 0.4);
const Complex UPPER(0.1, -0.5);
const Complex FIELD_VALUE(0.6, 0.8);
const Complex INITIAL_VELOCITY(0.7, 0.0);
constexpr double MASS_RATIO = 0.125;
constexpr double GAMMA = 0.1;

void require(bool condition, const char* message)
{
    if (!condition) throw std::runtime_error(message);
}
}

int main()
{
    // Compare both Thomas solvers with an independent dense Eigen solve.
    for (bool periodic : {false, true}) {
        MatrixXcd matrix = MatrixXcd::Zero(GRID_SIZE, GRID_SIZE);
        matrix.diagonal().setConstant(DIAGONAL);
        matrix.diagonal(-1).setConstant(LOWER);
        matrix.diagonal(1).setConstant(UPPER);
        if (periodic) {
            matrix(0, GRID_SIZE - 1) = LOWER;
            matrix(GRID_SIZE - 1, 0) = UPPER;
        }
        VectorXcd rhs(GRID_SIZE);
        for (int index = 0; index < GRID_SIZE; ++index)
            rhs(index) = DIAGONAL + static_cast<double>(index) * LOWER;
        VectorXcd actual(GRID_SIZE);
        if (periodic) {
            Thomas_tridiagonal_solver_PBC(actual, LOWER, DIAGONAL, UPPER, rhs, GRID_SIZE);
        } else {
            VectorXcd diagonal = matrix.diagonal();
            Thomas_tridiagonal_solver(actual, LOWER, diagonal, UPPER, rhs, GRID_SIZE);
        }
        VectorXcd expected = matrix.partialPivLu().solve(rhs);
        require((actual - expected).norm() < TOLERANCE, "Thomas solver differs from dense solve");
    }

    VectorXcd field = VectorXcd::Constant(GRID_SIZE, FIELD_VALUE);
    VectorXcd density(GRID_SIZE), delta(GRID_SIZE);
    Non_Linear_Matrix(density, field, GRID_SIZE);
    require((density - VectorXcd::Ones(GRID_SIZE)).norm() < TOLERANCE,
            "Nonlinear density must equal squared field modulus");
    Delta_Function(delta, GRID_SIZE);
    VectorXcd expected_delta = VectorXcd::Zero(GRID_SIZE);
    expected_delta(GRID_SIZE / 2 - 1) = 1.0;
    require((delta - expected_delta).norm() < TOLERANCE, "Delta location changed");
    require(std::abs(Velocity_calculation(field, field, GRID_SIZE, INITIAL_VELOCITY,
                                        MASS_RATIO, GAMMA) - INITIAL_VELOCITY) < TOLERANCE,
            "Uniform field must preserve impurity velocity");
    std::cout << "Numerical function checks passed\n";
}
