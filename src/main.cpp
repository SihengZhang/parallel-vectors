#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>

#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellTypes.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>

#include "interval.h"

struct CubicPoly {
    // Represents c3*x^3 + c2*x^2 + c1*x + c0
    double c0, c1, c2, c3;


    CubicPoly operator+(const CubicPoly& other) const {
        return {
            c0 + other.c0, 
            c1 + other.c1, 
            c2 + other.c2, 
            c3 + other.c3
        };
    }

    CubicPoly operator-(const CubicPoly& other) const {
        return {
            c0 - other.c0, 
            c1 - other.c1, 
            c2 - other.c2, 
            c3 - other.c3
        };
    }

    CubicPoly operator-() const {
        return {-c0, -c1, -c2, -c3};
    }
};

/**
 * Solves the cubic polynomial defined by the struct.
 * Returns a vector of real roots.
 */
std::vector<double> solveCubicReal(const CubicPoly& poly) {
    std::vector<double> roots;
    
    // Map struct members to algebraic names for readability
    // a*x^3 + b*x^2 + c*x + d = 0
    double a = poly.c3;
    double b = poly.c2;
    double c = poly.c1;
    double d = poly.c0;

    const double PI = 3.14159265358979323846;
    const double EPSILON = 1e-9;

    // 1. Handle non-cubic case (a = 0)
    if (std::abs(a) < EPSILON) {
        if (std::abs(b) < EPSILON) {
            // Linear: cx + d = 0
            if (std::abs(c) > EPSILON) roots.push_back(-d / c);
        } else {
            // Quadratic: bx^2 + cx + d = 0
            double quad_disc = c * c - 4 * b * d;
            if (quad_disc > EPSILON) {
                double sqrt_d = std::sqrt(quad_disc);
                roots.push_back((-c + sqrt_d) / (2 * b));
                roots.push_back((-c - sqrt_d) / (2 * b));
            } else if (std::abs(quad_disc) < EPSILON) {
                roots.push_back(-c / (2 * b));
            }
        }
        std::sort(roots.begin(), roots.end());
        return roots;
    }

    // 2. Normalize to x^3 + Ax^2 + Bx + C = 0
    double A = b / a;
    double B = c / a;
    double C = d / a;

    // 3. Depress to y^3 + py + q = 0
    double sqA = A * A;
    double p = B - (sqA / 3.0);
    double q = (2.0 * sqA * A / 27.0) - (A * B / 3.0) + C;

    // 4. Calculate Discriminant
    double discriminant = (q * q / 4.0) + (p * p * p / 27.0);
    double offset = A / 3.0;

    // --- CASE 1: 1 Real Root ---
    if (discriminant > EPSILON) {
        double sqrt_disc = std::sqrt(discriminant);
        double u_val = -q / 2.0 + sqrt_disc;
        double v_val = -q / 2.0 - sqrt_disc;

        double u = std::cbrt(u_val);
        double v = std::cbrt(v_val);
        
        roots.push_back(u + v - offset);
    } 
    // --- CASE 2: 3 Real Roots (Trigonometric Solution) ---
    else {
        double r_arg = -(p * p * p) / 27.0;
        double r = std::sqrt(std::max(0.0, r_arg)); 
        
        double phi_arg = std::max(-1.0, std::min(1.0, -q / (2.0 * r)));
        double phi = std::acos(phi_arg);
        
        double t = (r > 0.0) ? 2.0 * std::sqrt(-p / 3.0) : 0.0;

        roots.push_back(t * std::cos(phi / 3.0) - offset);
        roots.push_back(t * std::cos((phi + 2.0 * PI) / 3.0) - offset);
        roots.push_back(t * std::cos((phi + 4.0 * PI) / 3.0) - offset);
    }
    
    std::sort(roots.begin(), roots.end());
    return roots;
}

// compute det(A - lambda*B) expansion of 3x3 matrix
CubicPoly computeCharPoly(const std::vector<std::vector<double>>& A,
                          const std::vector<std::vector<double>>& B) {
    assert(A.size() == 3 && A[0].size() == 3);
    assert(B.size() == 3 && B[0].size() == 3);

    CubicPoly result;

    // Helper Lambda for 3x3 Determinant ---
    auto det3x3 = [](const std::vector<double>& r0, 
                     const std::vector<double>& r1, 
                     const std::vector<double>& r2) -> double {
        return r0[0] * (r1[1] * r2[2] - r1[2] * r2[1]) -
               r0[1] * (r1[0] * r2[2] - r1[2] * r2[0]) +
               r0[2] * (r1[0] * r2[1] - r1[1] * r2[0]);
    };

    // Calculate Coefficients for P(lambda) = c0 + c1*lam + c2*lam^2 + c3*lam^3
    // Constant term: det(A)
    result.c0 = det3x3(A[0], A[1], A[2]);

    // Linear term: - sum of dets with 1 row from B
    double term_1_1 = det3x3(B[0], A[1], A[2]);
    double term_1_2 = det3x3(A[0], B[1], A[2]);
    double term_1_3 = det3x3(A[0], A[1], B[2]);
    result.c1 = -(term_1_1 + term_1_2 + term_1_3);

    // Quadratic term: + sum of dets with 2 rows from B
    double term_2_1 = det3x3(A[0], B[1], B[2]);
    double term_2_2 = det3x3(B[0], A[1], B[2]);
    double term_2_3 = det3x3(B[0], B[1], A[2]);
    result.c2 = (term_2_1 + term_2_2 + term_2_3);

    // Cubic term: - det(B)
    result.c3 = -det3x3(B[0], B[1], B[2]);


    return result;                
}


int parallelVectorsInCell(const std::vector<std::vector<double>>& coordinates,
                           const std::vector<std::vector<double>>& vectorField_v,
                           const std::vector<std::vector<double>>& vectorField_w){

    assert(coordinates.size() == 4 && coordinates[0].size() == 3);
    assert(vectorField_v.size() == 4 && vectorField_v[0].size() == 3);
    assert(vectorField_w.size() == 4 && vectorField_w[0].size() == 3);

    // matrices for Q
    std::vector<std::vector<double>> matrix_A(3, std::vector<double>(3));
    std::vector<std::vector<double>> matrix_B(3, std::vector<double>(3));

    // matrices for Pi (i = 0, 1, 2, 3)
    std::vector<std::vector<double>> p0_A(3, std::vector<double>(3));
    std::vector<std::vector<double>> p0_B(3, std::vector<double>(3));
    std::vector<std::vector<double>> p1_A(3, std::vector<double>(3));
    std::vector<std::vector<double>> p1_B(3, std::vector<double>(3));
    std::vector<std::vector<double>> p2_A(3, std::vector<double>(3));
    std::vector<std::vector<double>> p2_B(3, std::vector<double>(3));
    std::vector<std::vector<double>> p3_A(3, std::vector<double>(3));
    std::vector<std::vector<double>> p3_B(3, std::vector<double>(3));

    std::vector<double> vector_a(3);
    std::vector<double> vector_b(3);

    
    // Fill matrices and vectors
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            matrix_A[i][j] = vectorField_v[i][j] - vectorField_v[3][j];
            matrix_B[i][j] = vectorField_w[i][j] - vectorField_w[3][j];
        }
    }
    for(int i = 0; i < 3; ++i) {p0_A[0][i] = vectorField_v[1][i]; p0_B[0][i] = vectorField_w[1][i];}
    for(int i = 0; i < 3; ++i) {p0_A[1][i] = vectorField_v[2][i]; p0_B[1][i] = vectorField_w[2][i];}
    for(int i = 0; i < 3; ++i) {p0_A[2][i] = vectorField_v[3][i]; p0_B[2][i] = vectorField_w[3][i];}
    for(int i = 0; i < 3; ++i) {p1_A[0][i] = vectorField_v[0][i]; p1_B[0][i] = vectorField_w[0][i];}
    for(int i = 0; i < 3; ++i) {p1_A[1][i] = vectorField_v[2][i]; p1_B[1][i] = vectorField_w[2][i];}
    for(int i = 0; i < 3; ++i) {p1_A[2][i] = vectorField_v[3][i]; p1_B[2][i] = vectorField_w[3][i];}
    for(int i = 0; i < 3; ++i) {p2_A[0][i] = vectorField_v[0][i]; p2_B[0][i] = vectorField_w[0][i];}
    for(int i = 0; i < 3; ++i) {p2_A[1][i] = vectorField_v[1][i]; p2_B[1][i] = vectorField_w[1][i];}
    for(int i = 0; i < 3; ++i) {p2_A[2][i] = vectorField_v[3][i]; p2_B[2][i] = vectorField_w[3][i];}
    for(int i = 0; i < 3; ++i) {p3_A[0][i] = vectorField_v[0][i]; p3_B[0][i] = vectorField_w[0][i];}
    for(int i = 0; i < 3; ++i) {p3_A[1][i] = vectorField_v[1][i]; p3_B[1][i] = vectorField_w[1][i];}
    for(int i = 0; i < 3; ++i) {p3_A[2][i] = vectorField_v[2][i]; p3_B[2][i] = vectorField_w[2][i];}

    vector_a = vectorField_v[3];
    vector_b = vectorField_w[3];

    // Calculate cubic polynomial
    CubicPoly Q = computeCharPoly(matrix_A, matrix_B);
    CubicPoly P0 = computeCharPoly(p0_A, p0_B);
    CubicPoly P1 = computeCharPoly(p1_A, p1_B);
    CubicPoly P2 = computeCharPoly(p2_A, p2_B);
    CubicPoly P3 = computeCharPoly(p3_A, p3_B);
    CubicPoly P0_ = Q - P0;
    CubicPoly P1_ = Q - P1;
    CubicPoly P2_ = Q - P2;
    CubicPoly P3_ = Q - P3;




    return EXIT_SUCCESS;



}


int main(int argc, char* argv[]) {

    // 1. Verify arguments
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename.vtu>" << std::endl;
        return EXIT_FAILURE;
    }
    std::string filename = argv[1];

    // 2. Create the reader
    auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    
    reader->SetFileName(filename.c_str());

    // 3. Update the reader to read the file
    try {
        reader->Update();
    } catch (...) {
        std::cerr << "Error reading file." << std::endl;
        return EXIT_FAILURE;
    }

    // 4. Get the output mesh
    vtkUnstructuredGrid* mesh = reader->GetOutput();

    // 5. Basic Verification
    std::cout << "Successfully read " << filename << std::endl;
    std::cout << "Number of Points: " << mesh->GetNumberOfPoints() << std::endl;
    std::cout << "Number of Cells:  " << mesh->GetNumberOfCells() << std::endl;
    

    vtkPointData* pd = mesh->GetPointData();
    std::cout << "Point Data Arrays: " << pd->GetNumberOfArrays() << std::endl;
    for(int i=0; i < pd->GetNumberOfArrays(); i++) {
        std::cout << " - " << pd->GetArrayName(i) << std::endl;
    }

    // Retrieve the pointers to your two arrays by name
    vtkDataArray* vectorArray1 = pd->GetArray("Velocity");
    vtkDataArray* vectorArray2 = pd->GetArray("Acceleration");

    if (!vectorArray1 || !vectorArray2) {
        std::cerr << "Error: Could not find required data arrays!" << std::endl;
        return EXIT_FAILURE;
    }

    vtkIdType numCells = mesh->GetNumberOfCells();

    for (vtkIdType cellId = 0; cellId < numCells; cellId++) {

        // Get the cell object
        vtkCell* cell = mesh->GetCell(cellId);

        // cell must be tet
        int cellType = cell->GetCellType();
        if(cellType != VTK_TETRA) {
            std::cout <<"The "<<cellId<<" cell type is: " << cellType << " (Not a Tet)" << std::endl;
            continue;
        }

        // Rows = Vertex Index (0-3), Cols = Component (x, y, z)
        std::vector<std::vector<double>> coordMatrix(4, std::vector<double>(3));
        std::vector<std::vector<double>> velocityMatrix(4, std::vector<double>(3));
        std::vector<std::vector<double>> accelerationMatrix(4, std::vector<double>(3));

        // Iterate over the points (vertices) of this specific tetrahedron
        vtkIdList* pointIds = cell->GetPointIds();

        for (vtkIdType i = 0; i < pointIds->GetNumberOfIds(); i++) {
            vtkIdType pointId = pointIds->GetId(i);

            // Get Geometry (Coordinates)
            double pt[3];
            mesh->GetPoint(pointId, pt);
        
            // Get vectors on this specific vertex
            double vVal[3]; 
            double aVal[3];
            vectorArray1->GetTuple(pointId, vVal);
            vectorArray2->GetTuple(pointId, aVal);


            for(int k=0; k<3; k++) {
                coordMatrix[i][k]        = pt[k];
                velocityMatrix[i][k]     = vVal[k];
                accelerationMatrix[i][k] = aVal[k];
            }
        }


        parallelVectorsInCell(coordMatrix, velocityMatrix, accelerationMatrix);

        if(cellId % 100000 == 0) {
            std::cout<<cellId<<std::endl;
        }
    

    }

    return EXIT_SUCCESS;
}