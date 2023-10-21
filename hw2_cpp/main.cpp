#include <igl/eigs.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/read_triangle_mesh.h>
#include <Eigen/Sparse>
#include <iostream>

Eigen::MatrixXd eigenvectors;
int selectedcolumn=0;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  std::cout<<"Key: "<<key<<" "<<(unsigned int)key<<std::endl;
  if (key == '1')
  {
    selectedcolumn = 1;
    viewer.data().set_data(eigenvectors.col(selectedcolumn));
  }
  else if (key == '2')
  {
    selectedcolumn = 2;
    viewer.data().set_data(eigenvectors.col(selectedcolumn));
  }
  else if (key == ' ')
  {
    selectedcolumn = (selectedcolumn + 1) % eigenvectors.cols();
    viewer.data().set_data(eigenvectors.col(selectedcolumn));
  }

  return false;
}



std::vector<std::set<int>> computeAdjacencies(Eigen::MatrixXi F, Eigen::MatrixXd V)
{
    std::vector<std::set<int>> adjacencies(V.rows());
    for (int i = 0; i < F.rows(); i++)
    {
        for (int j = 0; j < F.cols(); j++)
        {
            int v1 = F(i, j);
            int v2 = F(i, (j + 1) % F.cols());
            adjacencies[v1].insert(v2);
            adjacencies[v2].insert(v1);
        }
    }
    return adjacencies;
}

double calculateCotangentWeight(
    const Eigen::Vector3d& v1,
    const Eigen::Vector3d& v2,
    const Eigen::Vector3d& v3
) {
    // Calculate the edge vectors.
    Eigen::Vector3d e1 = v1 - v3;
    Eigen::Vector3d e2 = v2 - v3;

    // Calculate the cotangent weight using the dot product over the cross product formula.
    return e1.dot(e2) / (e1.cross(e2).norm());
}

// Define a function to compute the cotangent Laplacian operator using vertex adjacency.
Eigen::SparseMatrix<double> computeCotangentLaplacian(const Eigen::MatrixXd & vertices, const std::vector<std::set<int>>& vertexAdjacency) {
    int numVertices = vertices.rows();
    Eigen::SparseMatrix<double> CotangentLaplacian(numVertices, numVertices);

    // Create lists to store the triplet data for the sparse cotangent Laplacian matrix.
    std::vector<Eigen::Triplet<double>> tripletList;

    for (int i = 0; i < numVertices; ++i) {
        const std::set<int>& adjVertices = vertexAdjacency[i];

        // Compute cotangent weights for the current vertex.
        double diagonal = 0.0;

        for (int j : adjVertices) {
            // Calculate cotangent weight for the edge (i, j, k) where k is another adjacent vertex.
            for (int k : adjVertices) {
                if (k != i && k != j) {
                  //Check if k is adjacent to both 
                    if (vertexAdjacency[j].find(k) != vertexAdjacency[j].end()) {
                        // Calculate cotangent weight for the edge (i, j, k).
                        double cot_ijk = calculateCotangentWeight(vertices.row(i), vertices.row(j), vertices.row(k));
                        double w = cot_ijk /2.0;

                        // Update diagonal and off-diagonal entries.
                        diagonal += w;
                        tripletList.emplace_back(i, j, -w);
                    }
                }
            }
        }

        // Set the diagonal element in the Laplacian matrix.
        tripletList.emplace_back(i, i, diagonal);
    }

    // Create the sparse cotangent Laplacian matrix.
    CotangentLaplacian.setFromTriplets(tripletList.begin(), tripletList.end());

    return CotangentLaplacian;
}


int main(int argc, char * argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    if(!igl::read_triangle_mesh(
            argc>1?argv[1]: "../../input/bunny.obj",V,F))
    {
        std::cout<<"failed to load mesh"<<std::endl;
    }

    std::vector<std::set<int>> adjacencies = computeAdjacencies(F, V);

    Eigen::SparseMatrix<double> L = computeCotangentLaplacian(V, adjacencies);



    Eigen::SparseMatrix<double> L_eigen;

    igl::cotmatrix(V, F, L_eigen);

    // Find eigen vectors of L_eigen 
    // Eigen::MatrixXd V_eigen;
    // Eigen::VectorXd S_eigen;

    // // Calculate the mass matrix 
    // Eigen::SparseMatrix<double> M;
    // igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

    // eigs(L_eigen,M, 5, igl::EIGS_TYPE_SM, V_eigen, S_eigen);
    // std::cout << "Eigen values: " << std::endl;
    // std::cout << S_eigen << std::endl;
    // std::cout << "Eigen vectors: " << std::endl;
    // std::cout << V_eigen << std::endl;



    std::cout << "Laplacian matrix:" << std::endl;
    std::cout << L.row(90) << std::endl;
    std::cout << "Official Laplacian matrix for vertex "<< ":" << std::endl;
    std::cout << L_eigen.row(90) << std::endl;
    std::cout << "Is approx: ? " << L.isApprox(-L_eigen) << std::endl;



    





    // eigenvectors = V.block(0,0,V.rows(), 3);

    // std::cout << V << std::endl;
    // std::cout << F << std::endl;


    // igl::opengl::glfw::Viewer viewer;
    // igl::opengl::glfw::imgui::ImGuiPlugin plugin;
    // viewer.plugins.push_back(&plugin);
    // igl::opengl::glfw::imgui::ImGuiMenu menu;
    // plugin.widgets.push_back(&menu);
    
    // viewer.data().set_mesh(V,F);
    // viewer.callback_key_down = &key_down; // setting the callback
    // viewer.data().show_lines = true;
    // viewer.launch();
}