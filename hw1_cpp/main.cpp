#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include "trimesh.h"
#include <cmath>
/* 


Each iteration of subdivision scheme usually consists of the following steps:
• insert new vertices
• modify mesh connectivity (remove some of the old connections, add
new ones)
• relocate “old” (existing before this iteration) vertices, e.g., set each to
weighted average of its neighbors.
• relocate “new” vertices.

*/

void subdivide( trimesh::trimesh_t& mesh, Eigen::MatrixXd& oldV, Eigen::MatrixXi& oldF, std::vector< trimesh::edge_t >& edges , Eigen::MatrixXd& newV, Eigen::MatrixXi& newF)
{

    // 1. insert new vertices
    // 2. modify mesh connectivity (remove some of the old connections, add new ones)
    // 3. relocate “old” (existing before this iteration) vertices, e.g., set each to weighted average of its neighbors.
    // 4. relocate “new” vertices.
    
    // Eigen::MatrixXi oldF;
    // oldF = mesh.get_faces();
    // Eigen::MatrixXi NewF;
    // std::cout<< "\n" << oldF << "\n";

    int num_edges = edges.size();
    // Eigen::MatrixXd newV(oldV.rows() + num_edges, oldV.cols());
    // Eigen::MatrixXi newF(oldF.rows()*4, oldF.cols());
    std::vector<int> old_to_new(oldV.rows(), -1);
    std::map<std::pair<int,int>, int> edge_to_new;
    int indexcount = 0;
    int scheme = 2;
    for (int i=0; i<oldF.rows(); ++i){
        Eigen::Vector3i face = oldF.row(i);
        // std::cout<< "\n" << face << "\n";

        Eigen::Vector3d pt_1 = oldV.row(face(0));
        Eigen::Vector3d pt_3 = oldV.row(face(1));
        Eigen::Vector3d pt_5 = oldV.row(face(2));

        Eigen::Vector3d pt_2 = (pt_1 + pt_3)/2;
        Eigen::Vector3d pt_4 = (pt_3 + pt_5)/2;
        Eigen::Vector3d pt_6 = (pt_5 + pt_1)/2;

        if (old_to_new[face(0)] == -1){
            old_to_new[face(0)] = indexcount;
            indexcount++;
        }
        newV.row(old_to_new[face(0)]) = pt_1.transpose();
        //check if the edge is already in the map
        if (edge_to_new.find(std::make_pair(std::min(face(0), face(1)),std::max(face(0),face(1)))) == edge_to_new.end()){
            edge_to_new[std::make_pair(std::min(face(0), face(1)),std::max(face(0),face(1)))] = indexcount;
            indexcount++;
        }
        newV.row(edge_to_new[std::make_pair(std::min(face(0), face(1)),std::max(face(0),face(1)))]) = pt_2.transpose();
        
        if (old_to_new[face(1)] == -1){
            old_to_new[face(1)] = indexcount;
            indexcount++;
        }
        newV.row(old_to_new[face(1)]) = pt_3.transpose();

        if (edge_to_new.find(std::make_pair(std::min(face(1), face(2)),std::max(face(1),face(2)))) == edge_to_new.end()){
            edge_to_new[std::make_pair(std::min(face(1), face(2)),std::max(face(1),face(2)))] = indexcount;
            indexcount++;
        }
        newV.row(edge_to_new[std::make_pair(std::min(face(1), face(2)),std::max(face(1),face(2)))]) = pt_4.transpose();

        if (old_to_new[face(2)] == -1){
            old_to_new[face(2)] = indexcount;
            indexcount++;
        }
        newV.row(old_to_new[face(2)]) = pt_5.transpose();

        if (edge_to_new.find(std::make_pair(std::min(face(2), face(0)),std::max(face(2),face(0)))) == edge_to_new.end()){
            edge_to_new[std::make_pair(std::min(face(2), face(0)),std::max(face(2),face(0)))] = indexcount;
            indexcount++;
        }
        newV.row(edge_to_new[std::make_pair(std::min(face(2), face(0)),std::max(face(2),face(0)))]) = pt_6.transpose();

        // New faces

        Eigen::Vector3i new_face_1 ;
        new_face_1 << old_to_new[face(0)], edge_to_new[std::make_pair(std::min(face(0), face(1)),std::max(face(0),face(1)))], edge_to_new[std::make_pair(std::min(face(2), face(0)),std::max(face(2),face(0)))];
        Eigen::Vector3i new_face_2 ;
        new_face_2 << edge_to_new[std::make_pair(std::min(face(0), face(1)),std::max(face(0),face(1)))], old_to_new[face(1)], edge_to_new[std::make_pair(std::min(face(1), face(2)),std::max(face(1),face(2)))];
        Eigen::Vector3i new_face_3 ;
        new_face_3 << edge_to_new[std::make_pair(std::min(face(1),face(2)),std::max(face(1),face(2)))], old_to_new[face(2)], edge_to_new[std::make_pair(std::min(face(2), face(0)),std::max(face(2),face(0)))];
        Eigen::Vector3i new_face_4 ;
        new_face_4 << edge_to_new[std::make_pair(std::min(face(0), face(1)),std::max(face(0),face(1)))], edge_to_new[std::make_pair(std::min(face(1), face(2)),std::max(face(1),face(2)))], edge_to_new[std::make_pair(std::min(face(2), face(0)),std::max(face(2),face(0)))];

        newF.row(4 * i+0) = new_face_1.transpose();
        // std::cout << "Added face "<< i+0 << ": " << new_face_1.transpose() << "\n";
        newF.row(4 * i+1) = new_face_2.transpose();
        //std::cout << "Added face: "<< i+1 << ": " << new_face_2.transpose() << "\n";
        newF.row(4 * i+2) = new_face_3.transpose();
        //std::cout << "Added face: "<< i+2 << ": " << new_face_3.transpose() << "\n";
        newF.row(4 * i+3) = new_face_4.transpose();
        //std::cout << "Added face: "<< i+3 << ": " << new_face_4.transpose() << "\n";

        // Relocate old vertices

        if(scheme == 1){
        for(int j=0; j<3; ++j){
            Eigen::Vector3d pt = oldV.row(face(j));
            std::vector< trimesh::index_t > neighs;
            mesh.vertex_vertex_neighbors( face(j), neighs );
            int num_neighbours = neighs.size();
            double weight_n = (64 * num_neighbours / (40 - pow((3 + 2* cos(2 * M_PI / num_neighbours)),2))) - num_neighbours;
            Eigen::Vector3d new_pt(0,0,0);
            for (int k=0; k<neighs.size(); ++k){
                new_pt += oldV.row(neighs[k]);
            }
            new_pt += weight_n * pt;
            new_pt = new_pt / (num_neighbours + weight_n);
            newV.row(old_to_new[face(j)]) = new_pt.transpose();
        }
        }


        //Relocate new vertices
        //Get neighbours of both vertices and get the intersection
        // scheme = 1 for loop subdivision, scheme = 2 for butterfly subdivision
        if (scheme == 1){
            for (int j = 0; j < 3; ++j) {
            std::vector< trimesh::index_t > neighs_1;
            std::vector< trimesh::index_t > neighs_2;

            mesh.vertex_vertex_neighbors(face(j), neighs_1);
            std::sort(neighs_1.begin(), neighs_1.end());
            mesh.vertex_vertex_neighbors(face((j + 1) % 3), neighs_2);
            std::sort(neighs_2.begin(), neighs_2.end());


            std::vector< trimesh::index_t > neighs_1_2;
            std::set_intersection(neighs_1.begin(), neighs_1.end(), neighs_2.begin(), neighs_2.end(), std::back_inserter(neighs_1_2));
            int num_neighbours_1_2 = neighs_1_2.size();

            Eigen::Vector3d pt_1_2 = 1/8.0 * (3 * oldV.row(face(j)) + 3 * oldV.row(face((j + 1) % 3)) + oldV.row(neighs_1_2[0]) + oldV.row(neighs_1_2[1]));
            newV.row(edge_to_new[std::make_pair(std::min(face(j), face((j + 1) % 3)), std::max(face(j), face((j + 1) % 3)))]) = pt_1_2.transpose();
        }
    } else{
        // Butterfly subdivision scheme
        for (int j = 0; j < 3; ++j) {
            std::vector< trimesh::index_t > neighs_1;
            std::vector< trimesh::index_t > neighs_2;

            mesh.vertex_vertex_neighbors(face(j), neighs_1);
            std::sort(neighs_1.begin(), neighs_1.end());
            mesh.vertex_vertex_neighbors(face((j + 1) % 3), neighs_2);
            std::sort(neighs_2.begin(), neighs_2.end());


            std::vector< trimesh::index_t > neighs_1_2;
            std::set_intersection(neighs_1.begin(), neighs_1.end(), neighs_2.begin(), neighs_2.end(), std::back_inserter(neighs_1_2));

            std::vector< trimesh::index_t > neighs_3;
            std::vector< trimesh::index_t > neighs_4;

            mesh.vertex_vertex_neighbors(neighs_1_2[0], neighs_3);
            std::sort(neighs_3.begin(), neighs_3.end());
            mesh.vertex_vertex_neighbors(neighs_1_2[1], neighs_4);
            std::sort(neighs_4.begin(), neighs_4.end());

            std::vector< trimesh::index_t > neighs_0_f1;
            std::set_intersection(neighs_1.begin(), neighs_1.end(), neighs_3.begin(), neighs_3.end(), std::back_inserter(neighs_0_f1));

            std::vector< trimesh::index_t > neighs_0_f2;
            std::set_intersection(neighs_2.begin(), neighs_2.end(), neighs_3.begin(), neighs_3.end(), std::back_inserter(neighs_0_f2));

            std::vector< trimesh::index_t > neighs_1_f1;
            std::set_intersection(neighs_1.begin(), neighs_1.end(), neighs_4.begin(), neighs_4.end(), std::back_inserter(neighs_1_f1));

            std::vector< trimesh::index_t > neighs_1_f2;
            std::set_intersection(neighs_2.begin(), neighs_2.end(), neighs_4.begin(), neighs_4.end(), std::back_inserter(neighs_1_f2));


            Eigen::Vector3d pt =  -(oldV.row(neighs_0_f1[0]) + oldV.row(neighs_0_f1[1]) + oldV.row(neighs_1_f1[0]) + oldV.row(neighs_1_f1[1]) - 2 * oldV.row(face((j + 1) % 3)) + \
                                  oldV.row(neighs_0_f2[0]) + oldV.row(neighs_0_f2[1])  + oldV.row(neighs_1_f2[0]) + oldV.row(neighs_1_f2[1]) - 2 * oldV.row(face(j)));


            Eigen::Vector3d pt_1_2 = pt.transpose() + (8.0 * oldV.row(face(j)) + 8.0 * oldV.row(face((j + 1) % 3)) + 2.0 * oldV.row(neighs_1_2[0]) + 2.0 * oldV.row(neighs_1_2[1]));
            pt_1_2 = pt_1_2 / 16.0;

            newV.row(edge_to_new[std::make_pair(std::min(face(j), face((j + 1) % 3)), std::max(face(j), face((j + 1) % 3)))]) = pt_1_2.transpose();
    }
        }
        
    
    // std::cout<< "Vertices \n" << "\n" << newV << "\n";
    // std::cout<< "Faces \n" << newF << "\n";

    }

}



int main(int argc, char *argv[])
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::read_triangle_mesh("../../input/bunny.obj", V,F);

    // half-edges example
    // std::vector< trimesh::triangle_t > triangles;

    // int kNumVertices = V.rows();
    // int kNumFaces = F.rows();
    // triangles.resize( kNumFaces );
    // for (int i=0; i<kNumFaces; ++i){
    //     triangles[i].v[0] = F(i,0);
    //     triangles[i].v[1] = F(i,1);
    //     triangles[i].v[2] = F(i,2);
    // }

    // std::vector< trimesh::edge_t > edges;
    // trimesh::unordered_edges_from_triangles( triangles.size(), &triangles[0], edges );
    
    // trimesh::trimesh_t mesh;
    // mesh.build( kNumVertices, triangles.size(), &triangles[0], edges.size(), &edges[0] );

    // std::vector< trimesh::index_t > neighs;
    // for( int vi = 0; vi < kNumVertices; ++vi )
    // {
    //     mesh.vertex_vertex_neighbors( vi, neighs );

    //     std::cout << "neighbors of vertex " << vi << ": ";
    //     for( int i = 0; i < neighs.size(); ++i )
    //     {
    //         std::cout << ' ' << neighs.at(i);
    //     }
    //     std::cout << '\n';
    // }

    
    //! Has to implement the subdivision scheme here 
    for (int i=0; i<1 ; ++i){

        int kNumVertices = V.rows();
        int kNumFaces = F.rows();
        std::vector< trimesh::triangle_t > triangles;
        triangles.resize( kNumFaces );
        for (int i=0; i<kNumFaces; ++i){
            triangles[i].v[0] = F(i,0);
            triangles[i].v[1] = F(i,1);
            triangles[i].v[2] = F(i,2);
        }

        std::vector< trimesh::edge_t > edges;
        trimesh::unordered_edges_from_triangles( triangles.size(), &triangles[0], edges );

        trimesh::trimesh_t mesh;
        mesh.build( kNumVertices, triangles.size(), &triangles[0], edges.size(), &edges[0] );
        Eigen::MatrixXd newV(V.rows() + edges.size(), V.cols());
        Eigen::MatrixXi newF(F.rows()*4, F.cols());
        subdivide(mesh, V, F,edges, newV, newF);
        V.resize(newV.rows(), newV.cols());
        F.resize(newF.rows(), newF.cols());
        V = newV ;
        F = newF;
    }
    

    // Eigen::MatrixXi newF;
    // newF = mesh.get_faces();
    // std::cout<< "\n" << newF << "\n";


        // Plot the mesh
    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);
    viewer.data().set_face_based(true);
    // const Eigen::RowVector3d red = Eigen::RowVector3d(1,0,0);
    // add vertices highlights
    // viewer.data().point_size = 5;
    // viewer.data().add_points(V, red);
    // add vertices index
    // for (int i=0; i<V.rows(); ++i){
    //     viewer.data().add_label(V.row(i)+Eigen::RowVector3d(0.005, 0.005, 0),std::to_string(i));
    // }
    viewer.data().show_custom_labels = true;
    // launch viewer
    viewer.launch();

    
    // newV = mesh.get_vertices();
    // std::cout<< "\n" << newV << "\n";
    // // output the mesh
    // igl::writeOBJ("../output/cube.obj", newV, newF);

    // Plot the mesh
    // igl::opengl::glfw::Viewer viewer;
    // viewer.data().set_mesh(V, newF);
    // viewer.data().set_face_based(true);

    // const Eigen::RowVector3d red = Eigen::RowVector3d(1,0,0);
    // // add vertices highlights
    // viewer.data().point_size = 5;
    // viewer.data().add_points(V, red);
    // // add vertices index
    // for (int i=0; i<V.rows(); ++i){
    //     viewer.data().add_label(V.row(i)+Eigen::RowVector3d(0.005, 0.005, 0),std::to_string(i));
    // }
    // viewer.data().show_custom_labels = true;
    // // launch viewer
    // viewer.launch();
}
