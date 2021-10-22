#include <iostream>
#include <open3d/Open3D.h>

int main(int argc, char* argv[]){
    auto sphere = open3d::geometry::TriangleMesh::CreateSphere(1.0);
    sphere->ComputeVertexNormals();
    sphere->PaintUniformColor({0.0, 1.0, 0.0});
    open3d::visualization::DrawGeometries({sphere});
    return 0;
}