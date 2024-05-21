#include <gmsh.h>

int main()
{
  gmsh::initialize();
  gmsh::model::add("disk");

  const double c[] = {0.0, 0.0, 0.0};
  const double r = 1;
  gmsh::model::occ::addDisk(c[0], c[1], c[2], r, r);
  gmsh::model::occ::synchronize();

  const double meshsize = 0.1;
  gmsh::option::setNumber("Mesh.MeshSizeMax", meshsize);
  
  gmsh::model::mesh::generate(2);
  gmsh::write("disk.mesh");
  gmsh::write("disk.m");
}
