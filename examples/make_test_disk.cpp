#include <gmsh.h>

int main()
{
  gmsh::initialize();
  gmsh::model::add("disk");

  const double c[] = {0.0, 0.0, 0.0};
  const double r = 1;
  gmsh::model::occ::addDisk(c[0], c[1], c[2], r, r);
  gmsh::model::occ::synchronize();

  const double meshsize = 0.4;
  gmsh::option::setNumber("Mesh.MeshSizeMax", meshsize);
  gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);
  
  gmsh::model::mesh::generate(2);
  gmsh::write("disk.msh");
  gmsh::write("disk.m");
}
