#include <opm/geomech/DiscreteDisplacement.hpp>
namespace ddm{
double fractureK1(double dist,double width, double E, double nu){
    double K1;
    double mu = E/(2*(1+nu));//??
    K1 = (mu*sqrt(M_PI)/(2*std::sqrt(dist)*(1.0-nu)))*width;
    K1 /= 1.834; //factor found numerically //@@ Odd: changed this from '*' to '/'
    return K1;
}
/*
//formulas for share ned rotation to find share in correct directions
double fractureK2(double dist,double edist, double E, double nu){
    double K1;
    K1 = mu*sqrt(M_PI)/(2*sqrt(dist).*(1-nu))*emid;
    K1 *= 1.834; //factor found numerically
    return K1;
}

double fractureK3(double dist,double edist, double E, double nu){
    double K1;
    K1 = mu*sqrt(M_PI)/(2*sqrt(dist).*(1-nu))*edist*(1-nu);
    K1 *= 1.834; //factor found numerically
    return K1;
}
*/
double traceSymTensor(const Dune::FieldVector<double,6> symtensor){
  double trace = 0;
  for(int i=0; i < 3; ++i){
    trace += symtensor[i];
  }
  return trace;
}

Dune::FieldVector<double,6>
strainToStress(const double E,const double nu,const Dune::FieldVector<double,6> strain){
    Dune::FieldVector<double,6> stress;
    double volume_strain = traceSymTensor(strain);
    double mu = E/(2*(1+nu));//??
    double lambda = 2 * mu * nu / (1 - 2 * nu);
    for(int i=0; i < 3; ++i){
        stress[i]   =  2*mu*strain[i] + (lambda * volume_strain);
        stress[i+3] += 2*mu*strain[i+3]; // [xy,xz,yz]
    }
    return stress;
}

Dune::FieldMatrix<double,3,3>
symTensor2Matrix(const Dune::FieldVector<double,6> symtensor){
    Dune::FieldMatrix<double,3,3> mat;
    for(int i=0; i < 3; ++i){
        mat[i][i] = symtensor[i];
    }
    mat[2][1] = symtensor[3]; mat[1][2] = mat[2][1];
    mat[2][0] = symtensor[4]; mat[0][2] = mat[2][0];
    mat[1][0] = symtensor[5]; mat[0][1] = mat[1][0];
    return mat;
}


double
tractionSymTensor(const Dune::FieldVector<double,6> symtensor, Dune::FieldVector<double,3> normal){
    const Dune::FieldVector<double,6> stress;
    assert(std::abs(normal.two_norm()-1) < 1e-13);
    double traction=0.0;
    for(int i=0; i < 3; ++i){
        traction += symtensor[i]*normal[i]*normal[i];
    }
    traction += 2*symtensor[3]*normal[0]*normal[1];// xy*nx*ny;
    traction += 2*symtensor[4]*normal[0]*normal[2];// xz*nx*nz
    traction += 2*symtensor[5]*normal[1]*normal[2];// yz*ny*nz
    return traction;
}
void
//assembleMatrix(Dune::DynamicMatrix<Dune::FieldMatrix<double,1,1>>& matrix, const double E, const double nu, const Dune::FoamGrid<2, 3>& grid)
assembleMatrix(Dune::DynamicMatrix<double>& matrix, const double E, const double nu, const Dune::FoamGrid<2, 3>& grid)

{
    using Grid = Dune::FoamGrid<2, 3>;
    using GridView = typename Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
    for (auto elem1 : elements(grid.leafGridView())) {
        int idx1 = mapper.index(elem1);
        auto geom = elem1.geometry();
        auto center = geom.center();
        auto normal = normalOfElement(elem1);
        for (auto elem2 : elements(grid.leafGridView())) {
            int idx2 = mapper.index(elem2);
            // check if this is defined in relative coordinates
            Dune::FieldVector<double,3> slip;// = make3(1.0,0.0, 0.0);
            slip[0]=1;slip[1]=0;slip[2]=0;

            // symmetric stress voit notation
            Dune::FieldVector<double,6> strain = TDStrainFS(center, elem2, slip, nu);
            Dune::FieldVector<double,6> stress = strainToStress(E,nu,strain);

            double ntraction = tractionSymTensor(stress,normal);
            // matrix relate to pure traction not area weighted
            matrix[idx1][idx2] = ntraction;
        }
    }
}

}
