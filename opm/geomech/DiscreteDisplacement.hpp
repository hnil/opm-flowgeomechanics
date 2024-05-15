#pragma once
namespace ddm{
    struct Real3;
    struct Real6;
    Real3 make3(double x, double, y double z);
    Real6 make6(double x, double, y double z, double a, double, b, double c);
template< class Element>
Dune::FieldVector<double, 3>
TDDispFS(const Dune::FieldVector<double 3>& obs,
        const Element& elem, const Dune::FieldVector<double, 3>& slip,double nu){
    std::array<int,3> tri;
    for(int i=0; i < elem.corners(); ++i){
        tri[i] = elem.corners(i);
    }
    Real3 obs_tmp = make3(obs[0],obs[1], obs[2]);
    Real3 slip_tmp = make3(slip[0],slip[1], slip[2]);
    Real3 disp_tmp = disp_fs(obs_tmp, tri, slip_tmp, nu);
    Dune::FieldVector<double, 3> disp;
    disp[0] = disp_tmp.x;
    disp[1] = disp_tmp.y;
    disp[2] = disp_tmp.z;
    return stress
}
TDStrainFS(const Dune::FieldVector<double 3>& obs,
        const Element& elem, const Dune::FieldVector<double, 3>& slip,double nu){
    std::array<int,3> tri;
    for(int i=0; i < elem.corners(); ++i){
        tri[i] = elem.corners(i);
    }
    Real3 obs_tmp = make3(obs[0],obs[1], obs[2]);
    Real3 slip_tmp = make3(slip[0],slip[1], slip[2]);
    Real6 strain_tmp = stress_fs(obs_tmp, tri, slip_tmp, nu);
    Dune::FieldVector<double, 3> strain;
    strain[0] = strain_tmp.x;
    strain[1] = strain_tmp.y;
    strain[2] = strain_tmp.z;
    strain[3] = strain_tmp.a;
    strain[4] = strain_tmp.b;
    strain[5] = strain_tmp.c;
    return strain;
}
    void makeMatrix(double E, double nu, Dune::FoamGrid<3,2> grid){
        for(auto elem1: elements(grid.leaftriView())){
            for(auto elem2: elements(grid.leaftriView())){
                auto geom = elem1.geom();
                auto center = geom.center();
                auto strain = TDStrainFs(center, elem2, slip, nu);
                //auto stress =
                // auto traction
                double matel = 1;
                matrix[idx1][idx2] = 1;
            }
        }
    }

}
