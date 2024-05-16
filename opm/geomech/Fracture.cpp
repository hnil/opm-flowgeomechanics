#include "config.h"
#include <opm/geomech/Fracture.hpp>
#include <opm/geomech/Math.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/geomech/DiscreteDisplacement.hpp>
namespace Opm
{
void
Fracture::init(std::string well, int perf, int well_cell, Fracture::Point3D origo, Fracture::Point3D normal)
{
    wellinfo_ = WellInfo({well, perf, well_cell});
    origo_ = origo;
    axis_[2] = normal;
    axis_[0] = Point3D({std::copysign(normal[2], normal[0]),
                        std::copysign(normal[2], normal[1]),
                        -std::copysign(std::abs(normal[0]) + std::abs(normal[1]), normal[2])});
    axis_[1] = crossProduct(axis_[2], axis_[0]);
    for (int i = 0; i < 3; ++i) {
        axis_[i] /= axis_[i].two_norm();
        axis_[i] *= 20;
    }
    layers_ = 0;
    initFracture();
    grow(1, 0);
    nlinear_ = layers_;
    //grow(10, 1);
    grid_->grow();
    grid_->postGrow();
    vtkwriter_ = std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(grid_->leafGridView(), Dune::VTK::nonconforming);
    //

}
void Fracture::setupPressureSolver(){
        Opm::FlowLinearSolverParameters p;
    prmpressure_ = Opm::setupPropertyTree(p, true, true);
    {
        std::size_t pressureIndex; // Dummy
        const std::function<Vector()> weightsCalculator; // Dummy
        auto pressure_operator = std::make_unique<PressureOperatorType>(*pressure_matrix_);
        // using FlexibleSolverType = Dune::FlexibleSolver<SeqOperatorType>;
        pressure_operator_ = std::move(pressure_operator);
        auto psolver
            = std::make_unique<FlexibleSolverType>(*pressure_operator_, prmpressure_, weightsCalculator, pressureIndex);

        pressure_solver_ = std::move(psolver);
    }

}
std::string
Fracture::name() const
{
    std::string name = "Fracure_on_" + wellinfo_.name + "_perf_" + std::to_string(wellinfo_.perf);
    return name;
}
void
Fracture::initFracture()
{
    Dune::GridFactory<Grid> factory; // Dune::FoamGrid<2,3>> factory;
    size_t N = 6;
    double radius = 1;
    std::vector<unsigned int> inner_indices;
    std::vector<Dune::FieldVector<double, 3>> vertices;
    {
        // Dune::FieldVector<double, 3> vertex({0, 0, 0});
        Point3D vertex = surfaceMap(0.0, 0.0);
        vertices.push_back(vertex);
    }

    for (size_t i = 0; i < N; ++i) {

        inner_indices.push_back(i + 1);
        double theta = (i * 2 * M_PI) / N;
        double x = radius * cos(theta);
        double y = radius * sin(theta);
        Point3D vertex = surfaceMap(x, y);
        vertices.push_back(vertex);
        // assume the first 6 cells has source
        well_source_.push_back(i);
    }
    for (size_t i = 0; i < vertices.size(); i++) {
        factory.insertVertex(vertices[i]);
    }
    std::vector<std::vector<unsigned int>> cornerIDs;
    for (size_t i = 0; i < N; ++i) {
        unsigned int next = inner_indices[(i + 1) % N];
        std::vector<unsigned int> cornerID({unsigned(0), unsigned(inner_indices[i]), next});
        cornerIDs.push_back(cornerID);
    }
    for (size_t i = 0; i < N; ++i) {
        factory.insertElement(Dune::GeometryTypes::simplex(2), cornerIDs[i]); // std::move(mapping));
    }
    out_indices_ = inner_indices;
    grid_ = factory.createGrid();
    // grid_ = factory.createGrid();
}

void
Fracture::write() const
{
    if (reservoir_cells_.size() > 0) {
        vtkwriter_->addCellData(reservoir_cells_, "ReservoirCell");
    }
    if (reservoir_perm_.size() > 0) {
        vtkwriter_->addCellData(reservoir_perm_, "ReservoirPerm");
    }
    if (reservoir_dist_.size() > 0) {
        vtkwriter_->addCellData(reservoir_dist_, "ReservoirDist");
    }
    if (fracture_pressure_.size() > 0) {
        vtkwriter_->addCellData(fracture_pressure_, "FracturePressure");
    }
    if (reservoir_pressure_.size() > 0) {
        vtkwriter_->addCellData(fracture_pressure_, "ReservoirPressure");
    }
    if (fracture_width_.size() > 0) {
        vtkwriter_->addCellData(fracture_width_, "FractureWidth");
    }
    vtkwriter_->write(this->name().c_str());
};
void
Fracture::grow(int layers, int method)
{
    while (layers_ < layers) {
        std::vector<unsigned int> inner_indices = out_indices_;
        out_indices_.resize(0);
        if (method == 0) {
            this->insertExp(inner_indices);
        } else {
            this->insertLinear(inner_indices);
            nlinear_ += 1;
        }
        ++layers_;
        inner_indices = out_indices_;
    }
}
void
Fracture::insertLinear(const std::vector<unsigned int>& inner_indices)
{
    double radius = 1.0;
    size_t N = inner_indices.size();
    for (size_t i = 0; i < N; ++i) {
        double theta = (i * 2 * M_PI) / (N);
        theta += (layers_ - nlinear_) * 0.5 * (2 * M_PI / N);
        double out_radius = radius + layers_ + 1;
        double x = out_radius * cos(theta);
        double y = out_radius * sin(theta);
        Point3D vertex = surfaceMap(x, y);
        int new_ind = grid_->insertVertex(vertex);
        out_indices_.push_back(new_ind);
    }
    for (size_t i = 0; i < N; ++i) {
        {
            std::vector<unsigned int> cornerID(
                {inner_indices[i % N], out_indices_[(i) % (N)], inner_indices[(i + 1) % (N)]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }

        {
            std::vector<unsigned int> cornerID(
                {inner_indices[(i + 1) % N], out_indices_[(i) % (N)], out_indices_[(i + 1) % (N)]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }
    }
}

void
Fracture::insertExp(const std::vector<unsigned int>& inner_indices)
{
    size_t N = inner_indices.size();
    double radius = 1.0;

    for (size_t i = 0; i < N * 2; ++i) {
        double theta = (i * 2 * M_PI) / (N * 2);
        double out_radius = radius + layers_ + 1;
        double x = out_radius * cos(theta);
        double y = out_radius * sin(theta);
        Point3D vertex = surfaceMap(x, y);
        int new_ind = grid_->insertVertex(vertex);
        out_indices_.push_back(new_ind);
    }


    for (size_t i = 0; i < N; ++i) {
        {
            std::vector<unsigned int> cornerID(
                {inner_indices[i], out_indices_[(2 * i) % (N * 2)], out_indices_[(2 * i + 1) % (N * 2)]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }
        {
            std::vector<unsigned int> cornerID(
                {inner_indices[i], out_indices_[(2 * i + 1) % (N * 2)], inner_indices[(i + 1) % N]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }
        {
            std::vector<unsigned int> cornerID(
                {inner_indices[(i + 1) % N], out_indices_[(2 * i + 1) % (N * 2)], out_indices_[(2 * i + 2) % (N * 2)]});
            grid_->insertElement(Dune::GeometryTypes::simplex(2), cornerID);
        }
    }
}

Dune::FieldVector<double, 3>
Fracture::surfaceMap(double x, double y)
{
    Point3D vec(0.0);
    vec += x * axis_[0];
    vec += y * axis_[1];
    // for (int i = 0; i < 2; ++i) {
    // }
    vec += origo_;
    return vec;
}
template <class Grid3D>
void
Fracture::updateReservoirCells(const external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree, Grid3D& grid3D)
{
    reservoir_cells_.resize(grid_->leafGridView().size(0));
    using GridView = typename Grid::LeafGridView;
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    ElementMapper elemMapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    int tri_divide = 0;
    int tri_outside = 0;
    for (auto& element : elements(grid_->leafGridView())) {
        const auto elemIdx = elemMapper.index(element);
        auto geom = element.geometry();
        external::cvf::BoundingBox bb;
        using Vec3d = external::cvf::Vec3d;
        // if(false)
        // {
        // Point3D center = geom.center();
        // double eps = 1e-3;
        // Vec3d p1(center[0],center[1],center[2]);
        // external::cvf::Vec3d  p2(center[0]+eps,center[1]+eps,center[2]+eps);
        // bb.add(p1);
        // bb.add(p2);
        // }else
        {
            std::array<Vec3d, 3> corners;
            for (int i = 0; i < geom.corners(); ++i) {
                auto vertex = geom.corner(i);
                Vec3d point(vertex[0], vertex[1], vertex[2]);
                bb.add(point);
            }
        }

        std::vector<size_t> cells = external::findCloseCellIndices(cellSearchTree, bb);
        if (cells.size() > 0) {
            reservoir_cells_[elemIdx] = cells[0];
            if (cells.size() > 1) {
                // std::vector<HexIntersectionInfo> intersections;
                // for ( const auto& globalCellIndex : cells )
                // {
                //     cvf::Vec3d    hexCorners[8];  //resinsight numbering, see RigCellGeometryTools.cpp
                //     RigEclipseWellLogExtractor::hexCornersOpmToResinsight( hexCorners, globalCellIndex);
                //     RigHexIntersectionTools::lineHexCellIntersection( p1, p2, hexCorners, globalCellIndex,
                //     &intersections );
                // }
                tri_divide += 1;
            }
        } else {
            tri_outside += 1;
            reservoir_cells_[elemIdx] = -1;
            // std::cout << "Fallback Search" << std::endl;
            // int cell = Opm::findCell(grid3D,geom.center());
        }
    }
    std::cout << "For Fracture : " << this->name() << " : " << tri_divide << " triangles should be devided"
              << std::endl;
    std::cout << "For Fracture : " << this->name() << " : " << tri_outside << " triangles outside" << std::endl;
    std::cout << "Total triangles: " << grid_->leafGridView().size(0) << std::endl;
}


void
Fracture::updateReservoirProperties()
{
    double perm = 1e-13;
    size_t nc = grid_->leafGridView().size(0);
    assert(reservoir_cells_.size() == nc);
    reservoir_perm_.resize(nc, perm);
    reservoir_dist_.resize(nc, 10.0);
    reservoir_pressure_.resize(nc, 1.0);
    this->initFractureWidth();
}
void
Fracture::solve()
{
    this->solveFractureWidth();
    // grow fracture
    this->setSource();
    this->solvePressure();
}
void
Fracture::setSource()
{
    if(rhs_pressure_.size() == 0){
        size_t nc = grid_->leafGridView().size(0);
        rhs_pressure_.resize(nc);
    }
    rhs_pressure_ = 0;
    double scale = well_source_.size();
    for (auto cell : well_source_) {
        rhs_pressure_[cell] += 1.0 / scale;
    }
}
void
Fracture::solvePressure()
{
    size_t nc = grid_->leafGridView().size(0);
    fracture_pressure_.resize(nc);
    fracture_pressure_ = 1e5;
    if(!pressure_matrix_){
        this->initPressureMatrix();
    }
    this->assemblePressure();
    this->setSource(); // probably include reservoir pressure
    try {
        if(!pressure_solver_){
            this->setupPressureSolver();
        }
        Dune::InverseOperatorResult r;
        fracture_pressure_.resize(rhs_pressure_.size());
        fracture_pressure_ = 0;
        pressure_solver_->apply(fracture_pressure_, rhs_pressure_, r);
    } catch (Dune::ISTLError& e) {
        std::cerr << "exception thrown " << e << std::endl;
    }
}
void
Fracture::solveFractureWidth()
{
    this->assembleFracture();
    fracture_matrix_->solve(fracture_width_,rhs_width_);
    for(int i=0; i < fracture_width_.size(); ++i){
        assert(fracture_width_[i]>0);
    }
}

void
Fracture::initFractureWidth()
{
    size_t nc = grid_->leafGridView().size(0);
    fracture_width_.resize(nc);
    fracture_width_ = 1e-3;
    ElementMapper elemMapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto& element : elements(grid_->leafGridView())) {
        const auto elemIdx = elemMapper.index(element);
        auto geom = element.geometry();
        auto vec_origo = geom.center() - origo_;
        double dist_origo = vec_origo.two_norm();
        fracture_width_[elemIdx] *= dist_origo;
    }
}
void
Fracture::initPressureMatrix()
{
    // size_t num_columns = 0;
    //  index of the neighbour
    size_t nc = grid_->leafGridView().size(0);
    leakof_.resize(nc,0.0);
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto& element : Dune::elements(grid_->leafGridView())) {
        int eIdx = mapper.index(element);
        auto geom = element.geometry();
        auto eCenter = geom.center();
        // iterator over all intersections
        for (auto& is : Dune::intersections(grid_->leafGridView(),element)) {
            if (!is.boundary()) {
                int nIdx = mapper.index(is.outside());
                if (eIdx < nIdx) {
                    // calculate distance between the midpoints
                    auto nCenter = is.outside().geometry().center();
                    auto isCenter = is.geometry().center();
                    eCenter -= isCenter;
                    nCenter -= isCenter;
                    // probably should use projected distenace
                    auto igeom = is.geometry();
                    double area = igeom.volume();
                    double h1 = eCenter.two_norm() / area;
                    double h2 = nCenter.two_norm() / area;

                    {
                        Htrans matel(nIdx, eIdx, h1, h2);
                        htrans_.push_back(matel);
                    }
                }
            }
        }
        {
            double value = (reservoir_perm_[eIdx] * geom.volume()) ;
            value /= reservoir_dist_[eIdx];
            leakof_[eIdx] = value;
        }
    }

    // not need if build mode is implicit
    // std::sort(transes.begin(),transes.end(),sortcsr);
    //  build matrix
    //if (pressure_matrix_.bu == 0){
    pressure_matrix_ = std::make_unique<Matrix>();
        auto& matrix = *pressure_matrix_;
        matrix.setBuildMode(Matrix::implicit);
        // map from dof=3*nodes at a cell (ca 3*3*3) to cell
        matrix.setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
        //size_t nc = grid_->leafGridView().size(0);
        matrix.setSize(nc, nc);
        for (auto matel : htrans_) {
            size_t i = std::get<0>(matel);
            size_t j = std::get<1>(matel);
            matrix.entry(i, j) = 0;
            matrix.entry(j, i) = 0;
            matrix.entry(j, j) = 0;
            matrix.entry(i, i) = 0;
        }
        matrix.compress();
    //}
}
void
Fracture::assemblePressure()
{
    auto& matrix = *pressure_matrix_;
    for (auto matel : htrans_) {
        size_t i = std::get<0>(matel);
        size_t j = std::get<1>(matel);
        double t1 = std::get<2>(matel);
        double t2 = std::get<3>(matel);
        double h1 = fracture_width_[i];
        double h2 = fracture_width_[j];
        // harmonic mean of surface flow
        double value = 12. / (h1 * h1 * t1) + 1. / (h2 * h2 * t2);
        value = 1 / value;

        // matrix.entry(i, j) -= value;
        // matrix.entry(j, i) -= value;
        // //
        // matrix.entry(i, i) += value;
        // matrix.entry(j, j) += value;

        matrix[i][j] -= value;
        matrix[j][i] -= value;
        //
        matrix[i][i] += value;
        matrix[j][j] += value;
    }
    for (size_t i = 0; i < leakof_.size(); ++i) {
        //matrix.entry(i, i) += leakof_[i];
        matrix[i][i] += leakof_[i];
    }
}

void  Fracture::assembleFracture(){
    size_t nc = grid_->leafGridView().size(0);
    rhs_width_.resize(nc);
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto elem : elements(grid_->leafGridView())) {
        int idx = mapper.index(elem);
        //auto geom = elem.geometry();
        rhs_width_[idx] = reservoir_pressure_[idx];//*geom.volume();
    }

    if(!fracture_matrix_){
        fracture_matrix_ = std::make_unique<DynamicMatrix>();
    }
    fracture_matrix_->resize(nc,nc);
    double E = 1;
    double nu = 0.3;
    ddm::assembleMatrix(*fracture_matrix_,E,nu,*grid_);
}

// bool sortcsr(const tuple<size_t,size_t,double>& a,
//              const tuple<size_t,size_t,double>& b){
//     if(get<0>(a) < get<0>(b)){
//         return true;
//     }else{
//         if(get<1>(a) < get<1>(b)){
//             return true;
//         }else{
//             return false;
//         }
//     }
// }

// void buildRowWise(const std::tuple<size_t,size_t,double>& transes,Matrix& matrix){
//     matrix.setBuildMode(Matrx::row_wise);
//     matrix.setSize(1,nc,nc);
//     size_t row_index = 0;
//     auto trans_it = tranes.begin();
//     for (auto row = matrix.createbegin(),
//              end = matrix.createend(); row != end; ++row) {
//         while(row.index() == get<0>(*trans_it)){
//             row.insert(get<1>(*trans_it))
//                 ++trans_it;
//         }else{
//             ++row;
//         }
//     }
// }



template void
Fracture::updateReservoirCells<Dune::CpGrid>(const external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree,
                                             Dune::CpGrid& grid3D);
} // namespace Opm
