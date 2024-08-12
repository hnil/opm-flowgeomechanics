#include "config.h"
#include <opm/grid/polyhedralgrid.hh>
#include <opm/geomech/Fracture.hpp>
#include <opm/geomech/Math.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/geomech/DiscreteDisplacement.hpp>
#include <dune/common/filledarray.hh> // needed for printSparseMatrix??
#include <dune/istl/io.hh> // needed for printSparseMatrix??
#include <cmath>
#include <iostream>
#include <sstream>
#include<string>
#include <dune/common/fmatrixev.hh>
#include <dune/grid/utility/persistentcontainer.hh>

#include <opm/grid/polyhedralgrid.hh>

namespace Opm
{
void
Fracture::init(std::string well,
               int perf,
               int well_cell,
               Fracture::Point3D origo,
               Fracture::Point3D normal,
               Opm::PropertyTree prm
    )
{
    prm_ = prm;
    wellinfo_ = WellInfo({well, perf, well_cell});
    origo_ = origo;
    axis_[2] = normal;
    axis_[0] = Point3D({std::copysign(normal[2], normal[0]),
                        std::copysign(normal[2], normal[1]),
                        -std::copysign(std::abs(normal[0]) + std::abs(normal[1]), normal[2])});
    axis_[1] = crossProduct(axis_[2], axis_[0]);
    double init_scale = prm_.get<double>("config.axis_scale");
    for (int i = 0; i < 3; ++i) {
        axis_[i] /= axis_[i].two_norm();
        axis_[i] *= init_scale;
    }
    layers_ = 0;
    nlinear_ = 0;
    initFracture();
    grow(1, 0);
    nlinear_ = layers_;
    //grow(4, 1);
    grid_->grow();
    grid_->postGrow();

    size_t nc = grid_->leafGridView().size(0);
    reservoir_cells_ = std::vector<int>(nc, -1);
    fracture_pressure_.resize(nc); fracture_pressure_ = 1e5;
    this->initFractureWidth();
    //fracture_width_.resize(nc); fracture_width_ = 1e-3;
    //
    this->resetWriters();
}

void Fracture::resetWriters(){
    // nead to be reseat if grid is changed ??
    vtkwriter_ = std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(grid_->leafGridView(), Dune::VTK::nonconforming);

    std::string outputdir = prm_.get<std::string>("outputdir");
    std::string simName = prm_.get<std::string>("casename") + this->name();
    std::string multiFileName =  "";
    vtkmultiwriter_ = std::make_unique< Opm::VtkMultiWriter<Grid::LeafGridView, VTKFormat > >(/*async*/ false,
                                                                                              grid_->leafGridView(),
                                                                                              outputdir,
                                                                                              simName,
                                                                                              multiFileName
        );

}

void Fracture::setupPressureSolver(){
    Opm::FlowLinearSolverParameters p;
    p.linsolver_ = prm_.get<std::string>("pressuresolver");
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
void Fracture::removeCells(){
    // copy all to presistent container
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout()); // used id sets interally
    Dune::PersistentContainer<Grid,double> fracture_width(*grid_,0);
    fracture_width.resize();
    for(const auto& elem: elements(grid_->leafGridView())){
        size_t eIdx = mapper.index(elem);
        fracture_width[elem] = fracture_width_[eIdx];
    }

    //const auto& indexSet = foamGridLeafView.indexSet();// for indices
    //const auto& indexSet = grid.localIdSet();// presitent numbering
    for(const auto& elem: elements(grid_->leafGridView())){
        size_t eIdx = mapper.index(elem);
        if(reservoir_cells_[eIdx] < 0){
            grid_->removeElement(elem);
        }
    }
    grid_->grow();
    grid_->postGrow();
    // resize array

    fracture_width_.resize(grid_->leafGridView().size(0));
    // copy back from presistent contatiner
    for(const auto& elem: elements(grid_->leafGridView())){
        size_t eIdx = mapper.index(elem);
        fracture_width_[eIdx] = fracture_width[elem];
    }
    this->resetWriters();
}
std::string
Fracture::name() const
{
    std::string name = "Fracure_on_" + wellinfo_.name + "_perf_" + std::to_string(wellinfo_.perf) + "_nr" ;
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

std::vector<double> Fracture::stressIntensityK1() const{
        size_t nc = grid_->leafGridView().size(0);
        std::vector<double> stressIntensityK1(nc, std::nan("0"));
        ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
        for(const auto& elem: elements(grid_->leafGridView())){
            //bool isboundary = false;
            for (auto& is : Dune::intersections(grid_->leafGridView(),elem)) {
                if (is.boundary()) {
                    int nIdx = mapper.index(elem);
                    auto isCenter = is.geometry().center();
                    auto elCenter = elem.geometry().center();
                    auto vecC = isCenter-elCenter;
                    auto distC = vecC.two_norm();
                    double K1 = ddm::fractureK1(distC, fracture_width_[nIdx], E_, nu_);
                    stressIntensityK1[nIdx] = K1;
                }
            }

        }
        return stressIntensityK1;
    }

void
Fracture::write(int reportStep) const
{
    std::vector<double> K1;//need to be in scope until written
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
    if (reservoir_mobility_.size() > 0) {
        vtkwriter_->addCellData(reservoir_mobility_, "ReservoirMobility");
    }
    if (fracture_width_.size() > 0) {
        vtkwriter_->addCellData(fracture_width_, "FractureWidth");
        K1 = this->stressIntensityK1();
        vtkwriter_->addCellData(K1, "stressIntensityK1");
    }
    //std::stringstream ss(this->name());
    std::string outputdir = prm_.get<std::string>("outputdir");
    std::string simName = prm_.get<std::string>("casename");
    std::string filename = outputdir + "/Static_" + simName + this->name();
    if(reportStep > 0){
        filename = filename + "_step_" +  std::to_string(reportStep);
    }
    vtkwriter_->write(filename.c_str());
};

void Fracture::writemulti(double time) const
{
    //vtkmultiwriter_->gridChanged(); need to be called if grid is changed
    // need to have copies in case of async outout (and interface to functions)
    std::vector<double> K1 = this->stressIntensityK1();
    std::vector<double> fracture_pressure(fracture_pressure_.size(),0.0);
    std::vector<double> reservoir_pressure = reservoir_pressure_;//.size(),0.0);
    std::vector<double> reservoir_dist = reservoir_dist_;//.size(),0.0);
    std::vector<double> reservoir_perm = reservoir_perm_;//.size(),0.0);
    std::vector<double> reservoir_mobility = reservoir_mobility_;//.size(),0.0);
    std::vector<double> reservoir_traction(reservoir_stress_.size(),0);
    std::vector<double> reservoir_cells(reservoir_cells_.size(),0.0);
    std::vector<double> fracture_width(fracture_width_.size(),0.0);
    std::vector<double> rhs_width(rhs_width_.size(),0.0);
    std::vector<double> well_index(rhs_pressure_.size(),0.0);
    std::map<int,double> wellIndMap;
    std::vector<double> leakofrate = leakOfRate();
    //assert(leakofrate.size() == fracture_pressure_.size());
    {
        //make map to do it easy
        auto wellIndices = this->wellIndices();
        for(const auto& wind: wellIndices){
            int res_cell = std::get<0>(wind);
            double WI = std::get<1>(wind);
            wellIndMap[res_cell] = WI;
        }
    }
    // loop for need things only with fracture
    for(size_t i=0; i < rhs_width_.size(); ++i){
        rhs_width[i] = rhs_width_[i][0];
    }

    for(size_t i=0; i < fracture_width_.size(); ++i){
        fracture_width[i] = fracture_width_[i][0];
        fracture_pressure[i] = fracture_pressure_[i][0];
        reservoir_cells[i] = reservoir_cells_[i];// only converts to double
        reservoir_traction[i] = ddm::tractionSymTensor(reservoir_stress_[i],cell_normals_[i]);
    }

    for(size_t i=0; i < rhs_pressure_.size(); ++i){
         // maybe slow
        well_index[i] = wellIndMap[reservoir_cells_[i]];
    }
    
    vtkmultiwriter_->beginWrite(time);


    if (reservoir_cells.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_cells, "ReservoirCell");
    }
    if (reservoir_perm.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_perm, "ReservoirPerm");
    }
    if (reservoir_dist.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_dist, "ReservoirDist");
    }
    if (fracture_pressure.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(fracture_pressure, "FracturePressure");
    }
    if (reservoir_pressure.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_pressure, "ReservoirPressure");
    }
    if (reservoir_traction.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_traction, "ReservoirTraction");
    }
    if (reservoir_mobility.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(reservoir_mobility, "ReservoirMobility");
    }
    if (rhs_width.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(rhs_width, "ForceOnFracture");
    }
    if (leakofrate.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(leakofrate, "LeakOfRate");
    }

    if (well_index.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(well_index, "WellIndex");
    }
    if (fracture_width.size() > 0) {
        vtkmultiwriter_->attachScalarElementData(fracture_width, "FractureWidth");

        vtkmultiwriter_->attachScalarElementData(K1, "stressIntensityK1");
    }
    vtkmultiwriter_->endWrite();
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
Fracture::updateReservoirCells(const external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree,const Grid3D& grid3D)
{
    cell_normals_.resize(grid_->leafGridView().size(0));
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
            /*
            std::array<Vec3d, 3> corners;
            for (int i = 0; i < geom.corners(); ++i) {
                auto vertex = geom.corner(i);
                Vec3d point(vertex[0], vertex[1], vertex[2]);
                bb.add(point);
            }
            */
           //only lock at centroid
           auto vertex = geom.center();
           Vec3d point(vertex[0], vertex[1], vertex[2]);
           bb.add(point);
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
        auto normal = ddm::normalOfElement(element);
        cell_normals_[elemIdx] = normal;

    }
    std::cout << "For Fracture : " << this->name() << " : " << tri_divide << " triangles should be devided"
              << std::endl;
    std::cout << "For Fracture : " << this->name() << " : " << tri_outside << " triangles outside" << std::endl;
    std::cout << "Total triangles: " << grid_->leafGridView().size(0) << std::endl;
    auto it = std::find(reservoir_cells_.begin(),reservoir_cells_.end(),-1);
    auto extended_fractures = prm_.get<bool>("extended_fractures");
    if(!(it == reservoir_cells_.end()) && !(extended_fractures) ){
        std::cout << "Remove fracture outside of model" << std::endl;
        // remove fracture outside of model
        this->removeCells();
        this->updateReservoirCells(cellSearchTree, grid3D);
    }
}


void Fracture::setFractureGrid(std::unique_ptr<Fracture::Grid>& gptr)
{
  grid_ = std::move(gptr);
  layers_ = 0;
  nlinear_ = 0;

  size_t nc = grid_->leafGridView().size(0);
  fracture_pressure_.resize(nc); fracture_pressure_ = 1e5;
  // reservoir_pressure_.resize(nc, 1.0e5); @@ will be done in 'updateReservoirProperties' below
  updateReservoirProperties();
  this->assembleFracture();
}

void
Fracture::updateReservoirProperties()
{
    // updater for standalone test
    double perm = prm_.get<double>("reservoir.perm");
    double dist = prm_.get<double>("reservoir.dist");
    double mobility = prm_.get<double>("reservoir.mobility");
    size_t nc = grid_->leafGridView().size(0);
    //assert(reservoir_cells_.size() == nc);
    reservoir_perm_.resize(nc, perm);
    reservoir_dist_.resize(nc, dist);
    reservoir_mobility_.resize(nc, 1000);
    reservoir_pressure_.resize(nc, 100.0e5);
    //reservoir_stress_.resize(nc);
    //for (size_t i = 0; i != nc; ++i) reservoir_stress_[i] = Dune::FieldVector<double, 6> {0, 0, 0, 0, 0, 0};
    //cell_normals_.resize(grid_->leafGridView().size(0)); // @@
    nu_ = 0.25;
    E_ = 1e9;
    this->initFractureWidth();
}
void
Fracture::solve()
{
    std::string method = prm_.get<std::string>("solver.method");
    if(method == "simple"){
        this->solveFractureWidth();
        this->solvePressure();
    }else if(method == "only_pressure"){
        this->solvePressure();
    }else if(method == "only_width"){
        this->solveFractureWidth();
    }else if(method == "iterative"){
        int max_it = prm_.get<int>("max_iter");
        int it=0;
         bool changed = true;
        while(changed && (it < max_it)){
            auto fracture_width = fracture_width_;
            auto fracture_pressure = fracture_pressure_;
            this->solveFractureWidth();
            // grow fracture
            this->solvePressure();
            it +=1;
            double tol = prm_.get<double>("solver.max_change");
            double max_change=0;
            for(int i=0;fracture_width_.size(); ++i){
                double diff_width = fracture_width_[i] - fracture_width[i];
                double diff_press = fracture_pressure_[i] - fracture_pressure[i];
                max_change = std::max(max_change,diff_width/1e-2);
                max_change = std::max(max_change,diff_press/1e5);
            }
            changed = (max_change < tol);
        }

    } else if (method == "if") {
      // ensure system matrices exist 
      if (!pressure_matrix_) initPressureMatrix();
      assembleFracture();

      // iterate full nonlinear system until convergence
      fracture_width_ = 1e-2;   // Ensure not completely closed
      fracture_pressure_ = 0.0;
      const int max_iter = 100; // @@ move elsewhere
      const double tol = 1e-5;
      int iter = 0;

      while (!fullSystemIteration(tol) && iter++ < max_iter) {};

      // report on result
      if (iter >= max_iter)
        std::cout << "WARNING: Did not converge in " << max_iter << " iterations." << std::endl;
      else
        std::cout << "System converged in " << iter << " iterations." << std::endl;

      
      // const auto ctrl = prm_.get_child("control");
      // const std::string ctrl_type(ctrl.get<std::string>("type"));

      // const std::vector<size_t> ratecells = (ctrl_type == "rate") ?
      //   std::vector<size_t>(well_source_.begin(), well_source_.end()) : std::vector<size_t>();
      // const double rate = ratecells.empty() ? 0.0 : ctrl.get<double>("rate") / ratecells.size();

      // const std::vector<size_t> bhpcells = (ctrl_type == "rate") ?
      //   std::vector<size_t>() : std::vector<size_t>(well_source_.begin(), well_source_.end());
      // const double bhp = (ctrl_type == "pressure") ? ctrl.get<double>("pressure") : perf_pressure_;

      // // fracture_pressure_ and fracture_width_ are the two unknowns, which will
      // // be computed/updated inside the function
      // if (!pressure_matrix_) initPressureMatrix();

      // auto fmat = *fracture_matrix_;
      // fmat *= -1;
      // solve_fully_coupled(fracture_pressure_, fracture_width_, // the two unknowns
      //                     *pressure_matrix_, fmat,
      //                     htrans_, rate, ratecells, bhp, bhpcells, leakof_);
    }else{
        OPM_THROW(std::runtime_error,"Unknowns solution method");
    }
}
void
Fracture::setSource()
{
    if (rhs_pressure_.size() == 0) {
        size_t nc = grid_->leafGridView().size(0);
        rhs_pressure_.resize(nc);
    }
    rhs_pressure_ = 0;

    for(size_t i=0; i< reservoir_pressure_.size(); ++i){
        rhs_pressure_[i] += leakof_[i]*reservoir_pressure_[i];
    }
    auto control=prm_.get_child("control");
    std::string control_type = control.get<std::string>("type");
    if (control_type == "rate") {
        double scale = well_source_.size();
        double rate = control.get<double>("rate");
        for (auto cell : well_source_) {
            rhs_pressure_[cell] += rate / scale;
        }
    } else if (control_type == "pressure") {
        double pressure = control.get<double>("pressure");
        for (const auto& perfinj : perfinj_) {
            int cell = std::get<0>(perfinj);
            double value = std::get<1>(perfinj);
            rhs_pressure_[cell] += value * pressure;
        }
    } else if (control_type == "perf_pressure") {
        for (const auto& perfinj : perfinj_) {
            int cell = std::get<0>(perfinj);
            double value = std::get<1>(perfinj);
            rhs_pressure_[cell] += value * perf_pressure_;
        }
    } else {
        OPM_THROW(std::runtime_error, "Unknowns control");
    }
}

double Fracture::injectionPressure() const{
    std::string control_type = prm_.get<std::string>("control.type");
    if (control_type == "rate") {
        double bhp = 0.0;
        double scale = well_source_.size();
        // could have corrected for WI
        for (auto cell : well_source_) {
            bhp += fracture_pressure_[cell] / scale;
        }
        return bhp;
    // } else if (control_type == "bhp") {
    //     double bhp = prm_.get<double>("control.bhp");
    } else if (control_type == "pressure") {
        double bhp = prm_.get<double>("control.pressure");
        return bhp;
    } else if (control_type == "perf_pressure") {
        double bhp = perf_pressure_;
        return bhp;
    } else {
        OPM_THROW(std::runtime_error, "Unknowns control");
    }
    return 0.0;
}

std::vector<double> Fracture::leakOfRate() const{
    if(leakof_.size()==0){
        std::vector<double> rate;
        return rate;
    }
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    std::vector<double> leakofrate(grid_->leafGridView().size(0),0);
    for (auto& element : Dune::elements(grid_->leafGridView())) {
        int eIdx = mapper.index(element);
        auto geom = element.geometry();
        double dp = (fracture_pressure_[eIdx]-reservoir_pressure_[eIdx]);
        double q = leakof_[eIdx]*dp;
        leakofrate[eIdx] = q/geom.volume();
    }
    return leakofrate;
}

std::vector<std::tuple<int,double,double>> Fracture::wellIndices() const{
    // find unique reservoir cells
    if(leakof_.size() == 0){
        // if pressure is not assembled return empty
        std::vector<std::tuple<int,double, double>> wellind;
        return wellind;
    }
    std::vector<int> res_cells = reservoir_cells_;
    std::sort(res_cells.begin(),res_cells.end());
    auto last = std::unique(res_cells.begin(),res_cells.end());
    res_cells.erase(last, res_cells.end());
    std::vector<double> q_cells(res_cells.size(),0.0);
    std::vector<double> p_cells(res_cells.size(),0.0);
    double q_prev = 0;
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto& element : Dune::elements(grid_->leafGridView())) {
        int eIdx = mapper.index(element);
        double dp = (fracture_pressure_[eIdx]-reservoir_pressure_[eIdx]);
        double q = leakof_[eIdx]*dp;
        int res_cell = reservoir_cells_[eIdx];
        // just to search
        auto it = std::find(res_cells.begin(), res_cells.end(), res_cell);
        int ind_wellIdx =  it - res_cells.begin();
        assert(!(it == res_cells.end()));
        if( q_prev*q < 0 ){
            OPM_THROW(std::runtime_error,"Cross flow in fracture ??");
        }
        q_cells[ind_wellIdx] += q;
        p_cells[ind_wellIdx] = reservoir_pressure_[eIdx];// is set multiple times
    }
    std::vector<std::tuple<int,double, double>> wellIndices(res_cells.size());
    double inj_press = injectionPressure();
    for(size_t i=0; i < res_cells.size(); ++i){
        double dp = inj_press - p_cells[i];
        // simplest possible approach
        // assumes leakof is assumed to be calculated with reservoir cell as reference
        double well_index = q_cells[i]/dp;
        wellIndices[i] = {res_cells[i],well_index, origo_[3]};
        //wellIndices[i] = std::tuple<int,double>({res_cells[i],well_index});
    }
    return wellIndices;
}

void
Fracture::writePressureSystem() const{
    if(prm_.get<bool>("write_pressure_system")){
        Dune::storeMatrixMarket(*pressure_matrix_, "pressure_matrix");
        Dune::storeMatrixMarket(rhs_pressure_, "pressure_rhs");
    }
}

void
Fracture::writeFractureSystem() const{
    if(prm_.get<bool>("write_fracture_system")){
        //Dune::storeMatrixMarket(*fracture_matrix_, "fracture_matrix");
        Dune::storeMatrixMarket(rhs_width_, "rhs_width");
    }
}

void
Fracture::solvePressure() {
    size_t nc = grid_->leafGridView().size(0);
    fracture_pressure_.resize(nc);
    fracture_pressure_ = 1e5;
    if(!pressure_matrix_){
        this->initPressureMatrix();
    }
    this->assemblePressure();
    this->setSource(); // probably include reservoir pressure
    this->writePressureSystem();
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

// based on the (given) values of fracture_pressure__, compute rhs_width_ and
// fracture_width_
void Fracture::solveFractureWidth()
{
    this->assembleFracture();
    fracture_matrix_->solve(fracture_width_,rhs_width_);
    double max_width = prm_.get<double>("solver.max_width");
    double min_width = prm_.get<double>("solver.min_width");
    for(int i=0; i < fracture_width_.size(); ++i){
        assert(std::isfinite(fracture_width_[i]));
        if(fracture_width_[i]> max_width){
            std::cout << "Limit Fracture width" << std::endl;
            fracture_width_[i] = max_width;
        }
        if(fracture_width_[i] < min_width){
            std::cout << "Remove small Fracture width" << std::endl;
            fracture_width_[i] = min_width;
        }
        assert(std::isfinite(fracture_width_[i]));
    }
}

void Fracture::initFractureStates()
{
    this->initFractureWidth();
    this->initFracturePressureFromReservoir();
}

void
Fracture::initFractureWidth()
{
    size_t nc = grid_->leafGridView().size(0);
    fracture_width_.resize(nc);
    fracture_width_ = prm_.get<double>("config.initial_fracture_width");
    // ElementMapper elemMapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    // for (auto& element : elements(grid_->leafGridView())) {
    //     const auto elemIdx = elemMapper.index(element);
    //     auto geom = element.geometry();
    //     auto vec_origo = geom.center() - origo_;
    //     double dist_origo = vec_origo.two_norm();
    //     fracture_width_[elemIdx] *= dist_origo;
    // }
}
void
Fracture::initFracturePressureFromReservoir(){
    size_t nc = reservoir_cells_.size();
    fracture_pressure_.resize(nc);
    for(size_t i=0; i < nc; ++i){
        fracture_pressure_[i] = reservoir_pressure_[i];
    }
}
void
Fracture::initPressureMatrix()
{
    // size_t num_columns = 0;
    //  index of the neighbour
    //
    double fWI = prm_.get<double>("fractureWI");
    for(int cell : well_source_){
        perfinj_.push_back({cell,fWI});
    }

    size_t nc = grid_->leafGridView().size(0);
    leakof_.resize(nc,0.0);
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    for (auto& element : Dune::elements(grid_->leafGridView())) {
        int eIdx = mapper.index(element);
        auto geom = element.geometry();
        // iterator over all intersections
        for (auto& is : Dune::intersections(grid_->leafGridView(),element)) {

            if (!is.boundary()) {
                auto eCenter = geom.center();
                int nIdx = mapper.index(is.outside());
                if (eIdx < nIdx) {
                    // calculate distance between the midpoints
                    auto nCenter = is.outside().geometry().center();
                    auto isCenter = is.geometry().center();
                    auto d_inside = eCenter - isCenter;
                    auto d_outside = nCenter - isCenter;
                    // probably should use projected distenace
                    auto igeom = is.geometry();
                    double area = igeom.volume();
                    double h1 = area/d_inside.two_norm();
                    double h2 = area/d_outside.two_norm();

                    {
                        Htrans matel(nIdx, eIdx, h1, h2);
                        htrans_.push_back(matel);
                    }
                }
            }
        }
        {
            //auto normal = ddm::normalOfElement(element);
            //auto permmatrix =  reservoir_perm_[eIdx];
            //auto pn = permmatrix.mv(normal);
            //double permval  = pn.dot(normal);

            // keap reservoir perm as n'K'n
            double value = reservoir_mobility_[eIdx]*(reservoir_perm_[eIdx] * geom.volume()) ;
            value /= reservoir_dist_[eIdx];
            leakof_[eIdx] = value;
        }
    }

    // not need if build mode is implicit
    // std::sort(transes.begin(),transes.end(),sortcsr);
    //  build matrix
    //if (pressure_matrix_.bu == 0){
    //size_t nc = grid_->leafGridView().size(0);
    pressure_matrix_ = std::make_unique<Matrix>(nc,nc, 4, 0.4, Matrix::implicit);
        auto& matrix = *pressure_matrix_;
        //matrix.setBuildMode(Matrix::implicit);
        // map from dof=3*nodes at a cell (ca 3*3*3) to cell
        //matrix.setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
        //size_t nc = grid_->leafGridView().size(0);
        //matrix.setSize(nc, nc);
        for (auto matel : htrans_) {
            size_t i = std::get<0>(matel);
            size_t j = std::get<1>(matel);
            double zero_entry = 0.0; //1e-11;
            matrix.entry(i, j) = zero_entry; //0;
            matrix.entry(j, i) = zero_entry; //0;
            matrix.entry(j, j) = zero_entry; //0;
            matrix.entry(i, i) = zero_entry; //0;
        }
        matrix.compress();
    //}
}
void
Fracture::assemblePressure()
{
    auto& matrix = *pressure_matrix_;
    matrix = 0.0;
    double mobility=1e4; // @@ 1.0

    for (auto matel : htrans_) {
        size_t i = std::get<0>(matel);
        size_t j = std::get<1>(matel);
        double t1 = std::get<2>(matel);
        double t2 = std::get<3>(matel);
        double h1 = fracture_width_[i];
        double h2 = fracture_width_[j];
        // harmonic mean of surface flow
        double value = 12. / (h1 * h1 * h1 * t1) + 12. / (h2 * h2 * h2 * t2);

        value = 1 / value;
        value *= mobility;
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
    auto control = prm_.get_child("control");
    std::string control_type = control.get<std::string>("type");
    for (size_t i = 0; i < leakof_.size(); ++i) {
            // matrix.entry(i, i) += leakof_[i];
            matrix[i][i] += leakof_[i];
    }
    if (control_type == "rate") {
        // no extra tings in matrix
    } else if (control_type == "pressure" || "perf_pressure") {
        for (const auto& perfinj : perfinj_) {
            int cell = std::get<0>(perfinj);
            double value = std::get<1>(perfinj);
            matrix[cell][cell] += value;
        }
    }else{
        OPM_THROW(std::runtime_error,"Unknown control of injection into Fracture");
    }
}

void  Fracture::assembleFracture(){
    size_t nc = grid_->leafGridView().size(0);
    //rhs_width_.resize(nc);
    rhs_width_ = fracture_pressure_;
    ElementMapper mapper(grid_->leafGridView(), Dune::mcmgElementLayout());
    if (reservoir_stress_.size() > 0) {
      for (auto elem : elements(grid_->leafGridView())) {
        int idx = mapper.index(elem);
        //auto geom = elem.geometry();
        rhs_width_[idx] -= reservoir_pressure_[idx];//*geom.volume();
        // using compressible stress ??
        rhs_width_[idx] -= ddm::tractionSymTensor(reservoir_stress_[idx],cell_normals_[idx]);
        // temperature
      }
    }
    // Do we need to ad thermal "forces"

    if(!fracture_matrix_){
        fracture_matrix_ = std::make_unique<DynamicMatrix>();
    }
    fracture_matrix_->resize(nc,nc);
    *fracture_matrix_ = 0.0;
    ddm::assembleMatrix(*fracture_matrix_,E_, nu_,*grid_);
}

void Fracture::printPressureMatrix() const // debug purposes
{
  Dune::printSparseMatrix(std::cout, *pressure_matrix_, "matname", "linameo");
}

void Fracture::printMechMatrix() const // debug purposes
{
  Dune::printmatrix(std::cout, *fracture_matrix_, "matname", "linameo");
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
                                             const Dune::CpGrid& grid3D);
template void
Fracture::updateReservoirCells<Dune::PolyhedralGrid<3,3,double>>(const external::cvf::ref<external::cvf::BoundingBoxTree>& cellSearchTree,
                                             const Dune::PolyhedralGrid<3,3,double>& grid3D);
} // namespace Opm
