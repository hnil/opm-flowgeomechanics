/*
  Copyright 2025 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/draw_polygon_2.h>
#include <CGAL/draw_polygon_with_holes_2.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <cmath>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K_inexact;
typedef CGAL::Exact_predicates_exact_constructions_kernel K_exact;
typedef K_inexact::Point_3 Point_3;
typedef K_exact::Point_2 Point_2;
typedef K_inexact::Vector_3 Vector_3;
typedef K_exact::Vector_2 Vector_2;
typedef K_inexact::Segment_3 Segment_3;
typedef K_inexact::Plane_3 Plane_3;

typedef CGAL::Surface_mesh<Point_3> Mesh;
typedef CGAL::Polygon_2<K_exact> Polygon_2;
typedef K_exact::Point_2 Point_2;
// -----------------------------------------------------------
// Compute hex ∩ plane using segment-plane intersections
// -----------------------------------------------------------
struct PlaneBasis {
    Point_3 origin;
    Vector_3 u, v;
};

std::vector<Point_3> intersect_hex_with_plane(const Mesh& mesh, const Plane_3& plane);
PlaneBasis make_basis(const Plane_3& plane, const Point_3& origin);
Point_2 project(const Point_3& p, const PlaneBasis& B);
void sortPoints(Polygon_2& pts);
Mesh make_hex(const std::vector<std::array<double,3>>& hex);
double area_of_intersection(const std::vector<std::array<double,3>>& hex_org,
                            const std::vector<std::array<double,3>>& tri_org);
