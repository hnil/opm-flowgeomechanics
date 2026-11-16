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

#include "config.h"

#include "CgalHelper.hpp"
#include <CGAL/version.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <variant>
#include <optional>
template <typename T> struct typeOf;
std::vector<Point_3> intersect_hex_with_plane(const Mesh& mesh, const Plane_3& plane)
{
    std::vector<Point_3> pts;

    // Compute intersection of each mesh edge with plane
    for (auto e : mesh.edges())
    {
        auto h = mesh.halfedge(e);
        Point_3 p = mesh.point(mesh.source(h));
        Point_3 q = mesh.point(mesh.target(h));

        Segment_3 seg(p, q);
        auto result = CGAL::intersection(seg, plane);


        if (result)
        {
          // CGAL 6.0 replaced Boost optional/variant with std equivalents.
#if CGAL_VERSION_NR >= 1060000000
      const auto& vresult = result.value();
	  if (std::holds_alternative<Point_3>(vresult))
	    pts.push_back(std::get<Point_3>(vresult));
	  else if (std::holds_alternative<Segment_3>(vresult)) {
	    // An edge lying IN the plane intersects as a segment, not a point.
	    // Dropping it loses real cross-section and can leave fewer than three
	    // points for a cut that does exist - the 44-46 degree seeding abort.
	    const auto& s = std::get<Segment_3>(vresult);
	    pts.push_back(s.source());
	    pts.push_back(s.target());
	  }
#else //old version
	  if (const Point_3* ip = boost::get<Point_3>(&*result))
	         pts.push_back(*ip);
	  else if (const Segment_3* is = boost::get<Segment_3>(&*result)) {
	    // see the comment in the branch above: in-plane edge
	    pts.push_back(is->source());
	    pts.push_back(is->target());
	  }
#endif
        }
    }

    // remove duplicates
    std::sort(pts.begin(), pts.end(), [](auto& a, auto& b){
        if (a.x() != b.x()) return a.x() < b.x();
        if (a.y() != b.y()) return a.y() < b.y();
        return a.z() < b.z();
    });

    pts.erase(std::unique(pts.begin(), pts.end(), [](auto& a, auto& b){
        return CGAL::squared_distance(a,b) < 1e-14;
    }), pts.end());


    return pts;
}

Mesh make_hex(const std::vector<std::array<double,3>>& hex)
{
  // order should probably be mor correct to be used as mesh
    Mesh m;
    std::vector<Mesh::Vertex_index> v;
    for(const auto& point: hex){
      v.push_back(m.add_vertex(Point_3(point[0],point[1],point[2])));
    }

    auto add_quad = [&](auto a, auto b, auto c, auto d){
        m.add_face(a,b,c);
        m.add_face(a,c,d);
    };

    add_quad(v[0],v[1],v[3],v[2]);// top
    add_quad(v[4],v[5],v[7],v[6]);// bottom
    add_quad(v[0],v[1],v[5],v[4]);//fornt
    add_quad(v[1],v[3],v[7],v[5]);//right
    add_quad(v[3],v[2],v[6],v[7]);//back
    add_quad(v[2],v[0],v[4],v[6]);//left

    return m;

    return m;
}

double area_of_intersection(const std::vector<std::array<double,3>>& hex_org,
                            const std::vector<std::array<double,3>>& tri_org){
    Mesh hex = make_hex(hex_org);
    assert(tri_org.size() == 3);
    std::vector<Point_3> t;
    for(const auto& p: tri_org){
      t.push_back(Point_3(p[0],p[1],p[2]));
    }
    Plane_3 plane(t[0],t[1],t[2]); // z = 0.5
    // intersect of plane and hex
    auto poly = intersect_hex_with_plane(hex, plane);
    // The plane can miss the cell entirely, or meet it only at a corner or
    // along one edge. Fewer than three distinct points is no cross-section -
    // and indexing poly[0] on an empty vector is what aborted every deck
    // seeded at 44-46 degrees to the grid.
    if (poly.size() < 3) {
        return 0.0;
    }
    // make basis
    PlaneBasis B = make_basis(plane, poly[0]);
    // project polygon to 2D
    Polygon_2 P;
    for (auto& p : poly) P.push_back(project(p, B));
    // project triangle
    Polygon_2 T;
    T.push_back(project(t[0],B));
    T.push_back(project(t[1],B));
    T.push_back(project(t[2],B));
    // order points in 2D may be double problems
    sortPoints(P);
    sortPoints(T);
    // do intersection
    typedef CGAL::Polygon_with_holes_2<K_exact> Polygon_with_holes_2;
    std::vector<Polygon_with_holes_2> intersection_polygons;
    CGAL::intersection(P,T, std::back_inserter(intersection_polygons));
    if(intersection_polygons.size()==0){
      return 0.0;
    }
    double area = 0.0;
    for (auto& pp : intersection_polygons) {
        assert(pp.number_of_holes() == 0);
        auto outer = pp.outer_boundary();
        area += to_double(outer.area());
    }
    return area;
}

double area_of_intersection_fast(const std::vector<std::array<double,3>>& hex_org,
                                 const std::vector<std::array<double,3>>& tri_org)
{
    using V3 = std::array<double, 3>;
    const auto sub = [](const V3& a, const V3& b) { return V3{a[0] - b[0], a[1] - b[1], a[2] - b[2]}; };
    // local origin: field coordinates are large against cell sizes
    std::vector<V3> hex(hex_org.size()), tri(tri_org.size());
    for (size_t k = 0; k < hex.size(); ++k) hex[k] = sub(hex_org[k], tri_org[0]);
    for (size_t k = 0; k < tri.size(); ++k) tri[k] = sub(tri_org[k], tri_org[0]);
    const auto dot = [](const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
    const auto cross = [](const V3& a, const V3& b) {
        return V3{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    };
    assert(hex.size() == 8 && tri.size() == 3);
    V3 n = cross(sub(tri[1], tri[0]), sub(tri[2], tri[0]));
    const double nn = std::sqrt(dot(n, n));
    if (!(nn > 0.0))
        return 0.0;
    for (auto& c : n) c /= nn;
    const double d0 = dot(n, tri[0]);

    // the edges of make_hex: each quad (a,b,c,d) split into (a,b,c), (a,c,d)
    static constexpr int quads[6][4] = {{0,1,3,2},{4,5,7,6},{0,1,5,4},{1,3,7,5},{3,2,6,7},{2,0,4,6}};
    std::array<std::array<int, 2>, 30> edges{};
    int ne = 0;
    const auto add_edge = [&](int a, int b) {
        if (a > b) std::swap(a, b);
        for (int k = 0; k < ne; ++k)
            if (edges[k][0] == a && edges[k][1] == b) return;
        edges[ne++] = {a, b};
    };
    for (const auto& q : quads) {
        add_edge(q[0], q[1]); add_edge(q[1], q[2]); add_edge(q[2], q[0]);
        add_edge(q[2], q[3]); add_edge(q[3], q[0]);
    }
    std::array<double, 8> sd;
    for (int k = 0; k < 8; ++k) sd[k] = dot(n, hex[k]) - d0;
    std::vector<V3> pts;
    for (int k = 0; k < ne; ++k) {
        const int a = edges[k][0], b = edges[k][1];
        const double sa = sd[a], sb = sd[b];
        if (sa == 0.0 && sb == 0.0) { pts.push_back(hex[a]); pts.push_back(hex[b]); }
        else if (sa == 0.0) pts.push_back(hex[a]);
        else if (sb == 0.0) pts.push_back(hex[b]);
        else if ((sa < 0.0) != (sb < 0.0)) {
            const double t = sa / (sa - sb);
            const V3 e = sub(hex[b], hex[a]);
            pts.push_back(V3{hex[a][0] + t * e[0], hex[a][1] + t * e[1], hex[a][2] + t * e[2]});
        }
    }
    std::sort(pts.begin(), pts.end());
    pts.erase(std::unique(pts.begin(), pts.end(), [&](const V3& a, const V3& b) {
        const V3 r = sub(a, b);
        return dot(r, r) < 1e-14;
    }), pts.end());
    if (pts.size() < 3)
        return 0.0;

    V3 u = cross(n, V3{1, 0, 0});
    if (dot(u, u) < 1e-12) u = cross(n, V3{0, 1, 0});
    const double un = std::sqrt(dot(u, u));
    for (auto& c : u) c /= un;
    const V3 v = cross(n, u);
    using P2 = std::array<double, 2>;
    const auto proj = [&](const V3& p) { const V3 r = sub(p, pts[0]); return P2{dot(r, u), dot(r, v)}; };
    std::vector<P2> poly;
    double cx = 0, cy = 0;
    for (const auto& p : pts) { poly.push_back(proj(p)); cx += poly.back()[0]; cy += poly.back()[1]; }
    cx /= poly.size(); cy /= poly.size();
    std::sort(poly.begin(), poly.end(), [&](const P2& a, const P2& b) {
        return std::atan2(a[1] - cy, a[0] - cx) < std::atan2(b[1] - cy, b[0] - cx);
    });
    std::array<P2, 3> T{proj(tri[0]), proj(tri[1]), proj(tri[2])};
    const auto orient = [](const P2& a, const P2& b, const P2& c) {
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    };
    if (orient(T[0], T[1], T[2]) < 0.0) std::swap(T[1], T[2]);

    // Sutherland-Hodgman: clip the cross-section by the (convex) triangle
    std::vector<P2> out;
    for (int e = 0; e < 3 && !poly.empty(); ++e) {
        const P2& a = T[e];
        const P2& b = T[(e + 1) % 3];
        out.clear();
        for (size_t k = 0; k < poly.size(); ++k) {
            const P2& p = poly[k];
            const P2& q = poly[(k + 1) % poly.size()];
            const double sp = orient(a, b, p), sq = orient(a, b, q);
            if (sp >= 0.0) out.push_back(p);
            if ((sp >= 0.0) != (sq >= 0.0)) {
                const double t = sp / (sp - sq);
                out.push_back(P2{p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])});
            }
        }
        poly.swap(out);
    }
    double area = 0.0;
    for (size_t k = 0; k < poly.size(); ++k) {
        const P2& p = poly[k];
        const P2& q = poly[(k + 1) % poly.size()];
        area += p[0] * q[1] - q[0] * p[1];
    }
    return 0.5 * std::abs(area);
}

// -----------------------------------------------------------
// Build a HEX manually as a Surface_mesh
// -----------------------------------------------------------

PlaneBasis make_basis(const Plane_3& plane, const Point_3& origin)
{
    Vector_3 n = plane.orthogonal_vector();
    n = n / sqrt(n.squared_length());

    Vector_3 u = CGAL::cross_product(n, Vector_3(1,0,0));
    if (u.squared_length() < 1e-12)
        u = CGAL::cross_product(n, Vector_3(0,1,0));

    u = u / sqrt(u.squared_length());
    Vector_3 v = CGAL::cross_product(n, u);

    return {origin, u, v};
}

Point_2 project(const Point_3& p, const PlaneBasis& B)
{
    Vector_3 d = p - B.origin;
    return Point_2(d * B.u, d * B.v);
}

void sortPoints(Polygon_2& pts){
    double cx=0, cy=0;
    for (auto& p : pts) {
      cx += CGAL::to_double(p.x()); cy += CGAL::to_double(p.y());
    }
    Point_2 C(cx / pts.size(), cy / pts.size());
    auto proj = [&](const Point_2& p){
        Vector_2 d = p - C;
        double y = CGAL::to_double(d.y());
        double x = CGAL::to_double(d.x());
        return std::atan2(y, x);
    };
    // should probably remove identical
    std::sort(pts.begin(), pts.end(),
              [&](const Point_2& a, const Point_2& b){ return proj(a) < proj(b); });
}
