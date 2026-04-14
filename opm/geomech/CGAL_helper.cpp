#include "CGAL_helper.hpp"
#include <variant>
#include <optional>
// namespace Opm{
//   namespace CGAL{
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

	//typeOf<decltype(result)> typeOfResult;

        if (result)
        {
	  // old cgal
	  //if (const Point_3* ip = boost::get<Point_3>(&*result))
	  //       pts.push_back(*ip);
	  const auto& vresult = result.value();
	  if (std::holds_alternative<Point_3>(vresult))
	    pts.push_back(std::get<Point_3>(vresult));
	  //if (const Point_3* ip = std::get<Point_3>(*result))
	  //    pts.push_back(ip);
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

    // if (pts.size() < 3)
    //     return pts;

    // // Sort in a circular order
    // // centroid
    // K::FT cx=0, cy=0, cz=0;
    // for (auto& p : pts) {
    //     cx += p.x(); cy += p.y(); cz += p.z();
    // }
    // Point_3 C(cx / pts.size(), cy / pts.size(), cz / pts.size());

    // Vector_3 n = plane.orthogonal_vector();
    // n = n / sqrt(n.squared_length());

    // Vector_3 u = std::abs(n.x()) > 0.9 ? Vector_3(0,1,0) : Vector_3(1,0,0);
    // u = CGAL::cross_product(n, u);
    // u = u / sqrt(u.squared_length());
    // Vector_3 v = CGAL::cross_product(n, u);

    // auto proj = [&](const Point_3& p){
    //     Vector_3 d = p - C;
    //     return std::atan2(d * v, d * u);
    // };

    // std::sort(pts.begin(), pts.end(),
    //           [&](const Point_3& a, const Point_3& b){ return proj(a) < proj(b); });

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
        //std::cout << "Numholes " << pp.number_of_holes() << std::endl;
        auto outer = pp.outer_boundary();
        area += to_double(outer.area());
    }
    return area;
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
//  }
//}
