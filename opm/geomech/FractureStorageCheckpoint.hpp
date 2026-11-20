/* Copyright 2026. SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef OPM_FRACTURE_STORAGE_CHECKPOINT_HPP
#define OPM_FRACTURE_STORAGE_CHECKPOINT_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Opm {

// A geometric snapshot, independent of both Dune cell indices and the current
// mesh resolution. The rollback state remains untouched. Aperture is piecewise
// constant on checkpoint triangles, so integrating its overlap with a new cell
// transfers old fluid volume, rather than interpolating an intensive variable.
class FractureStorageCheckpoint {
public:
    using Point = std::array<double, 3>;
    using Triangle = std::array<Point, 3>;

    void capture(const std::vector<Triangle>& triangles,
                 const std::vector<double>& widths)
    {
        if (triangles.size() != widths.size())
            throw std::invalid_argument("Fracture storage checkpoint size mismatch");
        source_.clear();
        cached_triangles_.clear();
        cached_volumes_.clear();
        if (triangles.empty()) return;
        origin_ = triangles.front()[0];
        axis_ = subtract(triangles.front()[1], origin_);
        normalize(axis_);
        normal_ = cross(axis_, subtract(triangles.front()[2], origin_));
        normalize(normal_);
        second_axis_ = cross(normal_, axis_);
        scale_ = 1.0;
        for (const auto& triangle : triangles)
            for (const auto& p : triangle)
                scale_ = std::max(scale_, std::sqrt(dot(subtract(p, origin_), subtract(p, origin_))));
        for (std::size_t i = 0; i < triangles.size(); ++i) {
            if (!std::isfinite(widths[i]))
                throw std::invalid_argument("Non-finite fracture checkpoint aperture");
            auto t = project(triangles[i]);
            if (!(area(t) > 0.0))
                throw std::invalid_argument("Degenerate fracture checkpoint triangle");
            source_.push_back({t, bounds(t), widths[i]});
        }
        std::sort(source_.begin(), source_.end(), [](const auto& a, const auto& b) {
            return a.box[0] < b.box[0];
        });
    }

    const std::vector<double>& volumes(const std::vector<Triangle>& triangles) const
    {
        if (triangles == cached_triangles_ && cached_volumes_.size() == triangles.size())
            return cached_volumes_;
        std::vector<double> volumes(triangles.size(), 0.0);
        if (!source_.empty()) {
            for (std::size_t j = 0; j < triangles.size(); ++j) {
                const auto target = project(triangles[j]);
                const auto box = bounds(target);
                const auto end = std::upper_bound(source_.begin(), source_.end(), box[1],
                    [](double x, const auto& s) { return x < s.box[0]; });
                for (auto it = source_.begin(); it != end; ++it) {
                    if (it->box[1] < box[0] || it->box[3] < box[2] || box[3] < it->box[2])
                        continue;
                    volumes[j] += overlap(target, it->triangle) * it->width;
                }
            }
        }
        // Only mark a geometry as cached after a successful complete transfer.
        cached_triangles_ = triangles;
        cached_volumes_ = std::move(volumes);
        return cached_volumes_;
    }

private:
    using Point2 = std::array<double, 2>;
    using Triangle2 = std::array<Point2, 3>;
    struct Source { Triangle2 triangle; std::array<double, 4> box; double width; };
    std::vector<Source> source_;
    Point origin_{}, axis_{}, second_axis_{}, normal_{};
    double scale_ = 1.0;
    mutable std::vector<Triangle> cached_triangles_;
    mutable std::vector<double> cached_volumes_;

    static Point subtract(const Point& a, const Point& b)
    { return {a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
    static double dot(const Point& a, const Point& b)
    { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
    static Point cross(const Point& a, const Point& b)
    { return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]}; }
    static void normalize(Point& a) {
        const double length = std::sqrt(dot(a,a));
        if (!(length > 0.0) || !std::isfinite(length))
            throw std::invalid_argument("Invalid fracture storage plane");
        for (auto& x : a) x /= length;
    }
    static double side(const Point2& a, const Point2& b, const Point2& p)
    { return (b[0]-a[0])*(p[1]-a[1]) - (b[1]-a[1])*(p[0]-a[0]); }
    Triangle2 project(const Triangle& triangle) const {
        Triangle2 result;
        for (std::size_t i = 0; i < 3; ++i) {
            const auto p = subtract(triangle[i], origin_);
            const double length = std::sqrt(dot(p,p));
            if (!std::isfinite(length) || std::abs(dot(p,normal_)) > 1e-9*std::max(scale_, length))
                throw std::invalid_argument("Fracture storage requires coplanar checkpoint and current mesh");
            result[i] = {dot(p,axis_), dot(p,second_axis_)};
        }
        if (side(result[0], result[1], result[2]) < 0.0)
            std::swap(result[1], result[2]);
        return result;
    }
    static double area(const Triangle2& t)
    { return 0.5*std::abs(side(t[0],t[1],t[2])); }
    static std::array<double, 4> bounds(const Triangle2& t) {
        return {std::min({t[0][0],t[1][0],t[2][0]}), std::max({t[0][0],t[1][0],t[2][0]}),
                std::min({t[0][1],t[1][1],t[2][1]}), std::max({t[0][1],t[1][1],t[2][1]})};
    }
    static double overlap(const Triangle2& target, const Triangle2& source) {
        // Sutherland-Hodgman clipping of two convex triangles. Shared edges have
        // zero measure; no aperture floor or area renormalisation is introduced.
        std::vector<Point2> polygon(target.begin(), target.end());
        for (std::size_t edge = 0; edge < 3 && !polygon.empty(); ++edge) {
            const auto& a = source[edge]; const auto& b = source[(edge+1)%3];
            std::vector<Point2> clipped;
            auto previous = polygon.back(); double dp = side(a,b,previous);
            for (const auto& current : polygon) {
                const double dc = side(a,b,current);
                if ((dp >= 0.0) != (dc >= 0.0)) {
                    const double fraction = dp/(dp-dc);
                    clipped.push_back({previous[0]+fraction*(current[0]-previous[0]),
                                       previous[1]+fraction*(current[1]-previous[1])});
                }
                if (dc >= 0.0) clipped.push_back(current);
                previous = current; dp = dc;
            }
            polygon = std::move(clipped);
        }
        if (polygon.size() < 3) return 0.0;
        double twice_area = 0.0;
        for (std::size_t i = 1; i+1 < polygon.size(); ++i)
            twice_area += side(polygon[0],polygon[i],polygon[i+1]);
        return 0.5*std::abs(twice_area);
    }
};
} // namespace Opm
#endif
