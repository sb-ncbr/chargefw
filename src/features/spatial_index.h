#pragma once

#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>

#include <nanoflann.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace chargefw::features::detail {

class ConformerPointCloud {
  public:
    explicit ConformerPointCloud(const ConformerFeatures& geometry) : geometry_{&geometry} {
        if (geometry.first_nonfinite_atom_index().has_value()) {
            throw std::invalid_argument{"spatial index requires finite coordinates"};
        }
    }

    [[nodiscard]] auto kdtree_get_point_count() const -> std::size_t {
        return geometry_->molecule().atom_count();
    }

    [[nodiscard]] auto kdtree_get_pt(const std::size_t point_index,
                                     const std::size_t dimension) const -> double {
        const auto& position = geometry_->position(point_index);
        switch (dimension) {
        case 0:
            return position.x;
        case 1:
            return position.y;
        default:
            return position.z;
        }
    }

    template <class BoundingBox> auto kdtree_get_bbox(BoundingBox&) const -> bool {
        return false;
    }

  private:
    const ConformerFeatures* geometry_;
};

using KdTree =
    nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, ConformerPointCloud>,
                                        ConformerPointCloud, 3, std::size_t>;

class SourceIndexRadiusResultSet {
  public:
    using DistanceType = double;
    using IndexType = std::size_t;

    SourceIndexRadiusResultSet(const double radius_squared,
                               std::vector<std::size_t>& neighbor_indices)
        : radius_squared_{radius_squared}, neighbor_indices_{&neighbor_indices} {
        neighbor_indices_->clear();
    }

    [[nodiscard]] auto size() const noexcept -> std::size_t {
        return neighbor_indices_->size();
    }
    [[nodiscard]] auto full() const noexcept -> bool {
        return true;
    }
    [[nodiscard]] auto worstDist() const noexcept -> double {
        return radius_squared_;
    }

    auto addPoint(const double squared_distance, const std::size_t point_index) -> bool {
        if (squared_distance < radius_squared_) {
            neighbor_indices_->push_back(point_index);
        }
        return true;
    }

    auto sort() -> void {
        std::ranges::sort(*neighbor_indices_);
    }

  private:
    double radius_squared_;
    std::vector<std::size_t>* neighbor_indices_;
};

// A k-d tree over the atoms of one conformer, which must have finite coordinates. It refers to
// the conformer features, which must outlive it.
class ConformerSpatialIndex {
  public:
    explicit ConformerSpatialIndex(const ConformerFeatures& geometry)
        : points_{geometry}, tree_{3, points_} {}

    // Replaces the contents of neighbor_indices with the atoms within radius of center, including
    // atoms at exactly that distance, in no particular order.
    auto neighbor_indices_within(const core::Position& center, const double radius,
                                 std::vector<std::size_t>& neighbor_indices) const -> void {
        const std::array<double, 3> query{center.x, center.y, center.z};
        const auto inclusive_radius_squared =
            std::nextafter(radius * radius, std::numeric_limits<double>::infinity());

        SourceIndexRadiusResultSet result_set{inclusive_radius_squared, neighbor_indices};
        static_cast<void>(tree_.radiusSearchCustomCallback(
            query.data(), result_set, nanoflann::SearchParameters{0.0F, false}));
    }

  private:
    ConformerPointCloud points_;
    KdTree tree_;
};

} // namespace chargefw::features::detail
