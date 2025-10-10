#pragma once

#include <csf_pcl/cloth_simulation_filter.h>
#include <pcl/common/common.h>

#include <limits>
#include <utility>
#include <vector>

namespace csf_pcl
{
template <typename PointT>
ClothSimulationFilter<PointT>::Cloth::Cloth(
  const Eigen::Vector2f & min_point, const Eigen::Vector2f & max_point, const float initial_z,
  const float resolution, const float margin, const size_t rigidness, const float gravity,
  const float time_step)
: resolution_(resolution),
  margin_(margin),
  rigidness_(rigidness),
  gravity_(gravity),
  time_step_(time_step)
{
  origin_ = min_point - Eigen::Vector2f::Constant(margin_);
  size_ = ((max_point - min_point + Eigen::Vector2f::Constant(2 * margin_)) / resolution_)
            .cast<Eigen::Index>();

  current_height_map_ = Eigen::ArrayXXf::Constant(size_.x(), size_.y(), initial_z);
  previous_height_map_ = current_height_map_;
  movability_map_ =
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(size_.x(), size_.y(), true);
}

template <typename PointT>
void ClothSimulationFilter<PointT>::Cloth::step()
{
  const auto new_height_map = current_height_map_ + movability_map_.cast<float>() *
                                                      (current_height_map_ - previous_height_map_ +
                                                       gravity_ * std::pow(time_step_, 2));
  previous_height_map_ = current_height_map_;
  current_height_map_ = new_height_map;
}

template <typename PointT>
void ClothSimulationFilter<PointT>::Cloth::checkIntersection(
  const Eigen::ArrayXXf & intersection_height_map)
{
  movability_map_ = (current_height_map_ < intersection_height_map);
  current_height_map_ = movability_map_.select(current_height_map_, intersection_height_map);
}

template <typename PointT>
void ClothSimulationFilter<PointT>::Cloth::applyConstraint()
{
  for (size_t i = 0; i < rigidness_; ++i) {
    Eigen::ArrayXXf displacement_map = Eigen::ArrayXXf::Zero(size_.x(), size_.y());

    displacement_map.block(1, 0, size_.x() - 1, size_.y()) +=
      current_height_map_.block(0, 0, size_.x() - 1, size_.y());
    displacement_map.block(0, 0, size_.x() - 1, size_.y()) +=
      current_height_map_.block(1, 0, size_.x() - 1, size_.y());
    displacement_map.block(0, 1, size_.x(), size_.y() - 1) +=
      current_height_map_.block(0, 0, size_.x(), size_.y() - 1);
    displacement_map.block(0, 0, size_.x(), size_.y() - 1) +=
      current_height_map_.block(0, 1, size_.x(), size_.y() - 1);
    displacement_map.row(0) += current_height_map_.row(0);
    displacement_map.row(size_.x() - 1) += current_height_map_.row(size_.x() - 1);
    displacement_map.col(0) += current_height_map_.col(0);
    displacement_map.col(size_.y() - 1) += current_height_map_.col(size_.y() - 1);
    displacement_map = (displacement_map / 4 - current_height_map_) / 2;

    current_height_map_ += displacement_map * movability_map_.cast<float>();
  }
}

template <typename PointT>
float ClothSimulationFilter<PointT>::Cloth::getCurrentHeight(
  const GridCoordinates & grid_coords) const
{
  return current_height_map_(grid_coords.x(), grid_coords.y());
}

template <typename PointT>
Eigen::Vector2f ClothSimulationFilter<PointT>::Cloth::fromGridCoordinates(
  const GridCoordinates & grid_coords) const
{
  return origin_ + grid_coords.template cast<float>() * resolution_;
}

template <typename PointT>
ClothSimulationFilter<PointT>::ClothSimulationFilter(const bool extract_removed_indices)
: pcl::FilterIndices<PointT>(extract_removed_indices)
{
  this->filter_name_ = "ClothSimulationFilter";
}

template <typename PointT>
void ClothSimulationFilter<PointT>::applyFilter(pcl::Indices & indices)
{
  Eigen::Vector4f input_min_point_, input_max_point_;
  pcl::getMinMax3D(*this->input_, input_min_point_, input_max_point_);

  Cloth cloth(
    input_min_point_.head<2>(), input_max_point_.head<2>(),
    input_min_point_.z() - cloth_initial_z_offset_, cloth_resolution_, cloth_margin_,
    cloth_rigidness_, gravity_, time_step_);

  const auto [intersection_height_map, input_grid_coords] = generateIntersectionHeightMap(cloth);

  for (size_t i = 0; i < max_iterations_; ++i) {
    cloth.step();
    cloth.checkIntersection(intersection_height_map);
    cloth.applyConstraint();

    const auto max_height_variation =
      (cloth.currentHeightMap() - cloth.previousHeightMap()).abs().maxCoeff();
    if (max_height_variation < iteration_termination_threshold_) {
      break;
    }
  }

  indices.clear();
  indices.reserve(input_grid_coords.size());

  if (this->extract_removed_indices_) {
    this->removed_indices_->clear();
    this->removed_indices_->reserve(input_grid_coords.size());
  }

  for (size_t i = 0; i < input_grid_coords.size(); ++i) {
    if (!input_grid_coords[i]) {
      continue;
    }

    if (
      std::abs((*this->input_)[i].z - cloth.getCurrentHeight(*input_grid_coords[i])) <
      classification_threshold_) {
      if (!this->negative_) {
        indices.push_back(i);
      } else if (this->extract_removed_indices_) {
        this->removed_indices_->push_back(i);
      }
    } else {
      if (this->negative_) {
        indices.push_back(i);
      } else if (this->extract_removed_indices_) {
        this->removed_indices_->push_back(i);
      }
    }
  }

  indices.shrink_to_fit();
  if (this->extract_removed_indices_) {
    this->removed_indices_->shrink_to_fit();
  }
}

template <typename PointT>
std::pair<Eigen::ArrayXXf, std::vector<std::optional<Eigen::Vector<Eigen::Index, 2>>>>
ClothSimulationFilter<PointT>::generateIntersectionHeightMap(const Cloth & cloth) const
{
  Eigen::ArrayXXf intersection_height_map = Eigen::ArrayXXf::Constant(
    cloth.size().x(), cloth.size().y(), std::numeric_limits<float>::max());
  Eigen::ArrayXXf distance_from_cell_center_map = Eigen::ArrayXXf::Constant(
    cloth.size().x(), cloth.size().y(), std::numeric_limits<float>::max());

  std::vector<std::optional<GridCoordinates>> input_grid_coords;
  input_grid_coords.reserve(this->input_->size());

  for (const auto & point : *this->input_) {
    const auto grid_coords = cloth.toGridCoordinates(toEigen2f(point));
    input_grid_coords.push_back(grid_coords);
    if (!grid_coords) {
      continue;
    }

    const Eigen::Vector2f cell_center =
      cloth.fromGridCoordinates(*grid_coords) + Eigen::Vector2f::Constant(cloth_resolution_ / 2);
    const auto distance_from_cell_center = (toEigen2f(point) - cell_center).norm();
    if (
      distance_from_cell_center >
      distance_from_cell_center_map(grid_coords->x(), grid_coords->y())) {
      continue;
    }

    distance_from_cell_center_map(grid_coords->x(), grid_coords->y()) = distance_from_cell_center;
    intersection_height_map(grid_coords->x(), grid_coords->y()) = point.z;
  }

  return {intersection_height_map, input_grid_coords};
}
}  // namespace csf_pcl

#define PCL_INSTANTIATE_ClothSimulationFilter(T) \
  template class PCL_EXPORTS csf_pcl::ClothSimulationFilter<T>;
