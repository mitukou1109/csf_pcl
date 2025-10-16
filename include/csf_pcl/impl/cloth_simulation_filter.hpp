#pragma once

#include <csf_pcl/cloth_simulation_filter.h>
#include <pcl/common/common.h>

#include <algorithm>
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
void ClothSimulationFilter<PointT>::Cloth::applyPostProcessing(
  const Eigen::ArrayXXf & intersection_height_map, const float slope_fitting_threshold)
{
  std::vector<GridCoordinates> movable_particles;

  for (Eigen::Index i = 0; i < size_.x(); ++i) {
    for (Eigen::Index j = 0; j < size_.y(); ++j) {
      if (!movability_map_(i, j)) {
        continue;
      }
      movable_particles.emplace_back(i, j);
    }
  }

  constexpr std::array<std::array<Eigen::Index, 2>, 4> neighbors = {
    {{-1, 0}, {0, -1}, {1, 0}, {0, 1}}};

  while (true) {
    std::vector<std::pair<GridCoordinates, std::vector<GridCoordinates>>>
      edge_particles_with_unmovable_neighbors;

    for (auto it = movable_particles.begin(); it != movable_particles.end();) {
      if (!movability_map_(it->x(), it->y())) {
        it = movable_particles.erase(it);
        continue;
      }

      bool is_on_edge = false;
      std::vector<GridCoordinates> unmovable_neighbors;

      for (size_t n = 0; n < neighbors.size(); ++n) {
        const auto ni = it->x() + neighbors[n][0];
        const auto nj = it->y() + neighbors[n][1];
        if (ni < 0 || ni >= size_.x() || nj < 0 || nj >= size_.y()) {
          continue;
        }
        if (!movability_map_(ni, nj)) {
          is_on_edge = true;
          unmovable_neighbors.emplace_back(ni, nj);
        }
      }

      if (is_on_edge) {
        edge_particles_with_unmovable_neighbors.emplace_back(*it, unmovable_neighbors);
        it = movable_particles.erase(it);
        continue;
      }
      ++it;
    }

    if (edge_particles_with_unmovable_neighbors.empty()) {
      break;
    }

    for (const auto & [edge_particle, unmovable_neighbors] :
         edge_particles_with_unmovable_neighbors) {
      applySlopeFitting(
        intersection_height_map, edge_particle, unmovable_neighbors, slope_fitting_threshold);
    }
  }
}

template <typename PointT>
void ClothSimulationFilter<PointT>::Cloth::applySlopeFitting(
  const Eigen::ArrayXXf & intersection_height_map, const GridCoordinates & edge_particle,
  const std::vector<GridCoordinates> & unmovable_neighbors, const float slope_fitting_threshold)
{
  for (const auto & neighbor : unmovable_neighbors) {
    const auto height_difference = (intersection_height_map(edge_particle.x(), edge_particle.y()) -
                                    current_height_map_(edge_particle.x(), edge_particle.y())) -
                                   (intersection_height_map(neighbor.x(), neighbor.y()) -
                                    current_height_map_(neighbor.x(), neighbor.y()));

    if (std::abs(height_difference) > slope_fitting_threshold) {
      continue;
    }

    current_height_map_(edge_particle.x(), edge_particle.y()) =
      intersection_height_map(edge_particle.x(), edge_particle.y());
    movability_map_(edge_particle.x(), edge_particle.y()) = false;
    break;
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

  if (enable_post_processing_) {
    cloth.applyPostProcessing(intersection_height_map, slope_fitting_threshold_);
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
  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> valid_map =
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic>::Constant(
      cloth.size().x(), cloth.size().y(), false);

  std::vector<std::optional<GridCoordinates>> input_grid_coords;
  input_grid_coords.reserve(this->input_->size());

  for (const auto & point : *this->input_) {
    const auto grid_coords = cloth.toGridCoordinates(toEigen2f(point));
    input_grid_coords.push_back(grid_coords);
    if (!grid_coords) {
      continue;
    }

    intersection_height_map(grid_coords->x(), grid_coords->y()) =
      std::min(intersection_height_map(grid_coords->x(), grid_coords->y()), point.z);
    valid_map(grid_coords->x(), grid_coords->y()) = true;
  }

  intersection_height_map = valid_map.select(intersection_height_map, 0);

  for (size_t i = 0; i < intersection_height_interpolation_max_iterations_; ++i) {
    Eigen::ArrayXXf padded =
      Eigen::ArrayXXf::Zero(intersection_height_map.rows() + 2, intersection_height_map.cols() + 2);
    padded.block(1, 1, intersection_height_map.rows(), intersection_height_map.cols()) =
      intersection_height_map;
    intersection_height_map = valid_map.select(
      intersection_height_map,
      (padded.block(0, 1, intersection_height_map.rows(), intersection_height_map.cols()) +
       padded.block(2, 1, intersection_height_map.rows(), intersection_height_map.cols()) +
       padded.block(1, 0, intersection_height_map.rows(), intersection_height_map.cols()) +
       padded.block(1, 2, intersection_height_map.rows(), intersection_height_map.cols())) /
        4);
    if (
      (intersection_height_map -
       padded.block(1, 1, intersection_height_map.rows(), intersection_height_map.cols()))
        .abs()
        .maxCoeff() < intersection_height_interpolation_termination_threshold_) {
      break;
    }
  }

  return {intersection_height_map, input_grid_coords};
}
}  // namespace csf_pcl

#define PCL_INSTANTIATE_ClothSimulationFilter(T) \
  template class PCL_EXPORTS csf_pcl::ClothSimulationFilter<T>;
