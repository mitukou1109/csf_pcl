#pragma once

#include <pcl/filters/filter_indices.h>

#include <optional>
#include <utility>
#include <vector>

namespace csf_pcl
{
template <typename PointT>
class ClothSimulationFilter : public pcl::FilterIndices<PointT>
{
  using PointCloud = typename pcl::FilterIndices<PointT>::PointCloud;
  using GridCoordinates = Eigen::Vector<Eigen::Index, 2>;

public:
  explicit ClothSimulationFilter(const bool extract_removed_indices = false);

  void setClothResolution(const float cloth_resolution) { cloth_resolution_ = cloth_resolution; }
  void setClothMargin(const float cloth_margin) { cloth_margin_ = cloth_margin; }
  void setClothRigidness(const size_t cloth_rigidness) { cloth_rigidness_ = cloth_rigidness; }
  void setClothInitialZOffset(const float cloth_initial_z_offset)
  {
    cloth_initial_z_offset_ = cloth_initial_z_offset;
  }
  void setGravity(const float gravity) { gravity_ = gravity; }
  void setTimeStep(const float time_step) { time_step_ = time_step; }
  void setMaxIterations(const size_t max_iterations) { max_iterations_ = max_iterations; }
  void setIterationTerminationThreshold(const float iteration_termination_threshold)
  {
    iteration_termination_threshold_ = iteration_termination_threshold;
  }
  void enablePostProcessing(const bool enable) { enable_post_processing_ = enable; }
  void setSlopeFittingThreshold(const float slope_fitting_threshold)
  {
    slope_fitting_threshold_ = slope_fitting_threshold;
  }
  void setClassThreshold(const float classification_threshold)
  {
    classification_threshold_ = classification_threshold;
  }

  float getClothResolution() const { return cloth_resolution_; }
  float getClothMargin() const { return cloth_margin_; }
  size_t getClothRigidness() const { return cloth_rigidness_; }
  float getClothInitialZOffset() const { return cloth_initial_z_offset_; }
  float getGravity() const { return gravity_; }
  float getTimeStep() const { return time_step_; }
  size_t getMaxIterations() const { return max_iterations_; }
  float getIterationTerminationThreshold() const { return iteration_termination_threshold_; }
  bool isPostProcessingEnabled() const { return enable_post_processing_; }
  float getSlopeFittingThreshold() const { return slope_fitting_threshold_; }
  float getClassThreshold() const { return classification_threshold_; }

protected:
  void applyFilter(pcl::Indices & indices) override;

private:
  class Cloth
  {
  public:
    Cloth(
      const Eigen::Vector2f & min_point, const Eigen::Vector2f & max_point, const float initial_z,
      const float resolution, const float margin, const size_t rigidness, const float gravity,
      const float time_step);

    void step();

    void checkIntersection(const Eigen::ArrayXXf & intersection_height_map);

    void applyConstraint();

    void applyPostProcessing(
      const Eigen::ArrayXXf & intersection_height_map, const float slope_fitting_threshold);

    float getCurrentHeight(const GridCoordinates & grid_coords) const;

    template <typename T>
    std::optional<GridCoordinates> toGridCoordinates(const T & point_2d) const
    {
      const auto grid_coords = ((point_2d - origin_) / resolution_).template cast<Eigen::Index>();
      if ((grid_coords.array() < 0 || grid_coords.array() >= size_.array()).any()) {
        return std::nullopt;
      }
      return grid_coords;
    }

    Eigen::Vector2f fromGridCoordinates(const GridCoordinates & grid_coords) const;

    Eigen::Vector2f origin() const { return origin_; }
    GridCoordinates size() const { return size_; }
    Eigen::ArrayXXf currentHeightMap() const { return current_height_map_; }
    Eigen::ArrayXXf previousHeightMap() const { return previous_height_map_; }

  private:
    void applySlopeFitting(
      const Eigen::ArrayXXf & intersection_height_map, const GridCoordinates & edge_particle,
      const std::vector<GridCoordinates> & unmovable_neighbors,
      const float slope_fitting_threshold);

    float resolution_;
    float margin_;
    size_t rigidness_;
    float gravity_;
    float time_step_;

    Eigen::Vector2f origin_;
    GridCoordinates size_;

    Eigen::ArrayXXf current_height_map_;
    Eigen::ArrayXXf previous_height_map_;
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> movability_map_;
  };

  std::pair<Eigen::ArrayXXf, std::vector<std::optional<GridCoordinates>>>
  generateIntersectionHeightMap(const Cloth & cloth) const;

  auto toEigen2f(const PointT & point) const { return point.getVector3fMap().template head<2>(); }

  float cloth_resolution_{1.0};
  float cloth_margin_{2.0};
  size_t cloth_rigidness_{3};

  float cloth_initial_z_offset_{0.05};
  float gravity_{9.81};
  float time_step_{0.1};

  size_t max_iterations_{500};
  float iteration_termination_threshold_{0.005};

  bool enable_post_processing_{true};
  float slope_fitting_threshold_{0.3};

  float classification_threshold_{0.5};
};
}  // namespace csf_pcl
