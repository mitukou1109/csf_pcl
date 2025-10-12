# csf_pcl

A standalone **PCL (Point Cloud Library)** implementation of the _Cloth Simulation Filter (CSF)_<sup>[1]</sup> for ground and non-ground point cloud separation.

If you're using **ROS 2**, please check out the companion package:
👉 [**csf_ros**](https://github.com/mitukou1109/csf_ros.git)

## 📦 Installation

### Using **colcon** (recommended for ROS 2 workspaces)

Simply clone this repository into your workspace’s `src` directory and build it as usual.
See [here](https://github.com/mitukou1109/csf_ros?tab=readme-ov-file#usage) for instructions.

### Standalone build

For standalone (non-ROS) use, build and install the library as follows:

```bash
git clone https://github.com/mitukou1109/csf_pcl.git
cd csf_pcl
mkdir build && cd build
cmake ..
make -j$(nproc)
sudo make install
```

## ⚙️ Usage in CMake projects

You can easily integrate **csf_pcl** into your own CMake-based projects like the following example:

```cmake
cmake_minimum_required(VERSION 3.8)
project(my_project)

find_package(csf_pcl REQUIRED)

add_executable(my_target main.cpp)
target_link_libraries(my_target PRIVATE csf_pcl::csf_pcl)
```

## 🧩 Reference

[1] W. Zhang et al., “An Easy-to-Use Airborne LiDAR Data Filtering Method Based on Cloth Simulation,” Remote Sensing, vol. 8, no. 6, p. 501, June 2016, doi: 10.3390/rs8060501.

## 📝 License

This project is licensed under the [BSD License](LICENSE).
