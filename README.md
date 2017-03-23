# Mesh-generation-3D

Generating mesh for point cloud data without explicilty reconstructing the whole surface exploiting the fact that well sampled data inherently captures the topology and geometry of the scanned object.  


## Getting Started

Clone --> Build --> Input pointcloud --> visualize

### Prerequisites

* [CGAL](http://www.cgal.org/) - Computational Geometry and Algorithms library
* [BOOST](http://www.boost.org/)
* [CUDA](https://developer.nvidia.com/cuda-downloads)
* [Meshlab](www.meshlab.net)


### Building

```
$ cgal_create_cmake_script
$ cmake .
$ make
$ ./Improve_DT
```
* In case of build failure try clearing cache of cmake

### Progress pics 

* [Drive](https://drive.google.com/open?id=0BzT2RJmWrKsiT3dTWFdWV0p6Rmc)
* [Imgur](http://imgur.com/a/23XCR)



## To-do list 
- [x] Topology recover 
- [ ] Clean the code 
- [ ] Sliver exudation,Hole repairing 
- [ ] Mesh refinement 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details


