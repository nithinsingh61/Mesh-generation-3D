#ifndef GLOBAL_DATATYPES_H
#define GLOBAL_DATATYPES_H

//HEADER FILES
#include <CGAL/Point_3.h>
#include <CGAL/intersection_3.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
//#include <CGAL\Polyhedron_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_geom_traits_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Polyhedron_3.h>


#include <iostream>
#include <list>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <vector>
#include <stack>
#include <queue>
#include <unordered_map>
#include <sstream>
#include <time.h>
#include <cmath>


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>





#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/IO/write_xyz_points.h>
#include <CGAL/Union_find.h>



//#undef CGAL_CARTESIAN_H

using namespace std;

//new class for 2d triangle

//Define a new vertex base class for storing the vertices in 2D triangulation of each facet of the input

template <class Gt, class Vb = CGAL::Triangulation_vertex_base_with_info_3<unsigned int,Gt> >
class Delaunay_vertex_3D : public Vb
{
	typedef Vb											Base;

public:
	typedef typename Vb::Vertex_handle					Vertex_handle;
	typedef typename Vb::Cell_handle					Cell_handle;
	typedef typename Vb::Point							Point;


	//rebind mechanism
	template < typename TDS2>
	struct	Rebind_TDS {
		typedef	typename Vb::template Rebind_TDS<TDS2>::Other	Vb3;
		typedef Delaunay_vertex_3D<Gt, Vb3>						Other;
	};
	bool visited;
	int label;

public:
	Delaunay_vertex_3D():Base() {init();}
	Delaunay_vertex_3D( const Point& p) : Base(p) {init();}
	Delaunay_vertex_3D( const Point& p, Cell_handle c) : Base(p, c) {init();}
	Delaunay_vertex_3D( Cell_handle c): Base(c) {init();}
	~Delaunay_vertex_3D() {}
	int id;

//private:
	inline void init()
	{
		id=0;
		visited = false;
		label = 0;
	}
};

//Define a new face base class for storing the facets of the 2D triangulation of each facet of the input

template <class Gt, class Cb = CGAL::Triangulation_cell_base_with_info_3<unsigned int, Gt> >
class Delaunay_cell_3D : public Cb
{

	typedef Cb											Base;
	typedef typename Cb::Triangulation_data_structure   TDS;
	
public:
	typedef Gt											Geom_traits;
	typedef TDS											Triangulation_data_structure;
	typedef typename TDS::Vertex_handle					Vertex_handle;
	//typedef typename TDS::Facet_handle					Facet_handle;
	typedef typename TDS::Cell_handle					Cell_handle;

	//rebind mechanism
	template < typename TDS2>
	struct	Rebind_TDS {
		typedef	typename Cb::template Rebind_TDS<TDS2>::Other	Cb3;
		typedef Delaunay_cell_3D<Gt, Cb3>						Other;
	};

	bool cell_is_inside;
	bool rest_facets[4];
	bool large_edge[6];

	Delaunay_cell_3D():Base() {init();}
	Delaunay_cell_3D( Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) : Base(v0, v1, v2, v3) {init();}
	Delaunay_cell_3D( Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, Vertex_handle v3, Cell_handle n0, Cell_handle n1, Cell_handle n2, Cell_handle n3) : Base(v0, v1, v2, v3, n0, n1, n2, n3) {init();}
	~Delaunay_cell_3D() {}
	inline void init() {

		rest_facets[0] = false;
		rest_facets[1] = false;
		rest_facets[2] = false;
		rest_facets[3] = false;
		
		cell_is_inside = true;
	}
};

   

//typedef CGAL::Cartesian<double>											Rep;
typedef CGAL::Exact_predicates_inexact_constructions_kernel					Rep;
typedef CGAL::Filtered_kernel<Rep>											my_K;
typedef CGAL::Point_3<Rep>													Point_3d;
//typedef CGAL::Triangulation_geom_traits_3<Rep>							Gt;
typedef CGAL::Triangulation_geom_traits_3<my_K>								Gt;
//typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned int, Gt>		Vb;
//typedef CGAL::Triangulation_cell_base_with_info_3<unsigned int, Gt>		Cb;
typedef Delaunay_vertex_3D<my_K>											Vb;
typedef Delaunay_cell_3D<my_K>												Cb;
typedef CGAL::Triangulation_data_structure_3<Vb,Cb>							Tds;
typedef CGAL::Triangulation_3<Gt,Tds>										Triangulation;
typedef CGAL::Delaunay_triangulation_3<Gt,Tds>								Delaunay;
//typedef Delaunay::Point														Point_3;
typedef Delaunay::Cell_handle												Cell_handle;

//typedef CGAL::Facet_handle													Facet_handle;
typedef Delaunay::Vertex_handle												Vertex_handle;
//typedef CGAL::Segment_3<Rep>												line_3d;
typedef my_K::Point_3														Point;
//typedef Delaunay::Point														Point;
typedef Delaunay::Facet														Facet;
typedef Delaunay::Edge														Edge;
typedef Delaunay::Vertex_handle												vh;
typedef Delaunay::Facet_circulator											facet_circulator;
typedef Delaunay::Edge_iterator												e_it;

//typedef Delaunay::Edge_iterator												e_it;
//typedef Deaunay::Cell_circulator											Cell_circulator;

typedef my_K::Ray_3											Ray_3;
typedef my_K::Segment_3										Segment_3;
typedef CGAL::Segment_3<Rep>								Line_3;
typedef my_K::Vector_3		                     		    Vector_3;
typedef my_K::Direction_3									Direction_3;
typedef my_K::Iso_cuboid_3									Iso_cuboid_3;
typedef Delaunay::All_cells_iterator							all_cit;
//extern vector<int> adj[1000];
extern vector <vector <int> > adj; 

#endif //GLOBAL_DATATYPE_H
