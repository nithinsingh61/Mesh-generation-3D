#include "global_datatype_3D.h"
#include<numeric>
#include<functional>
/*
//HEADER FILES
#include<CGAL\Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL\Point_3.h>
#include<CGAL\intersection_3.h>
#include<CGAL\squared_distance_3.h>
#include <CGAL\basic.h>
#include <CGAL\Cartesian.h>
#include <CGAL\Triangulation_vertex_base_with_info_3.h>
#include <CGAL\Triangulation_cell_base_with_info_3.h>
//#include <CGAL\Polyhedron_3.h>
#include <CGAL\Triangulation_vertex_base_3.h>
#include <CGAL\Triangulation_data_structure_3.h>
#include <CGAL\Triangulation_geom_traits_3.h>
#include <CGAL\Triangulation_3.h>
#include <CGAL\Delaunay_triangulation_3.h>
#include <CGAL\Mesh_cell_base_3.h>

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

//typedef CGAL::Cartesian<double>												Rep;
typedef CGAL::Exact_predicates_inexact_constructions_kernel						Rep;
typedef CGAL::Filtered_kernel<Rep>											my_K;
typedef CGAL::Point_3<Rep>													Point_3d;
//typedef CGAL::Triangulation_geom_traits_3<Rep>								Gt;
typedef CGAL::Triangulation_geom_traits_3<my_K>								Gt;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned int, Gt>		Vb;
typedef CGAL::Triangulation_cell_base_with_info_3<unsigned int, Gt>			Cb;

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
typedef Delaunay::Vertex_handle												vh;
typedef Delaunay::Facet_circulator											facet_circulator;

typedef my_K::Ray_3											    Ray_3;
typedef my_K::Segment_3										Segment_3;
typedef CGAL::Segment_3<Rep>								Line_3;
typedef my_K::Vector_3		                     		    Vector_3;
typedef my_K::Direction_3									Direction_3;
typedef my_K::Iso_cuboid_3									Iso_cuboid_3;
*/


int MAX_NUM = 30;
int last_vertex_id = -1;
int int_num = 0;



struct MyObject
{
//public:
	double length;
	Edge edg;

//	MyObject(const double& length, Edge edg)
//		: length(length), edg(edg){}
} obj;

//template<typename A, typename B>
bool sort_it(const MyObject &a, const MyObject &b)
{
	return a.length < b.length;
}
//Class Node for Octree
class Node
{
public:
	Node();

	Node(const Point& bottom, const Point& top);
	Node(const Iso_cuboid_3& box);
	virtual ~Node() {}

	std::vector<Point> insidePoints;
	Iso_cuboid_3 cuboid;
	Point minCoor, maxCoor;

//private:
	std::shared_ptr<Node> parent;
	std::vector<std::shared_ptr<Node> > child;

};

Node::Node(): minCoor(0, 0, 0), maxCoor(0, 0, 0), cuboid(minCoor, maxCoor) {}
Node::Node(const Point& bottom, const Point& top): minCoor(bottom), maxCoor(top), cuboid(bottom, top) {}
Node::Node(const Iso_cuboid_3& box): cuboid(box), minCoor(box.min()), maxCoor(box.max()) {}

//Octree class
class Octree
{
public:
	Octree();
	Octree(const Point& bottom, const Point& top);
	Octree(const Iso_cuboid_3& box);
	virtual ~Octree() {}

	//bool is_leafNode();
	bool is_leafNode(const std::shared_ptr<Node>& node);
	void spilt(std::shared_ptr<Node>& node);
	void insertPoints(std::shared_ptr<Node>& node);
	int checkNumPoints(const Node& node);

	//Point leftDown,rightUp;
	std::shared_ptr<Node> rootNode;
	std::shared_ptr<Node> leafNode;
};

Octree::Octree() {}
Octree::Octree(const Point& bottom, const Point& top) {rootNode->maxCoor = top; rootNode->minCoor = Point(bottom);}
Octree::Octree(const Iso_cuboid_3& box): rootNode(new Node(box)) {}


bool Octree::is_leafNode(const std::shared_ptr<Node>& node)
{
	return node->child.size() == 0;
}

void Octree::spilt(std::shared_ptr<Node>& node) {

	//CounterClockWise arrangement of Children startting from left bottom
	Point midPoint((node->minCoor.x() + node->maxCoor.x()) / 2, (node->minCoor.y() + node->maxCoor.y()) / 2, (node->minCoor.z() + node->maxCoor.z()) / 2);


	Iso_cuboid_3 r1(node->minCoor.x(), midPoint.y(), node->minCoor.z(), midPoint.x(), node->maxCoor.y(), midPoint.z());
	Iso_cuboid_3 r2(midPoint.x(), midPoint.y(), node->minCoor.z(), node->maxCoor.x(), node->maxCoor.y(), midPoint.z());
	Iso_cuboid_3 r3(node->minCoor.x(), node->minCoor.y(), node->minCoor.z(), midPoint.x(), midPoint.y(), midPoint.z());
	Iso_cuboid_3 r4(midPoint.x(), node->minCoor.y(), node->minCoor.z(), node->maxCoor.x(), midPoint.y(), midPoint.z());
	Iso_cuboid_3 r5(node->minCoor.x(), midPoint.y(), midPoint.z(), midPoint.x(), node->maxCoor.y(), node->maxCoor.z());
	Iso_cuboid_3 r6(midPoint.x(), midPoint.y(), midPoint.z(), node->maxCoor.x(), node->maxCoor.y(), node->maxCoor.z());
	Iso_cuboid_3 r7(node->minCoor.x(), node->minCoor.y(), midPoint.z(), midPoint.x(), midPoint.y(), node->maxCoor.z());
	Iso_cuboid_3 r8(midPoint.x(), node->minCoor.y(), midPoint.z(), node->maxCoor.x(), midPoint.y(), node->maxCoor.z());

	node->child.push_back(std::shared_ptr<Node>(new Node(r1)));
	node->child.push_back(std::shared_ptr<Node>(new Node(r2)));
	node->child.push_back(std::shared_ptr<Node>(new Node(r3)));
	node->child.push_back(std::shared_ptr<Node>(new Node(r4)));
	node->child.push_back(std::shared_ptr<Node>(new Node(r5)));
	node->child.push_back(std::shared_ptr<Node>(new Node(r6)));
	node->child.push_back(std::shared_ptr<Node>(new Node(r7)));
	node->child.push_back(std::shared_ptr<Node>(new Node(r8)));

	for (int i = 0; i < node->child.size(); ++i) {
		node->child.at(i)->parent = node;
	}
}
//Function to insert points in the node of Octree
void Octree::insertPoints(std::shared_ptr<Node>& node)
{
	if (node->parent != NULL)
	{
		for (int i = 0; i < node->parent->insidePoints.size(); ++i) {
			if ((node->cuboid.has_on_boundary(node->parent->insidePoints.at(i))) || (node->cuboid.has_on_bounded_side(node->parent->insidePoints.at(i)))) {
				node->insidePoints.push_back(node->parent->insidePoints.at(i));
			}
		}
	}

	else {std::cout << "This is a root Node" << std::endl;}
}


//Funtion to check number of points in the node
int Octree::checkNumPoints(const Node& node)
{
	return node.insidePoints.size();
}

//Function to generate Octree
std::shared_ptr<Node> generate_Octree(Octree& tree, std::shared_ptr<Node>& node) {
	std::shared_ptr<Node> temp;
	if (node->insidePoints.size() >= MAX_NUM) {
		tree.spilt(node);

		for (int i = 0; i < node->child.size(); ++i) {
			tree.insertPoints(node->child.at(i));


			temp = node->child.at(i);
			while (temp->insidePoints.size() >= MAX_NUM) {
				temp = generate_Octree(tree, temp);
			}
		}
	}
	return temp;
}

//Printing out the node of the Octree
void printNode1(std::shared_ptr<Node>& node) {
	if (node != nullptr) {
		std::cout << node->cuboid.vertex(0) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << node->cuboid.vertex(1) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << node->cuboid.vertex(2) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << node->cuboid.vertex(3) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << node->cuboid.vertex(4) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << node->cuboid.vertex(5) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << node->cuboid.vertex(6) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << node->cuboid.vertex(7) << " " << "0" << " " << "255" << " " << "0" << std::endl;

	}
	else return;
}
void printNode2(std::shared_ptr<Node>& node) {
	if (node != nullptr) {
		std::cout << ++last_vertex_id << " ";
		std::cout << (++last_vertex_id) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << (last_vertex_id) << " ";
		std::cout << (++last_vertex_id) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << (last_vertex_id) << " ";
		std::cout << (++last_vertex_id) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << (last_vertex_id) << " ";
		std::cout << (last_vertex_id - 3) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << (++last_vertex_id) << " ";
		std::cout << (++last_vertex_id) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << (last_vertex_id) << " ";
		std::cout << (++last_vertex_id) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << (last_vertex_id) << " ";
		std::cout << (++last_vertex_id) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << (last_vertex_id) << " ";
		std::cout << (last_vertex_id - 3) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << (last_vertex_id - 7) << " ";
		std::cout << (last_vertex_id - 2) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << (last_vertex_id - 6) << " ";
		std::cout << (last_vertex_id - 1) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << (last_vertex_id - 5) << " ";
		std::cout << (last_vertex_id) << " " << "0" << " " << "255" << " " << "0" << std::endl;
		std::cout << (last_vertex_id - 4) << " ";
		std::cout << (last_vertex_id - 3) << " " << "0" << " " << "255" << " " << "0" << std::endl;
	}
	else return;
}

//Function to print the quadtree in form of segments of boxes to display using opengl
void printTree1(const Octree& tree) {
	std::queue<std::shared_ptr<Node>> nodeQueue;
	std::shared_ptr<Node> temp;
	if (tree.rootNode == nullptr) {return;}
	else {
		nodeQueue.push(tree.rootNode);
		while (!nodeQueue.empty()) {
			temp = nodeQueue.front();
			printNode1(temp);
			nodeQueue.pop();
			for (int i = 0; i < temp->child.size(); ++i) {
				nodeQueue.push(temp->child.at(i));
			}
		}
	}
}

void printTree2(const Octree& tree) {
	std::queue<std::shared_ptr<Node>> nodeQueue;
	std::shared_ptr<Node> temp;
	if (tree.rootNode == nullptr) {return;}
	else {
		nodeQueue.push(tree.rootNode);
		while (!nodeQueue.empty()) {
			temp = nodeQueue.front();
			printNode2(temp);
			nodeQueue.pop();
			for (int i = 0; i < temp->child.size(); ++i) {
				nodeQueue.push(temp->child.at(i));
			}
		}
	}
}

//Function to get the minimum and Maximum Coordinates of the Point Set
/*void MinMaxCoor(std::ifstream &input,Iso_cuboid_3& box)
{
	double maxX=0,minX=DBL_MAX,maxY=0,minY=DBL_MAX, maxZ=0,minZ=DBL_MAX;
	std::string str;
	while((std::getline(input,str))){
        Point point;
		std::istringstream stream(str);

        while((stream>>point)){
			if(minX>point.x()){minX=point.x();}
			if(maxX<point.x()){maxX=point.x();}
			if(minY>point.y()){minY=point.y();}
			if(maxY<point.y()){maxY=point.y();}
			if(minZ>point.z()){minZ=point.z();}
			if(maxZ<point.z()){maxZ=point.z();}
        }
    }
	//std::cout<<maxX<<" "<<minX<<" "<<maxY<<" "<<minY<<std::endl;
	maxX=maxX+0.05;maxY=maxY+0.05;maxZ=maxZ+0.05;minX=minX-0.05;minY=minY-0.05;minZ=minZ-0.05;
	Iso_cuboid_3 bbox(minX,minY,minZ,maxX,maxY,maxZ);
	box=bbox;
}*/

void MinMaxCoor(std::vector<Point> &input, Iso_cuboid_3& box) {
	//double maxX=0,minX=DBL_MAX,maxY=0,minY=DBL_MAX;
	double maxX = 0, minX = DBL_MAX, maxY = 0, minY = DBL_MAX, maxZ = 0, minZ = DBL_MAX;
	std::vector<Point>::iterator it = input.begin();
	Point point;
	//string str;
	//while((std::getline(input,str))){
	while (it != input.end()) {
		point = *it;
		//istringstream stream(str);
		// while((stream>>point)){
		if (minX > point.x()) {minX = point.x();}
		if (minY > point.y()) {minY = point.y();}
		if (maxX < point.x()) {maxX = point.x();}
		if (maxY < point.y()) {maxY = point.y();}
		if (minZ > point.z()) {minZ = point.z();}
		if (maxZ < point.z()) {maxZ = point.z();}
		++it;
	}
	//}
	//std::cout<<maxX<<" "<<minX<<" "<<maxY<<" "<<minY<<std::endl;
	maxX = maxX + 0.05; maxY = maxY + 0.05; maxZ = maxZ + 0.05; minX = minX - 0.05; minY = minY - 0.05; minZ = minZ - 0.05;
	Iso_cuboid_3 bbox(minX, minY, minZ, maxX, maxY, maxZ);
	box = bbox;
}
//Function to print the Octree in form of segments of boxes to display using opengl
/*void printTree(const Octree& tree) {
	std::queue<std::shared_ptr<Node>> nodeQueue;
	std::shared_ptr<Node> temp;
	if (tree.rootNode==nullptr){return;}
	else {
		nodeQueue.push(tree.rootNode);
		while(!nodeQueue.empty()){
			temp=nodeQueue.front();
			printNode(temp);
			nodeQueue.pop();
			for(int i=0;i<temp->child.size();++i){
				nodeQueue.push(temp->child.at(i));
			}
		}
	}
}*/

//Funtion to Create Delaunay Triangulation from a file of points
/*void create_Delaunay(Delaunay& dt,std::ifstream &input)
{
	std::istream_iterator <Point> begin(input);
   	std::istream_iterator <Point> end;
	dt.insert(begin,end);
}*/
void create_Delaunay(Delaunay& dt, std::vector<Point> &input)
{
	//std::istream_iterator <Point_2> begin(input);
	//std::istream_iterator <Point_2> end;
	dt.insert(input.begin(), input.end());

}
//To assign ids to the Delaunay Triangulation
void Assign_ids(Delaunay& dt)
{
	int id = 0;
	Delaunay::Finite_vertices_iterator	vit;
	for (vit = dt.finite_vertices_begin(); vit != dt.finite_vertices_end(); ++vit) {
		vit->info() = id;
		//cout<<vit->info()<<"id "<<endl;
		++id;
	}
	/*	Delaunay::Finite_facets_iterator	fit;
		for(fit =dt.finite_facets_begin(); fit != dt.finite_facets_end(); ++fit)
				{

			fit->first->vertex((fit->second+1)%4)->info()=id;
			id++;
				}*/
}
//funtion to Extend/Limit the rays to the Bounding Box(The first box in the QuadTree)
//template <typename Type1>

Segment_3 convToSeg(const Iso_cuboid_3& box, const Ray_3& ray) {

	CGAL::Object obj = CGAL::intersection(ray, box);
	//std::cout<<ray<<" Ray in convtoseg"<<std::endl;

	const Point* tempPoint = CGAL::object_cast<Point>(&obj);
	const Segment_3* tempSeg = CGAL::object_cast<Segment_3>(&obj);
	Segment_3 seg;

	if (tempPoint != nullptr) {
		//std::cout<<" In point convtoseg"<<std::endl;

		Segment_3 temp(ray.source(), *tempPoint);
		seg = temp;
		return seg;
	}
	if (tempSeg != nullptr) {
		//std::cout<<" In segment convtoseg"<<std::endl;
		seg = *tempSeg;
		return seg;
	}

	//std::cout<<seg<<" Ray in convtoseg"<<std::endl;
}


//funtion to get a bounding box given a ray(converted to a segment) or a segment for a given threshold value

std::vector<Point> cylinder(const Segment_3& edge)
{
	//double dist = thresh;
	//dist = std::sqrt(dist);

	//Segment_3 tempSeg(Point(0,0,0),Point(0,0,0));
	Point begin = Point (edge.source().x(), edge.source().y(), edge.source().z());
	Point end = Point (edge.end().x(), edge.end().y(), edge.end().z());

	std::vector<Point> X;
	X.reserve(2);
	X.push_back(begin);
	X.push_back(end);

	//X.push_back(CGAL::squared_distance(edge.source(), edge.end()));

	return X;
}

void thresh_Point(const Segment_3& seg, const Point point, const double& thresh, std::vector<Point>& threshPoints)
{
	std::vector<Point> threshSeg = cylinder(seg);
	//threshSeg.erase(threshSeg.begin()+2);
	Point pt1 = (threshSeg.at(0));
	Point pt2 = (threshSeg.at(1));


	bool insideThresh = false;

	//float CylTest_CapsFirst( const Vec3 & pt1, const Vec3 & pt2, float lengthsq, float radius_sq, const Vec3 & testpt )

	double dx, dy, dz;	// vector d  from line segment point 1 to point 2
	double pdx, pdy, pdz;	// vector pd from point 1 to test point
	double dot, dsq;


	dx = pt2.x() - pt1.x();	// translate so pt1 is origin.  Make vector from
	dy = pt2.y() - pt1.y();     // pt1 to pt2.  Need for this is easily eliminated
	dz = pt2.z() - pt1.z();

	pdx = point.x() - pt1.x();		// vector from pt1 to test point.
	pdy = point.y() - pt1.y();
	pdz = point.z() - pt1.z();

	// Dot the d and pd vectors to see if point lies behind the
	// cylinder cap at pt1.x, pt1.y, pt1.z

	dot = pdx * dx + pdy * dy + pdz * dz;

	// If dot is less than zero the point is behind the pt1 cap.
	// If greater than the cylinder axis line segment length squared
	// then the point is outside the other end cap at pt2.

	if ( dot < 0 || dot > (CGAL::squared_distanceC3(pt1.x(), pt1.y(), pt1.z(), pt2.x(), pt2.y(), pt2.z()) ))
	{
		insideThresh = false;
	}
	else
	{
		// Point lies within the parallel caps, so find
		// distance squared from point to line, using the fact that sin^2 + cos^2 = 1
		// the dot = cos() * |d||pd|, and cross*cross = sin^2 * |d|^2 * |pd|^2
		// Carefull: '*' means mult for scalars and dotproduct for vectors
		// In short, where dist is pt distance to cyl axis:
		// dist = sin( pd to d ) * |pd|
		// distsq = dsq = (1 - cos^2( pd to d)) * |pd|^2
		// dsq = ( 1 - (pd * d)^2 / (|pd|^2 * |d|^2) ) * |pd|^2
		// dsq = pd * pd - dot * dot / lengthsq
		//  where lengthsq is d*d or |d|^2 that is passed into this function

		// distance squared to the cylinder axis:

		dsq = ((pdx * pdx + pdy * pdy + pdz * pdz) - dot * dot / (CGAL::squared_distanceC3(pt1.x(), pt1.y(), pt1.z(), pt2.x(), pt2.y(), pt2.z())));

		if ( dsq > CGAL::square(thresh) )
		{
			insideThresh = false;
		}
		else
		{
			insideThresh = true;	// return distance squared to axis
		}
	}
	if (insideThresh == true)
	{
		threshPoints.push_back(point);
	}
}




//Function to get the insidePoints given tree, voronoi edge and threshold
std::vector<Point> insidePoints(const Octree& tree, const Segment_3& edge, const double& thresh) {

	std::vector<std::vector<Point> > tempPoints;
	tempPoints.reserve(5);
	std::queue<std::shared_ptr<Node>> nodeQueue;
	//std::unordered_map<Point,int> mapping;
	//int count=0;

	std::vector<Point> th_Points;
	th_Points.reserve(10);
	std::shared_ptr<Node> temp;
	//temp.reserve(4);


	//std::vector<Segment_3> tempSeg=bBox(edge,thresh);
	//std::vector<double> tempSeg=cylinder(edge);
	Segment_3 tempSeg = edge;
	if (tree.rootNode == nullptr) {return th_Points;}
	else {

		nodeQueue.push(tree.rootNode);
		while (!nodeQueue.empty()) {
			temp = nodeQueue.front();

			if (temp->child.size() == 0)
			{
				tempPoints.push_back(temp->insidePoints);

			}



			for (int i = 0; i < temp->child.size(); ++i) {

				bool boxInside = true;
				if ((CGAL::squared_distance(tempSeg, temp->child.at(i)->cuboid.vertex(0)) <= thresh) ||
				        (CGAL::squared_distance(tempSeg, temp->child.at(i)->cuboid.vertex(1)) <= thresh) ||
				        (CGAL::squared_distance(tempSeg, temp->child.at(i)->cuboid.vertex(2)) <= thresh) ||
				        (CGAL::squared_distance(tempSeg, temp->child.at(i)->cuboid.vertex(3)) <= thresh) ||
				        (CGAL::squared_distance(tempSeg, temp->child.at(i)->cuboid.vertex(4)) <= thresh) ||
				        (CGAL::squared_distance(tempSeg, temp->child.at(i)->cuboid.vertex(5)) <= thresh) ||
				        (CGAL::squared_distance(tempSeg, temp->child.at(i)->cuboid.vertex(6)) <= thresh) ||
				        (CGAL::squared_distance(tempSeg, temp->child.at(i)->cuboid.vertex(7)) <= thresh))
				{
					boxInside = true;
				}
				else
				{
					boxInside = false;
				}

				//if((CGAL::do_intersect(tempSeg.at(0),temp->child.at(i)->cuboid))||(CGAL::do_intersect(tempSeg.at(1),temp->child.at(i)->cuboid))||(CGAL::do_intersect(tempSeg.at(2),temp->child.at(i)->cuboid))||(boxInside==true)){
				if ((CGAL::do_intersect(tempSeg, temp->child.at(i)->cuboid.bbox())) || (boxInside == true)) {
					nodeQueue.push(temp->child.at(i));
				}

			}
			nodeQueue.pop();
		}
	}

	for (int i = 0; i < tempPoints.size(); ++i) {
		for (int j = 0; j < tempPoints[i].size(); ++j) {
			thresh_Point(edge, tempPoints[i].at(j), thresh, th_Points);
		}
	}
	std::vector<Point> th1_Points;
	th1_Points.reserve(10);
	//std::cout<<th_Points.size()<<" Original"<<std::endl;


	return th_Points;
}
/*
Vertex_handle Nearest_Point(const Segment_3& edge,std::vector<Vertex_handle>& points)
{
	double dist, dist_1;
	Vertex_handle p;
	dist = CGAL::squared_distance(edge, points.at(0)->point());
	for(int i=1; i< points.size(); ++i)
	{
		dist_1 = CGAL::squared_distance(edge, points.at(i)->point());
		if(dist < dist_1)
		{
			dist = dist_1;
			p = points.at(i);
		}
	}
	return p;
}*/
Vertex_handle Nearest_Point(const Segment_3& edge, std::vector<Vertex_handle>& points)
{
	double dist, dist_1;
	Vertex_handle p = points.at(0);
	dist = CGAL::squared_distance(edge, points.at(0)->point());

	//cout<<"points size"<<points.size()<<endl;

	if (points.size() > 1)
	{
		for (int i = 1; i < points.size(); ++i)
		{
			dist_1 = CGAL::squared_distance(edge, points.at(i)->point());
			if (dist < dist_1)
			{
				dist = dist_1;
				p = points.at(i);
			}
			//cout<<" "<<i<<" "<<points.at(i)->point()<<endl;
			//cout<<" "<<i<<" "<<p->point()<<endl;

		}
		//cut<< "he"<<endl;
		//cout<<p->point()<<endl;
		return p;
	}
	else return points.at(0);
}

vh projectionPoints(const Segment_3& edge, std::vector<vh>& Proj_vhandle)
{
	//std::cout<<"In projectionPoints"<<std::endl;

	std::vector<Point> threshSeg = cylinder(edge);
	//std::cout<<"got edges end points"<<std::endl;

	std::vector<Point> proj_Points;
	Vertex_handle farthestPoint;
	int farthestPoint_index;

	//Point proj_Point;
	//proj_Points.reserve(10);
	//threshSeg.erase(threshSeg.begin()+2);

	Point pt1 = (threshSeg.at(0));
	Point pt2 = (threshSeg.at(1));

	double dx, dy, dz;	// vector d  from line segment point 1 to point 2
	double pdx, pdy, pdz;	// vector pd from point 1 to test point
	double dot_1, dot_2;
	Point dot_f;
	double dist, dist_1;
	//std::cout<<"find Projection"<<std::endl;
	for (int l = 0; l < Proj_vhandle.size(); ++l)
	{
		Point point = Proj_vhandle.at(l)->point();

		dx = pt2.x() - pt1.x();	// translate so pt1 is origin.  Make vector from
		dy = pt2.y() - pt1.y();     // pt1 to pt2.
		dz = pt2.z() - pt1.z();

		pdx = point.x() - pt1.x();		// vector from pt1 to test point.
		pdy = point.y() - pt1.y();
		pdz = point.z() - pt1.z();

		dot_1 = pdx * dx + pdy * dy + pdz * dz;
		dot_2 = dx * dx + dy * dy + dz * dz;

		Point dot_f1 = Point((dot_1 / dot_2) * dx, (dot_1 / dot_2) * dy, (dot_1 / dot_2) * dz);
		dot_f = Point(pt1.x() + dot_f1.x(), pt1.y() + dot_f1.y(), pt1.z() + dot_f1.z());

		proj_Points.push_back(dot_f);

		//farthestPoint_index = findFarthest(edge, proj_Points);
		//nearest_point.push_back(Nearest_Point(edge, all_cluster.at(l)));
	}
	//std::cout<<"Projection found"<<std::endl;

	dist = CGAL::squared_distanceC3(pt1.x(), pt1.y(), pt1.z(), proj_Points.at(0).x(), proj_Points.at(0).y(), proj_Points.at(0).z());
	//std::cout<<"dist	"<<dist<<std::endl;
	farthestPoint_index = 0;
	for (int i = 1; i < proj_Points.size(); ++i)
	{
		dist_1 = CGAL::squared_distanceC3(pt1.x(), pt1.y(), pt1.z(), proj_Points.at(i).x(), proj_Points.at(i).y(), proj_Points.at(i).z());
		if (dist < dist_1)
		{
			dist = dist_1;
			farthestPoint_index = i;
			//std::cout<<"new dist	"<<dist<<std::endl;
		}
	}
	farthestPoint = Proj_vhandle.at(farthestPoint_index);

	//std::cout<<"............NEXT................"<<std::endl;

	return farthestPoint;
}
/*template<class InputIt, class T, class BinaryOperation>
T accumulate(InputIt first, InputIt last, T init,
             BinaryOperation op)
{
    for (; first != last; ++first) {
        init = op(init, *first);
    }
    return init;
}*/


/*template <typename Type>
void COCONE(Delaunay& dt,std::vector<Point>& sample,std::vector<vector <Point> >& Neighbors,std::vector<Segment_3>& seg,Type ray,int threshold )
{
	std::vector< vector <Point> > TestNeigh;

	for(vit = dt.finite_vertices_begin();vit!=dt.finite_vertices_end();++vit)
	{
		bool notThresh_point=true;
		int distThreshold=CGAL::squared_distance(ray,vit->point());

		if(distThreshold < threshold){
			notThresh_point=false;
		}
		else{
			notThresh_point=true;
		}

		if(notThresh_point==false)
		{
			std::list<Cell_handle> cells;
			dt.finite_incident_cells(vit,std::back_inserter(cells));
			std::list<Cell_handle>::iterator it;

			double max_dist = 0;
			Point pv, c_center;

			//Iterate over incident cells
			for(it = cells.begin(); it != cells.end(); ++it)
			{
				Cell_handle cell = *it;
				c_center = dt.dual(cell);

				//Distance between point and Voronoi vertex
				double dist = CGAL::squared_distanceC3(vit->point().x(),vit->point().y(),vit->point().z(), c_center.x(),c_center.y(), c_center.z());

				if(dist >= max_dist)
				{
					pv = c_center;
					max_dist = dist;
				}
			}


			line_3d Voronoi_edge;
			for(it = cells.begin(); it != cells.end(); ++it)
			{
				Cell_handle cell = *it;
				std::pair<line_3d, unsigned int> edge_table;
				unsigned int num = 0;

				//edges dual to the facets
				for(int i = 0; i < 4; ++i)
				{
					CGAL::Object o = T.dual(cell, (vit->info()+i)%4);
					if(const line_3d * edges = CGAL::object_cast<line_3d>(&o))
					{
						Voronoi_edge = *edges;
						double T1, T2;
						T1 = CGAL::min((CGAL::angle(pv, Voronoi_edge.source(), vit->point())), (CGAL::angle(pv, Voronoi_edge.end(), vit->point())));
						T2 = CGAL::max((CGAL::angle(pv, Voronoi_edge.source(), vit->point())), (CGAL::angle(pv, Voronoi_edge.end(), vit->point())));

						//std::cout<<"T1 is = "<<T1<<" "<<"T2 is = "<<T2<<std::endl;

						double I_max = 2.355, I_min = 0.785;  // theta is 45 degree
						if(((T1 <= I_max) && (T1 >= I_min)) || ((T2 <= I_max) && (T2 >= I_min))){
							num++;
							edge_table = std::make_pair (Voronoi_edge, num);
							//std::cout<<edge_table.first<<" "<<edge_table.second<<" "<<std::endl;
							if(edge_table.second == 3)
								//std::cout<<edge_table.first<<" "<<edge_table.second<<" "<<std::endl;

						}
						else{;}
					}

				}
				//std::cout<<"outside 1"<<std::endl;
			}
		}
	}
}*/

void labeling(std::vector<Edge> Cylinder_Edges, std::vector<int> compo, int n)
{
	std::vector<Edge>::iterator ed;

	for (ed = Cylinder_Edges.begin(); ed != Cylinder_Edges.end(); ++ed)
	{
		for (std::vector<int>::iterator i = compo.begin(); i != compo.end(); ++i)
		{
			if (ed->first->vertex(ed->second)->info() == *i)
				ed->first->vertex(ed->second)->label = n;
			else if (ed->first->vertex(ed->third)->info() == *i)
				ed->first->vertex(ed->third)->label = n;
		}
	}
}

//bool *visited = new bool[1000];
bool visited[100000] = {0};

//std::vector<int> *adj;
//std::vector<int> *adj;
vector <vector <int> > adj;

class Graph
{
	int V;    // No. of vertices
	//list<int> *adj;    // Pointer to an array containing adjacency lists
public:
	Graph(int V);  // Constructor
	void addEdge(int v, int w); // function to add an edge to graph
	//void BFS(std::vector<Vertex_handle> v_handle, int n);  // prints BFS traversal from a given source s
	void BFS(int s, int n);
	//boolean visited =false ;

};

int *level = new int[1000];
Graph::Graph(int V)
{

	this->V = V;
	//adj = new std::vector<int>[V];
	adj.resize(40000);
}

void Graph::addEdge(int v, int w)
{
	adj[v].push_back(w); //Add w to v’s list.
	//adj[w].push_back[v];
}
void Graph::BFS(int s, int n)
{
	std::vector<int > queue ;
	visited[s] = true;
	level[s] = n;
	queue.push_back(s);
	std::vector<int> ::iterator i;
	while (!queue.empty())
	{
		s = queue.front();
		queue.erase(queue.begin());
		//queue.pop_front();
		for (i = adj[s].begin(); i != adj[s].end(); ++i)
		{
			//cout<< *i<<" "<<adj[s].size()<<" "<< *adj[s].begin()<<" "<<s<<" "<<adj.size()<< endl;
			if (!visited[*i])
			{
				visited[*i] = true;
				level[*i] = n;
				queue.push_back(*i);
			}
		}
	}

}


void deCell(Delaunay& dt, const vh& point) // point local unmarking !!
{
	std::vector<Cell_handle> cells;
	dt.incident_cells(point, std::back_inserter(cells));

	std::vector<Cell_handle>::iterator c_it;

	for (c_it = cells.begin(); c_it != cells.end(); ++c_it)
	{
		Cell_handle Cc = *c_it;

		for (int i = 0; i < 4; ++i)
		{
			Cc->rest_facets[i] = false;
			Cell_handle opp_cell = Cc->neighbor(i);

			int opp_index = opp_cell->index(Cc);
			opp_cell->rest_facets[opp_index] = false;
		}
	}
}

void check_duplicate(std::vector<Point> &input)
{
	//std::vector<Point_2> th_pts;
	for (int j = 0; j < input.size(); ++j)
	{
		for (int k = j; k < input.size(); ++k)
		{
			if (j == k) {continue;}
			else
			{
				if (input.at(j) == input.at(k))
				{
					//ThreshPoints.erase(ThreshPoints.begin(),ThreshPoints.begin()+k);	--k;
					input.erase(input.begin() + k);	--k;
				}

			}
		}
	}
}


void unmarkall(Delaunay &dt)
{
	all_cit a_cit;
	a_cit = dt.all_cells_begin();
	for ( ; a_cit != dt.all_cells_end(); ++a_cit)
	{
		for (int i = 0; i < 4; i++)
		{
			//a_fit->correct_segments[i] = false;
			a_cit->rest_facets[i] = false;
		}

	}



}

// return tuple both point and num ka value bhi
pair<Point, int> pointnu(Delaunay& dt_sample_sub, const Facet& fit, Octree& tree)

{



	CGAL::Object o = dt_sample_sub.dual(fit);   // change no one
	const Segment_3 *s = CGAL::object_cast<Segment_3>(&o);
	const Ray_3 *r = CGAL::object_cast<Ray_3>(&o);


	std::vector<Point> ThreshPoints, NewThreshPoints;
	std::vector<vector <Point> > Neighbors;
	Segment_3* temp = new Segment_3;

	if (r)
	{
		//std::cout << "Its a Ray" << std::endl;
		//cout<<"Its a Ray"<<std::endl;
		if (tree.rootNode->cuboid.has_on_bounded_side((*r).source()))
		{
			*temp = convToSeg(tree.rootNode->cuboid, *r);
		}
	}
	if (s)
	{
		//std::cout << "Its a Segment" << std::endl;
		*temp = *s;
	}



	ThreshPoints = insidePoints(tree, *temp, 0.003);
	//cout << "thresh_1278" << " " << ThreshPoints.size() << endl;
	vector<Point> X;
	X	= cylinder(*temp);

	std::ofstream file1;
	file1.open("segs.txt", std::ios_base::app);
	file1 << X[0] << " " << X[1] << endl;
	file1.close();





	check_duplicate(ThreshPoints); // idk
	Delaunay dt_thresh;
	create_Delaunay(dt_thresh, ThreshPoints);
	Assign_ids(dt_thresh);// idk



	if (ThreshPoints.size() <= 3)
	{


		Point p;
		p = Point(0, 0, 0);
		std::pair <Point, int > pointnum(p, 0);
		delete temp;
		return pointnum;

	}


	if (ThreshPoints.size() > 3)
	{

		//cout << "inside thresh points size og do this" << endl;

		Delaunay T;
		create_Delaunay(T, ThreshPoints);
		Assign_ids(T);


		std::vector<Edge> Cylinder_Edges;
		std::vector<Edge>::iterator e_itr;
		std::vector<Point> Cylinder_Points;
		std::vector<Point>::iterator p_it;
		std::vector<Vertex_handle> Cylinder_vertex;
		std:: fill( visited, visited + sizeof(visited), false);  //li

		for (Delaunay::Finite_vertices_iterator vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
		{
			int count = 0;
			std::list<Edge> edges;
			std::list<Edge>::iterator it;
			std::vector<double> edge_length;
			std::vector<MyObject> edgLenvec;
			std::vector<MyObject>::iterator it_e;
			Cylinder_Points.push_back(vit->point());
			T.finite_incident_edges(vit, std::back_inserter(edges));
// ce
			if (edges.size() > 1)
			{
				//std::cout<<"1  "<<std::endl;
				for (it = edges.begin(); it != edges.end(); ++it)
				{
					double length = CGAL::squared_distance(it->first->vertex(it->second)->point(), it->first->vertex(it->third)->point());
					obj.length = length;
					obj.edg = *it;

					edgLenvec.push_back(obj);
				}

				std::sort(edgLenvec.begin(), edgLenvec.end(), sort_it);
				std::list<double> a;
				std::list<double>::iterator a_it;

				for (int i = 0; i < edgLenvec.size() - 1; ++i)
				{
					a.push_back(edgLenvec.at(i + 1).length - edgLenvec.at(i).length);
				}
				//Find the first order difference

				double average = 0;
				for (a_it = a.begin(); a_it != a.end(); ++a_it)
				{
					average = average + *a_it;
				}

				double avg = average / a.size();
				for (a_it = a.begin(); a_it != a.end(); ++a_it)
				{
					if (avg > *a_it)
					{
						++count;
						continue;
					}
					else
						break;
				}

				for (int i = 0; i <= (count); ++i)
				{
					Cylinder_vertex.push_back(edgLenvec.at(i).edg.first->vertex(edgLenvec.at(i).edg.third));
					Cylinder_Edges.push_back(edgLenvec.at(i).edg);
				}
			}
		}

		sort( Cylinder_vertex.begin(), Cylinder_vertex.end() );
		Cylinder_vertex.erase(unique(Cylinder_vertex.begin(), Cylinder_vertex.end()), Cylinder_vertex.end());

		sort( Cylinder_Edges.begin(), Cylinder_Edges.end() );
		Cylinder_Edges.erase(unique(Cylinder_Edges.begin(), Cylinder_Edges.end()), Cylinder_Edges.end());


		Graph gr(Cylinder_vertex.size());

		//std::cout<<"graph start"<<std::endl;
		for (e_itr = Cylinder_Edges.begin(); e_itr != Cylinder_Edges.end(); ++e_itr)
		{
			gr.addEdge(e_itr->first->vertex(e_itr->second)->info(), e_itr->first->vertex(e_itr->third)->info());
		}

		int num = 0;
		for (std::vector<Vertex_handle>::iterator v_itr = Cylinder_vertex.begin(); v_itr != Cylinder_vertex.end(); ++v_itr)
		{
			Vertex_handle vh = *v_itr;
			if (visited[vh->info()] == false)
			{
				//output<<"ll";
				++num;

				//gr.BFS(Cylinder_vertex, num);
				gr.BFS(vh->info(), num);

				//std::cout<<"vh num"<<vh->info()<<" "<<num<<std::endl;

			}
			vh->label = level[vh->info()];
		}

		if (num >= 1)
		{

			std::vector<Vertex_handle> Proj_vhandle;
			for (int i = 1; i <= num; ++i)
			{
				for (std::vector<Vertex_handle>::iterator v_itr = Cylinder_vertex.begin(); v_itr != Cylinder_vertex.end(); ++v_itr)
				{
					Vertex_handle vh = *v_itr;
					if (vh->label == i)
					{
						Proj_vhandle.push_back(vh);
						break;
					}
				}
			}



			//find projection of Proj_vhandle
			vh Farthest_point = projectionPoints(*temp, Proj_vhandle);

			//std::cout<<"Point "<<Farthest_point->point()<<"	Label "<<Farthest_point->label<<std::endl;
			std::vector<Vertex_handle> Farthest_component;
			for (std::vector<Vertex_handle>::iterator v_itr = Cylinder_vertex.begin(); v_itr != Cylinder_vertex.end(); ++v_itr)
			{
				Vertex_handle vh = *v_itr;
				if (vh->label == Farthest_point->label)
					Farthest_component.push_back(vh);
				//std::cout<<"Label	"<<vh->label<<"	Point	"<<vh->point()<<std::endl;
			}
			//std::cout<<"cylinder size	"<<Farthest_component.size()<<std::endl;

			Vertex_handle nearestPoint = Nearest_Point(*temp, Farthest_component);
			//std::cout << "Inserted Point in do this	" << nearestPoint->point() << std::endl;


			std::pair <Point, int > pointnum(nearestPoint->point(), num); //li point use karna ya vertex handle

			delete temp;
			return pointnum ;



		}


	}


}




void singularvertex (Delaunay& dt_sample , Delaunay::Finite_vertices_iterator & fvit, Octree& tree)
{

	std::vector<Facet> restFacets ;
	std::vector<Facet> allfacets;

	dt_sample.finite_incident_facets(fvit, std::back_inserter(allfacets));



	for (std::vector<Facet>::iterator all_facets = allfacets.begin(); (all_facets != allfacets.end()); ++all_facets)
	{

		if ( all_facets->first->rest_facets[all_facets->second] == true)
		{
			restFacets.push_back(*all_facets);
		}


	}



	CGAL::Union_find<Facet> facets;

	facets.insert(restFacets.begin(), restFacets.end());
	//dt_sample.incident_facets(fvit1, std::back_inserter(facets));

	typedef std::map<Vertex_handle, typename CGAL::Union_find<Facet>::handle>  Vertex_Set_map;
	typedef typename Vertex_Set_map::iterator Vertex_Set_map_iterator;

	Vertex_Set_map vsmap;


	for ( CGAL::Union_find<Facet>::iterator it = facets.begin();
	        it != facets.end();
	        ++it) {
		Cell_handle& ch = (*it).first;
		int& i = (*it).second;
		for (int j = 0; j < 3; ++j) {
			Vertex_handle w = ch->vertex(dt_sample.vertex_triple_index(i, j));
			Vertex_handle v = fvit;// same as above !!
			if (w != v) {
				Vertex_Set_map_iterator vsm_it = vsmap.find(w);
				if (vsm_it != vsmap.end()) {
					//cout<<"kk";
					facets.unify_sets(vsm_it->second, it);
					//cout<<"cand"<<endl;
				} else {
					//cout<<"rand"<<endl;
					vsmap.insert(std::make_pair(w, it));
				}
			}
		}
	}

	//int i = facets.size();  // we cast as it cannot be too many
	int j = facets.number_of_sets();
	cout << "number of sets" << j << endl;

	if (j > 1)
	{



		vector<Point> singularvsurfipoints;
		//Vertex_handle singularvfarthest ;
		Point singularvfarthest;
		// vertex hnadle
		for (std::vector<Facet>::iterator singularvfacets = restFacets.begin(); ( singularvfacets != restFacets.end()); ++singularvfacets)
		{

			pair<Point, int> pointnum;
			pointnum = pointnu(dt_sample, *singularvfacets, tree);


			if (pointnum.second >= 1)
				singularvsurfipoints.push_back(pointnum.first);


		}


		int largestdist = 0;
		int dist = 0;


		for (std::vector<Point>::iterator farthestiter = singularvsurfipoints.begin(); farthestiter != singularvsurfipoints.end(); ++farthestiter)
		{

			dist = squared_distance(fvit->point(), *farthestiter);



			if (dist >= largestdist)
			{

				singularvfarthest = *farthestiter;

			}



		}

		Vertex_handle Vhandle = dt_sample.insert(singularvfarthest);


	}
}







bool iterate = true;

int points_inserted = 0;


int main()
{

	Delaunay::Finite_vertices_iterator							vit;
	Delaunay::Finite_vertices_iterator							fvit;

	Delaunay::Finite_cells_iterator								cit;
	Delaunay::Finite_facets_iterator							fit;
	Delaunay::Finite_facets_iterator							ffit;

	Delaunay::Finite_edges_iterator								eit;
	Delaunay::Finite_edges_iterator								feit;

	Delaunay::Facet_circulator									fc;
	Delaunay::Cell_circulator									Cc;

	//ofstream output;

	//output.open("output.txt");

	std::clock_t t1, t2;
	t1 = clock();
	//int stop;
	std::vector<Point> OriginalSample, RandomSample, threshPoints;
	//std::srand(time(NULL));
	std::ifstream input("point.txt");

	std::ofstream outputRandomSample("RandomSample.txt");



	std::ofstream file1;
	file1.open("vert.txt");
	/*std::ofstream file2;
	file2.open("input_verti_info.txt");
	std::ifstream input_verti_info("input_verti_info.txt");
	*/







	int num_of_points = 0;
	std::string data;
	while (getline(input, data)) {
		Point original;
		std::istringstream stream(data);
		while (stream >> original)		{
			OriginalSample.push_back(original);
			++num_of_points;
		}
	}
	input.close();

	Iso_cuboid_3 boundingBox;
	//std::cout<<" Octree start//////////////////////////////////////////"<<std::endl;
	MinMaxCoor(OriginalSample, boundingBox);

	Octree tree(boundingBox);
	tree.rootNode->insidePoints = OriginalSample;
	generate_Octree(tree, tree.rootNode);



	Point p1, p2, p3, p4, p5, p6, p7, p8, p9, p10;
	Point p11, p12, p13, p14, p15, p16, p17, p18, p19, p20;



	//	const char* fname = "points.xyz"; // has all the points

	const char * fname = "allply/e10000.xyz" ; // has points only from sampling by poisson sampling



	std::vector<Point> points_;
	std::ifstream stream(fname);
	if (!stream ||
	        !CGAL::read_xyz_points(
	            stream, std::back_inserter(points_)))
	{
		std::cerr << "Error: cannot read file " << fname << std::endl;
		return EXIT_FAILURE;
	}



	for (unsigned int i = 0; i < points_.size(); i++)


	{


		RandomSample.push_back(points_[i]);


		points_inserted++;


	}





	Delaunay dt_sample;
	create_Delaunay(dt_sample, RandomSample);
	Assign_ids(dt_sample);



	int iterations = 0;


	bool iterate = true;
	while (iterate)
	{
		Assign_ids(dt_sample);



		std::ofstream file2;
		file2.open("input_verti_info.txt");
		std::ifstream input_verti_info("input_verti_info.txt");
		int number_of_skinny_triangle = 0;
		int skinny_bound_zero = 0;
		int skinny_bound_one = 0;
		int skinny_bound_two = 0;
		int skinny_bound_three = 0;
		int skinny_bound_is_two = 0;

		std::vector<std::vector <Point> > Neighbors;

		std::list<Segment_3> Voronoi_edges;
		std::list<Segment_3>::iterator Vor_it;

		//std::list<Facet> dt_facets;
		//std::list<Facet>::iterator facet_it;
		//std::cout<<"inside while"<<std::endl;
		//std::cout<<std::endl;
		//output<<std::endl;

		//std::cout<<"Total facets"<<dt_sample.number_of_facets()<<std::endl;
		std::cout << "Total facets" << dt_sample.number_of_facets() << std::endl;

		//std::cout<<"facet begin points	"<<dt_sample.finite_facets_begin()->first->vertex((dt_sample.finite_facets_begin()->second)%3)->point()<<" "<<dt_sample.finite_facets_begin()->first->vertex((dt_sample.finite_facets_begin()->second +1)%3)->point()<<" "<<dt_sample.finite_facets_begin()->first->vertex((dt_sample.finite_facets_begin()->second+2)%3)->point()<<std::endl;
		//std::cout<<"\n";
		//std::cout<<"facet end points	"<<dt_sample.finite_facets_end()->first->vertex((dt_sample.finite_facets_end()->second)%3)->point()<<" "<<dt_sample.finite_facets_end()->first->vertex((dt_sample.finite_facets_end()->second +1)%3)->point()<<" "<<dt_sample.finite_facets_end()->first->vertex((dt_sample.finite_facets_end()->second+2)%3)->point()<<std::endl;
		bool rays_ok = true;
		bool segs_ok = true;
		bool top_ok = true;


		for (fit = dt_sample.finite_facets_begin(); fit != dt_sample.finite_facets_end(); ++fit)
		{
			iterate = false;



			//std::cout << "Number of Point" << endl;

			pair<Point, int> pointnum;
			pointnum = pointnu(dt_sample, *fit, tree);




			if (pointnum.second > 1)
			{
				top_ok = false;
				Vertex_handle tempVhandle = dt_sample.insert(pointnum.first);
				deCell(dt_sample, tempVhandle);  // idk li
				iterate = true;
				++points_inserted;
				std::cout << "Point inserted" << std::endl;
				//std::cout<<"Number of Points	"<<dt_sample.number_of_vertices()<<std::endl;
				break;
			}



			if (pointnum.second == 0) //idk li
			{


				iterate = true;



			}

			if (pointnum.second == 1)
			{


				//std::cout << "inside youre " << std::endl;
				//std::cout<<"3 "<<fit->first->vertex((fit->second+1)%4)->point()<<" "<<fit->first->vertex((fit->second+2)%4)->point()<<" "<<fit->first->vertex((fit->second+3)%4)->info()<<endl;
				fit->first->rest_facets[fit->second] = true;
				Cell_handle opp_cell = fit->first->neighbor(fit->second);
				int opp_index = opp_cell->index(fit->first);
				opp_cell->rest_facets[opp_index] = true;

				file2 << "3 " << fit->first->vertex((fit->second + 1) % 4)->info() << " " << fit->first->vertex((fit->second + 2) % 4)->info() << " " << fit->first->vertex((fit->second + 3) % 4)->info() << endl;
				//std::cout << fit->first->vertex((fit->second + 1) % 4)->point() << " " << fit->first->vertex((fit->second + 2) % 4)->point() << " " << fit->first->vertex((fit->second + 3) % 4)->point() << endl;

				file1 << fit->first->vertex((fit->second + 1) % 4)->info() << " " << fit->first->vertex((fit->second + 1) % 4)->point() << " " << fit->first->vertex((fit->second + 2) % 4)->info() << " " << fit->first->vertex((fit->second + 2) % 4)->point() << " " << fit->first->vertex((fit->second + 3) % 4)->info() << " " << fit->first->vertex((fit->second + 3) % 4)->point() << endl;

				//rest_facets++;
				//file1<<"3 "<<fit->first->vertex((fit->second+1)%4)->info()<<" "<<fit->first->vertex((fit->second+2)%4)->info()<<" "<<fit->first->vertex((fit->second+3)%4)->info()<<endl;
				//file2<<fit->first->vertex((fit->second+1)%4)->point()<<" "<<fit->first->vertex((fit->second+2)%4)->point()<<" "<<fit->first->vertex((fit->second+3)%4)->point()<<endl;
			}
			//}//vit
			//if(iterate == true)
			//	break;

		}//facets

		//while





		if (top_ok == false)
		{

			unmarkall(dt_sample);
		}

		if (top_ok == true)
		{	//3





			if (iterations == 0)
			{

				return 0;

			}




			cout << "ALL iterations " << iterations++ << endl;

			cout << "-----------------Manifoldness --------" << endl;

			//  --------------- now in manifold -------------------------
			iterate = false;

			bool flag = true;

			while (flag)
			{
				//cout<< "startinf in flag ";
				//cout<<dt_sample.size()<<endl;
				cout << dt_sample.number_of_vertices() << endl;
				fvit = dt_sample.finite_vertices_begin();

				for (; fvit != dt_sample.finite_vertices_end(); ++fvit) // first break out
				{

					flag = false;
					map <int , int > counter ;
					//working
					//cout<< dt_sample.degree(fvit);
					//cout<<"indisde dsampel"<<endl;
					std::vector<Facet> Facetvec;
					//// cout << "the point is " << fvit->point() << endl;




					dt_sample.finite_incident_facets(fvit, std::back_inserter(Facetvec));

					// singular edge ka code ;;
					int k = 0;
					for (std::vector<Facet>::iterator Facetit2 = Facetvec.begin(); (Facetit2 != Facetvec.end()); ++Facetit2)
					{

						//cout<< "On a facet"<<endl;

						if (Facetit2->first->rest_facets[Facetit2->second] == true)
						{



							std::cout << "3 " << Facetit2->first->vertex((Facetit2->second + 1) % 4)->info() << " " << Facetit2->first->vertex((Facetit2->second + 2) % 4)->info() << " " << Facetit2->first->vertex((Facetit2->second + 3) % 4)->info() << endl;
							//file1 << Facetit2->first->vertex((Facetit2->second + 1) % 4)->info() << " " << Facetit2->first->vertex((Facetit2->second + 1) % 4)->point() << " " << Facetit2->first->vertex((Facetit2->second + 2) % 4)->info() << " " << Facetit2->first->vertex((Facetit2->second + 2) % 4)->point() << " " << Facetit2->first->vertex((Facetit2->second + 3) % 4)->info() << " " << Facetit2->first->vertex((Facetit2->second + 3) % 4)->point() << endl;


							std::cout << Facetit2->first->vertex((Facetit2->second + 1) % 4)->point() << " " << Facetit2->first->vertex((Facetit2->second + 2) % 4)->point() << " " << Facetit2->first->vertex((Facetit2->second + 3) % 4)->point() << endl;

							//// cout << fvit->info() << endl;


							if (fvit->info() != Facetit2->first->vertex((Facetit2->second + 1) % 4)->info())
							{
								counter[Facetit2->first->vertex((Facetit2->second + 1) % 4)->info()]++;
								int count = counter[Facetit2->first->vertex((Facetit2->second + 1) % 4)->info()];
								//cout << "Inside 1" << " " << count << endl;
								if (count > 2)
								{
									flag = true;


									//cout << "ram1";
									iterate = true;
									//	cout << "t_points " << " " << ThreshPoints_.size() << endl;
									//dothis(dt_sample, *Facetit2, tree, *temp_, flag);
									pair<Point, int> pointnum;

									pointnum = pointnu(dt_sample, *Facetit2, tree);
									//cout<< "Point bhi dekhun "<< pointnum.first<<endl;


									if (pointnum.second >= 1) //
									{
										int pt = dt_sample.number_of_vertices();
										Vertex_handle tempVhandle = dt_sample.insert(pointnum.first);

										if (pt == dt_sample.number_of_vertices())
										{

											flag = false;

										}
										deCell(dt_sample, tempVhandle);  // idk li
										iterate = true;
										++points_inserted;				 // point inserted
										break;//idk li
									}
								}

							}

							if (fvit->info() != Facetit2->first->vertex((Facetit2->second + 2) % 4)->info())
							{

								counter[Facetit2->first->vertex((Facetit2->second + 2) % 4)->info()]++;
								int count = counter[Facetit2->first->vertex((Facetit2->second + 2) % 4)->info()];
								//cout << "Inside 2" << " " << count << endl;
								if (count > 2)
								{
									flag = true;
									//cout << "ram2";
									iterate = true;
									//	cout << "t_points " << " " << ThreshPoints_.size() << endl;
									//dothis(dt_sample, *Facetit2, tree, *temp_, flag); // point inserted

									pair<Point, int> pointnum;

									pointnum = pointnu(dt_sample, *Facetit2, tree);
									//cout<< "Point bhi dekhun "<< pointnum.first<<endl;



									if (pointnum.second >= 1)
									{
										int pt = dt_sample.number_of_vertices();
										Vertex_handle tempVhandle = dt_sample.insert(pointnum.first);

										//Vertex_handle tempVhandle = dt_sample.insert(pointnum.first);
										if (pt == dt_sample.number_of_vertices())
										{

											flag = false; // change

										}
										deCell(dt_sample, tempVhandle);  // idk li
										iterate = true;
										++points_inserted;


										break;
									}

								}

							}

							if (fvit->info() != Facetit2->first->vertex((Facetit2->second + 3) % 4)->info())
							{

								counter[Facetit2->first->vertex((Facetit2->second + 3) % 4)->info()]++;
								int count = 	counter[Facetit2->first->vertex((Facetit2->second + 3) % 4)->info()];
								//cout << "Inside 3" << " " << count << endl;
								if (count > 2)
								{
									flag = true;
									//cout << "ram3";
									iterate = true;
									//	cout << "t_points " << " " << ThreshPoints_.size() << endl;
									//dothis(dt_sample, *Facetit2, tree, *temp_, flag); // point inserted

									pair<Point, int> pointnum;

									pointnum = pointnu(dt_sample, *Facetit2, tree);

									//cout<< "Point bhi dekhun "<< pointnum.first<<endl;

									if (pointnum.second >= 1)
									{
										int pt = dt_sample.number_of_vertices();
										Vertex_handle tempVhandle = dt_sample.insert(pointnum.first);

										if (pt == dt_sample.number_of_vertices())
										{

											flag = false;

										}

										//Vertex_handle tempVhandle = dt_sample.insert(pointnum.first);
										deCell(dt_sample, tempVhandle);  // idk li
										iterate = true;
										++points_inserted;
										break;
									}

								}
							}

						}

					}



					int pt = dt_sample.number_of_vertices();
					//Vertex_handle v = *fvit;
					singularvertex(dt_sample , fvit, tree);

					if (pt != dt_sample.number_of_vertices())
					{
						iterate = true;
						flag = true; //
						break;

					}

					if (flag == true)
					{
						unmarkall(dt_sample);


						counter.clear();

						break; // goes to while reintialzes vertex iteator due to invalidation

					}



				}

				if (flag == true)
				{



					//counter.clear();

					break; // goes to while reintialzes vertex iteator due to invalidation

				}
			} // while anaomly exists

			//}
		}//if rays == ok and seg ==ok


	}//while(iterate)
	cout << "  .......  END  .......  " << endl;
}//main

