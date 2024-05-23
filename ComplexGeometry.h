#pragma once
#ifndef _MESH_GENERATOR_H_
#define _MESH_GENERATOR_H_

#include "Struct.h"
#include <tuple>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <fstream>
#include "Struct.h"
#include <unordered_set>
#include <mpi.h>
#include <omp.h>
#include <map>
#include "GeometryWriter.h"
#include "GeometryReader.h"
#include "mpi_call.h"
#include "Utility.h"
#include "VelocitySet.h"

class Geometry {

public:

	Writer Writer;
	Reader Reader;
	Utility Utility;
	MPI_CALL<MPI_Datatype> mpi_caller;
	std::vector< std::vector< std::vector<std::shared_ptr<Point>>>> pointsBox;
	std::vector<Triangle>globalTriangles;
	std::vector<Triangle>localTriangles;
	std::string stlFileName;
	BoundingBox boundingbox;
	Octree octree;
	int maxDepth;
	std::vector<Point> boundaryPoints;
	std::vector<Point> InletboundaryPoints;
	std::vector<Point> OutletboundaryPoints;
	std::vector<Point> WallboundaryPoints;
	int globalNumtris;
	std::vector<int> GeomSize;
	double spacing;
	double rayLength;
	int myrank, num_procs;

	std::shared_ptr<VelocitySet> velSet_;

	Geometry()=default;
	Geometry(const std::string& FileName, std::shared_ptr<VelocitySet> velSet, double spacing, int mdepth);
	~Geometry();
	void GetBoundingbox();
	void DistributeTriangles();
	bool isPointInsideBoundingBox(double x, double y, double z, const BoundingBox& box);
	std::vector<Triangle> findTrianglesInSubBox(const BoundingBox& subBox);
	void buildOctree(Octree& octree, int depth);
	void generatePointsInBoundingBox();
	bool DoesLineBoxIntersect(BoundingBox& box, Line& line, double t0, double t1);
	std::vector<Octree> findDeepestIntersectedLeaves(Line& line, Octree& octree);
	bool doLinesIntersect(const Line& line1, const Line& line2);
	bool isPointInsideTriangle(double x, double y,const Triangle& triangle);
	int countIntersectinsTriangle(Line& line, std::vector<Triangle>& posintersectedtri);
	std::vector<Point> findBoundaryPoints();
	void adjustLineCoordinates(Line& line, BoundingBox& boundingbox, double rayOffset, int directionSwitch);
	void labelPoints();
	std::vector<Point> excludeInletOutlet();
	void getInOutBP();
	void JustifyInOutWall();
};

#endif
