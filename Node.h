#pragma once
#include <iostream>
#include "VelocitySet.h"
#include "include.h"
#include "Struct.h"
#include "Macro.h"
class Node
{



public:
	std::vector<double> m_distributions;
	std::vector<double> m_newDistributions;
	std::vector<double> m_eqdistributions;
	std::shared_ptr<Node> m_neighbours[_nLatNodes];
	int x_position;
	int y_position;
	int z_position;
	double density;
	std::vector<double> velocity;
	NodeType NodeType_;
	BoundaryType BoundaryType_;
	double WSS_x;
	double WSS_y;
	double WSS_z;
	double WSS_Mag;
	Node() = default;
	Node(int coord_x, int coord_y, int coord_z, NodeType nt);
	~Node();

	void Stream(int dir);
	void Initialize(std::shared_ptr<VelocitySet> velSet, int dir);
	double Equilibrium(std::shared_ptr<VelocitySet> velset, int dir);
	virtual double Density();
	virtual std::vector<double> Velocity(std::shared_ptr<VelocitySet> velset);
	virtual void MacroscopicValue(std::shared_ptr<VelocitySet> velSet);
	double Pressure(std::shared_ptr<VelocitySet> velSet);
	void BounceBack(int dir, int opp_dir);
	virtual void ApplyBc();


};

