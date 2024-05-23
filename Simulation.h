#pragma once

#include "VelocitySet.h"
#include "Subdomain.h"
#include "OutputFile.h"
#include "ZouHe.h"
#include <vector>
#include "pulsatile.h"
#include "GeometryWriter.h"


class Simulation
{

public:
	Writer Writer;
	int myrank_, num_procs_;
	SubDomain SubDomain_;
	Pulsatile Pulsatile_;
	double pressure_;
	double m_relaxation;
	double m_dt;
	double m_omega;
	std::vector<double> wall_shear_stress;
	std::vector<double> InletVelocity_;
	Simulation(SubDomain& SubDomain, double tau, double pressure, double del_t);
	~Simulation();
	void Collision();
	void Streaming();
	void ApplyZouHeBc();
	void ApplyOutletBc();
	void ApplyWallBc();
	void FaceCommunicate(int iter);
	void EdgeCommunicate();
	double getInletVelocity(double x, double period);
	void Iterate();
	void SetRelaxation(double tau);
	void WriteOutput(int iter);

};


