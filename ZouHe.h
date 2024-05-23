#pragma once
#include "Node.h"
#include "Macro.h"
#include "Utility.h"

class ZouHe : public Node
{

private:
	std::vector<double> m_inletVel;
	std::shared_ptr<VelocitySet> velSet;


public:
	Utility Utility;
	ZouHe(int coord_x, int coord_y, int coord_z, std::shared_ptr <VelocitySet> velset,NodeType nt);
	~ZouHe();
	double Density() override;
	std::vector<double> Velocity(std::shared_ptr<VelocitySet> velset) override;
	void SetInletVel(double ux, double uy,double uz);
    void ApplyBc() override;
	void MacroscopicValue(std::shared_ptr<VelocitySet> velSet) override;
	

};

