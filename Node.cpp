#include "Node.h"

Node::Node(int coord_x, int coord_y, int coord_z, NodeType nt)
	: x_position(coord_x), y_position(coord_y), z_position(coord_z), NodeType_(nt)
{
	NodeType_ = nt;
	m_eqdistributions.resize(_nLatNodes);
}


Node::~Node()
{
}

void Node::Stream(int dir)
{

	m_neighbours[dir]->m_distributions[dir] = m_newDistributions[dir];
	
}



void Node::Initialize(std::shared_ptr<VelocitySet> velSet, int dir)
{
	double rho = 1.0;   

	m_distributions.push_back(velSet->GetWeight(dir));
	m_newDistributions.push_back(0.0);
	velocity.resize(3);
}

double Node::Density() {

	double rho = 0.0;
	for (int dir = 0; dir < _nLatNodes; dir++) {

		rho += this->m_distributions[dir];
	}
	return rho;
}

std::vector<double> Node::Velocity(std::shared_ptr<VelocitySet>velset) {

	std::vector<double>vel;
	vel.resize(3, 0.0);

	for (int dir = 0; dir < _nLatNodes; dir++) {
		for (int dim = 0; dim < _dim; dim++) {
			vel[dim] += velset->GetDirection(dir)[dim] * this->m_distributions[dir];
		}
	}

	for (int dim = 0; dim < _dim; dim++) {
		vel[dim] /= Density();
	}
	return vel;
}

double Node::Equilibrium(std::shared_ptr<VelocitySet> velset, int dir) {

	double du = Velocity(velset)[0] * velset->GetDirection(dir)[0] + Velocity(velset)[1] * velset->GetDirection(dir)[1] + Velocity(velset)[2] * velset->GetDirection(dir)[2];
	double u_sqr = Velocity(velset)[0] * Velocity(velset)[0] + Velocity(velset)[1] * Velocity(velset)[1] + Velocity(velset)[2] * Velocity(velset)[2];

	return velset->GetWeight(dir) * Density() * (1 + 3.0 * du + 9.0 / 2.0 * du * du - 3.0 / 2.0 * u_sqr);

}



void Node::MacroscopicValue(std::shared_ptr<VelocitySet> velSet) {

	density = 0;
	velocity[0] = 0.0;
	velocity[1] = 0.0;
	velocity[2] = 0.0;

	for (int dir = 0; dir < _nLatNodes; dir++) {

		density += m_distributions[dir];
		velocity[0] += m_distributions[dir] * velSet->GetDirection(dir)[0];
		velocity[1] += m_distributions[dir] * velSet->GetDirection(dir)[1];
		velocity[2] += m_distributions[dir] * velSet->GetDirection(dir)[2];

	}

	velocity[0] /= density;
	velocity[1] /= density;
	velocity[2] /= density;
}


double Node::Pressure(std::shared_ptr<VelocitySet> velSet)
{
	double c_square = velSet->Get_c_square(0);
	double pressure = c_square * density;
	return pressure;
}


void Node::BounceBack(int dir, int opp_dir)
{

	m_distributions[opp_dir] = m_distributions[dir];
}


void Node::ApplyBc()
{

}