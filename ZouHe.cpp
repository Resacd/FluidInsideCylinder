#include "ZouHe.h"
#include <cmath>
#include <iostream>

ZouHe::ZouHe(int coord_x, int coord_y, int coord_z, std::shared_ptr <VelocitySet> velset, NodeType nt)
    : Node(coord_x, coord_y, coord_z, nt) 
{

    SetInletVel(0.01, 0.0,0.0);
    velSet = velset;

}

ZouHe::~ZouHe() {


}


double ZouHe::Density() {

    return 1;
}


std::vector<double> ZouHe::Velocity(std::shared_ptr<VelocitySet> velset) {

    return m_inletVel;

}


void ZouHe::SetInletVel(double ux, double uy, double uz)
{
    m_inletVel.push_back(ux);
    m_inletVel.push_back(uy);
    m_inletVel.push_back(uz);
}


void ZouHe::MacroscopicValue(std::shared_ptr<VelocitySet> velSet) {



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
    //density = 1;
    //velocity = m_inletVel;

}




#ifdef D3Q27

void ZouHe::ApplyBc()
{
    double rho = ((m_distributions[0] + m_distributions[3] + m_distributions[4]       +
                   m_distributions[5] + m_distributions[6] + m_distributions[15]      +
                   m_distributions[16] + m_distributions[17] + m_distributions[18])   +
                   2 * (m_distributions[2] + m_distributions[9] + m_distributions[10] +
                   m_distributions[12] + m_distributions[14] + m_distributions[20]    +
                   m_distributions[22] + m_distributions[24] + m_distributions[26])) / (1 - m_inletVel[0]);

    double Ny = (1.0 / 6) * (m_distributions[3] + m_distributions[15] + m_distributions[16]
        - m_distributions[4] - m_distributions[17] - m_distributions[18]) - (1.0 / 9) * (rho * m_inletVel[1]);

    double Nz = (1.0 / 6) * (m_distributions[5] + m_distributions[15] + m_distributions[17]
        - m_distributions[6] - m_distributions[16] - m_distributions[18]) - (1.0 / 9) * (rho * m_inletVel[2]);


    m_distributions[1] = m_distributions[2]   + (4 / 9.0) * rho  * m_inletVel[0];

    m_distributions[7] = m_distributions[10]  + (1 / 9.0) * rho  * (m_inletVel[0] + m_inletVel[1]) - Ny;
                                                                    
    m_distributions[8] = m_distributions[9]   + (1 / 9.0) * rho  * (m_inletVel[0] - m_inletVel[1]) + Ny;
                                                                    
    m_distributions[11] = m_distributions[14] + (1 / 9.0) * rho  * (m_inletVel[0] + m_inletVel[2]) - Nz;
                                                                    
    m_distributions[13] = m_distributions[12] + (1 / 9.0) * rho  * (m_inletVel[0] - m_inletVel[2]) + Nz;
                                                                    
    m_distributions[19] = m_distributions[20] + (1 / 36.0) * rho * (m_inletVel[0] + m_inletVel[1] + m_inletVel[2]) - Ny - Nz;
                                                         
    m_distributions[21] = m_distributions[22] + (1 / 36.0) * rho * (m_inletVel[0] + m_inletVel[1] - m_inletVel[2]) - Ny + Nz;
                                                         
    m_distributions[23] = m_distributions[24] + (1 / 36.0) * rho * (m_inletVel[0] - m_inletVel[1] + m_inletVel[2]) + Ny - Nz;
                                                         
    m_distributions[25] = m_distributions[26] + (1 / 36.0) * rho * (m_inletVel[0] - m_inletVel[1] - m_inletVel[2]) + Ny + Nz;


    MacroscopicValue(velSet);
}

#endif

#ifdef D3Q19

void ZouHe::ApplyBc()
{

    //double rhoN = 1.25;
    //m_inletVel[0] = 1 - (1.0/rhoN)*(m_distributions[2] + m_distributions[3] + m_distributions[4] + m_distributions[5] +
    //    m_distributions[14] + m_distributions[15] + m_distributions[16] + m_distributions[17] + m_distributions[18] +
    //    2 * (m_distributions[1] + m_distributions[10] + m_distributions[11] + m_distributions[12] + m_distributions[13]));
    
    double rhoN = (m_distributions[2] + m_distributions[3] + m_distributions[4] + m_distributions[5] +
        m_distributions[14] + m_distributions[15] + m_distributions[16] + m_distributions[17] + m_distributions[18] +
        2 * (m_distributions[1] + m_distributions[10] + m_distributions[11] + m_distributions[12] + m_distributions[13]))/(1 - m_inletVel[0]);

    double Ny = (1.0 / 2) * (m_distributions[2] + m_distributions[14] + m_distributions[15]
        - m_distributions[3] - m_distributions[16] - m_distributions[17]) - (1.0 / 3) * (rhoN * m_inletVel[1]);

    double Nz = (1.0 / 2) * (m_distributions[4] + m_distributions[16] + m_distributions[14]
        - m_distributions[5] - m_distributions[15] - m_distributions[17]) - (1.0 / 3) * (rhoN * m_inletVel[2]);


  m_distributions[0] = m_distributions[1]   + (1 / 3.0) * rhoN * m_inletVel[0];

  m_distributions[7] = m_distributions[10]  + (1 / 6.0) * rhoN * (m_inletVel[0] - m_inletVel[1]) + Ny;

  m_distributions[6] = m_distributions[11]  + (1 / 6.0) * rhoN * (m_inletVel[0] + m_inletVel[1]) - Ny;

  m_distributions[8] = m_distributions[13]  + (1 / 6.0) * rhoN * (m_inletVel[0] + m_inletVel[2]) - Nz;

  m_distributions[9] = m_distributions[12]  + (1 / 6.0) * rhoN * (m_inletVel[0] - m_inletVel[2]) + Nz;

 // MacroscopicValue(velSet);
}

#endif





