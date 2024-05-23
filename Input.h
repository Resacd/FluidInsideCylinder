#pragma once
#ifndef _INPUT_H_
#define _INPUT_H_

#include <vector>




double dx = 0.1;
double u_max = 0.1;
//u_max*dt/dx<<1
double Ma = 0.1;
double dt = Ma * (dx / u_max);
double nu = 1e-6;
//double tau = (3 * nu * dt) / (dx * dx)+0.5;
const double tau = 0.67;
int mdepth = 1;
double pressure = 1.0;


#endif 
