#ifndef CUSTOM_DIST_H
#define CUSTOM_DIST_H

/// Sample an element from a gaussian distribution using the Box-Muller Method
/// @param[in] std_deviation Standard deviation
/// @return A sample from a gaussian distribution using the box-muller transform
double dist_gauss_box_muller(double std_deviation);

/// Sample an element from a gaussian distribution using the Polar Method
/// @param[in] std_deviation Standard deviation
/// @return A sample from a gaussian distribution using the polar method
double dist_gauss_polar(double std_deviation);

/// Sample an element from a gaussian distribution using the Ziggurat Algorithm
/// @param[in] std_deviation Standard deviation
/// @return A sample from a gaussian distribution using the ziggurat algorithm
double dist_gauss_ziggurat(double std_deviation);

#endif //CUSTOM_DIST_H
