#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <algorithm>
#include <iostream>
#include <array>
#include <complex>
#include <random>
#include <functional>
#include <fstream>
#include <vector>
#include <random>
#include <cstdlib>
#include <Eigen/Geometry>
#include <numeric>


class Geometry{

public:

    using Vec2D = std::complex<double>;
    using Polyline = std::vector<Vec2D>;

    Vec2D rotate90(Vec2D a);

    double angleOf(Vec2D a);

    double randomG(double a, double b);

    double cross(const Vec2D &a, const Vec2D &b);

    double dot(const Vec2D& a, const Vec2D& b);

    double length(const Vec2D& a);

    Vec2D closestPoint( Vec2D x, Vec2D a, Vec2D b);

    double distancePolylines(Vec2D& x, const std::vector<Polyline>& P);

    bool isSilhouette(Vec2D x,  Vec2D a, Vec2D b, Vec2D c);

    double silhouetteDistancePolylines(Vec2D x, const std::vector<Polyline>& P);

    double rayIntersection(Vec2D x,Vec2D v,  Vec2D a,Vec2D b);

    Vec2D intersectPolylines(Vec2D x, Vec2D v, double r, const std::vector<Polyline> &P, Vec2D& n, bool &onBoundary);

    double WalkOnStars(Vec2D x0, std::vector<Polyline> boundaryDirichilet, std::vector<Polyline> boundaryNeumann, std::function<double(Vec2D)> g);

    double WalkOnSphere(Vec2D x0, std::vector<Polyline> boundaryDirichilet, std::function<double(Vec2D)> g);

};


#endif