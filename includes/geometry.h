
#include <algorithm>
#include <iostream>
#include <array>
#include <complex>
#include <cmath>
#include <random>
#include <functional>
#include <fstream>
#include <vector>
#include <random>
#include <cstdlib>
#include <numeric>
#include <limits>

using Vec2D = std::complex<double>;
using Polyline = std::vector<Vec2D>;
const double infinity = std::numeric_limits<double>::infinity();

Vec2D rotate90(Vec2D a);


double lines(Vec2D x);

double angleOf(Vec2D a);

double randomG(double a, double b);

double cross( Vec2D a, Vec2D b);

double dot(Vec2D a, Vec2D b);

double length(Vec2D a);

Vec2D closestPoint( Vec2D x, Vec2D a, Vec2D b);

double distancePolylines(Vec2D x, const std::vector<Polyline>& P);

bool isSilhouette(Vec2D x,  Vec2D a, Vec2D b, Vec2D c);

double silhouetteDistancePolylines(Vec2D x, const std::vector<Polyline>& P);

double calculateSphereArea(double r);

Vec2D randPointOnSphere(Vec2D &x, double r);

double G(Vec2D &x, Vec2D &y, double sphereR);

double gradientG(Vec2D &x, Vec2D &y, double sphereR);

Vec2D randPointInSphere(Vec2D &x, double r);

double rayIntersection(Vec2D x,Vec2D v,  Vec2D a,Vec2D b);

Vec2D intersectPolylines(Vec2D x, Vec2D v, double r, const std::vector<Polyline> &P, Vec2D& n, bool &onBoundary);

double WalkOnStars(Vec2D x0, std::vector<Polyline> boundaryDirichlet, std::vector<Polyline> boundaryNeumann, std::function<double(Vec2D)> g);

double WalkOnSphere(Vec2D x0, std::vector<Polyline> boundaryDirichlet, std::function<double(Vec2D)> g);

double signedAngle(Vec2D x, const std::vector<Polyline>& P);

bool insideDomain(Vec2D x, const std::vector<Polyline> &boundaryDirichlet, const std::vector<Polyline>& boundaryNeumann);

