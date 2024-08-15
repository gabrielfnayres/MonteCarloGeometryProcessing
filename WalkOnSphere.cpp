#include "geometry.h"


double Geometry::cross(const Vec2D &a, const Vec2D &b){
    return ((a.real()*b.real()) - (a.imag()*b.imag()));
}

double Geometry::dot(const Vec2D& a, const Vec2D& b){
    return ( a.real()*b.real() + a.imag()*b.imag());
}

double Geometry::length(const Vec2D& a){
    return sqrt(norm(a));
}

Geometry::Vec2D Geometry::closestPoint(Vec2D x, Vec2D a, Vec2D b){ // atomic closest
    Vec2D u = b - a;
    double t = std::clamp( dot(x - a, u)/dot(u,u), 0.0, 1.0);
    return(1.0 - t)*a.real() + t*b.real();
}

double Geometry::distancePolylines(Vec2D &x, const std::vector<Polyline>& P){

    double d = INFINITY;

    for(int i = 0; i < P.size(); i++){
        for(int j = 0; j < P[i].size(); j++){
            Vec2D y = closestPoint(x, P[i][j], P[i][j+1]);
            d = std::min(d, length(y));
        }
    }
    return d;
}

bool Geometry::isSilhouette(Vec2D x, Vec2D a, Vec2D b, Vec2D c){
    return cross(b-a, x-a) * cross(c-b, x-b) < 0;
}

double Geometry::silhouetteDistancePolylines(Vec2D x, const std::vector<Polyline>& P){
    double d =  INFINITY;
    for(int i = 0; i < P.size(); i++){
        for(int j = 0; j < P[i].size(); j++){
            if(isSilhouette(x, P[i][j - 1], P[i][j], P[i][j + 1])){
                Vec2D t = x - P[i][j];
                d = std::min(d, length(t));
            }
        }
    }
    return d;
}

double Geometry::rayIntersection(Vec2D x, Vec2D v, Vec2D a, Vec2D b){
    Vec2D u = b - a;
    Vec2D w = x - a;

    double d = cross(v, u);
    double s = cross(v, w) / d;
    double t = cross(u, w) / d;

    if(t > 0. && 0. <= s && s <= 1.) return t;
    return INFINITY;
}

