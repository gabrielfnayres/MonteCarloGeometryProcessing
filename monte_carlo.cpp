#include "geometry.h"


double Geometry::randomG(double a, double b){
    double n = ((double)rand() / (a - b + 1)) + b;
    return n;
}

double Geometry::angleOf(Vec2D a){
    return atan(a.imag()/a.real());
}

double Geometry::cross(const Vec2D &a, const Vec2D &b){
    return ((a.real()*b.real()) - (a.imag()*b.imag()));
}

double Geometry::dot(const Vec2D& a, const Vec2D& b){
    return ( a.real()*b.real() + a.imag()*b.imag());
}

double Geometry::length(const Vec2D& a){
    return sqrt(norm(a));
}

Geometry::Vec2D Geometry::rotate90(Vec2D a){
    return Vec2D(-a.imag(), a.real());
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

Geometry::Vec2D Geometry::intersectPolylines(Vec2D x, Vec2D v, double r, const std::vector<Polyline> &P, Vec2D& n, bool &onBoundary){


    double tm = r;
    n = Vec2D{0.0, 0.0};
    onBoundary = false;
    
    for(int i = 0; i < P.size(); i++){
        for(int j = 0; j < P[i].size() - 1; j++){
            const double c = 1e-5;
            double t = rayIntersection(x + c*v, v, P[i][j], P[i][j+1]);

            if(t < tm){
                tm = t;
                n = rotate90(P[i][j+1] - P[i][j]);
                n /= length(n);
                onBoundary = true;
            }
        }
    }
    return x + tm*v;
}

double Geometry::WalkOnStars(Vec2D x0, std::vector<Polyline> boundaryDirichilet, std::vector<Polyline> boundaryNeumann, std::function<double(Vec2D)> g){

    const double epsilon = 1e-4;
    const double rmin = 1e-4; // limits how small the steps will shrink near the silhouette
    const int nWalks = 65536;
    const int maxSteps = 65536;

    double sum = 0.0; // accumulate all values g we encounter at the boundary


    for(int i = 0; i < nWalks; i++){
        Vec2D x = x0;
        Vec2D n{0.0, 0.0};

        bool onBoundary = false;
        double r;
        double Dirichilet;
        double Silhouette;
        int steps = 0;

        do{

            Dirichilet = distancePolylines(x, boundaryDirichilet);
            Silhouette = silhouetteDistancePolylines(x, boundaryNeumann);
            r = std::max(rmin, std::min(Dirichilet, Silhouette));

            double theta = randomG(-M_PI, M_PI);
            
            if(onBoundary){
                theta = theta/2 + angleOf(n);
            }

            Vec2D v{cos(theta), sin(theta)};

            x = intersectPolylines(x, v, r, boundaryNeumann, n, onBoundary);

            steps++;
        } while(Dirichilet > epsilon && steps < maxSteps);

        if(steps >= maxSteps){
            std::cerr << "MAX STEPS!!!!!!" << std::endl;
        }

        sum += g(x);
    }
    return sum/nWalks;
}

