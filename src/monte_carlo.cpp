#include "../includes/geometry.h"

std::uniform_real_distribution<double> uniformDis(0.0, 1.0);

std::mt19937 rnd;
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


double Geometry::calculateSphereArea( double r){
    return  M_PI*(r*r);
}

Geometry::Vec2D Geometry::randPointOnSphere(Vec2D &x, double r){
    double theta = uniformDis(rnd)*2*M_PI;
    return Vec2D{x.real() + r*cos(theta), x.imag() + r*sin(theta)};
}

Geometry::Vec2D Geometry::randPointInSphere(Vec2D &x, double r){
    double radious = r*sqrt(uniformDis(rnd));
    double theta  = uniformDis(rnd)*2*M_PI;
    
    return Vec2D{x.real() + radious + cos(theta), x.imag() + radious + sin(theta)};
}

double Geometry::G(Vec2D &x, Vec2D &y, double sphereR){

    double r = abs((y.real() - x.real()) + (y.imag() - x.imag()));
    return ((1/2*(M_PI))*log(sphereR/r));
}

double Geometry::gradientG(Vec2D &x, Vec2D &y, double sphereR){
    double r = abs((y.real() - x.real()) + (y.imag() - x.imag()));
    double n = (y.real() - x.real()) + (y.real() - x.real());
    return ((n/2*M_PI)*((1/(r*r)) - (1/(sphereR*sphereR))));
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

double Geometry::WalkOnSphere(Vec2D x0, std::vector<Polyline> boudaryDirichilet, std::function<double(Vec2D)> g){
    
    double epsilon = 1e-4;
    int nWalks = 65523;
    int maxSteps = 65523;
    double rmin =  1e-4;

    double sum = 0.0;

    for(int i = 0; i < nWalks; i++){

        Vec2D x = x0;
        bool onBoundary = false;
        double Dirichilet;
        double Silhouette;
        int steps = 0;

        do{
            Dirichilet = distancePolylines(x, boudaryDirichilet);    
            

            double r =  std::max(rmin, Dirichilet);
            Silhouette = calculateSphereArea(r);
            
            Vec2D point = randPointOnSphere(x, r);

            steps++;
        }while(Dirichilet > epsilon && steps >= maxSteps);

        if(steps >= maxSteps){
            std::cerr << "MAX STEPS!!!!!!" << std::endl;
        }
        sum += g(x);
    }
    return sum/nWalks;
}

double Geometry::signedAngle(Vec2D x, const std::vector<Polyline> &P){
    
    double theta = 0;
    for(int i = 0; i < P.size(); i++){
        for(int j = 0; j < P[i].size(); j++){
            theta += std::arg((P[i][j+1] - x)/(P[i][j] - x));
        }
    }

    return theta;
}

bool Geometry::insideDomain(Vec2D x, const std::vector<Polyline> &boundaryDirichilet, const std::vector<Polyline>& boundaryNeumann){

    double theta = signedAngle(x, boundaryDirichilet) + signedAngle(x, boundaryNeumann);

    const double delta = 1e-4;
    return std::abs(theta-2*M_PI) < delta;
}
