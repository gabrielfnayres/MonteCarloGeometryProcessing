#include <algorithm>
#include <array>
#include <complex>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <limits>


using Vec2D = std::complex<double>;
using Polyline = std::vector<Vec2D>;
const double infinity = std::numeric_limits<double>::infinity();

std::uniform_real_distribution<double> uniformDis(0.0, 1.0);
std::default_random_engine rnd;

double randomG(double a, double b){
    const double rmax = 1.0/(double)RAND_MAX;
    double u = rmax*(double)rand();
    return u*(b - a) + a;
}

double angleOf(Vec2D a){
    return arg(a);
}

double cross( Vec2D a,  Vec2D b){
    return ((real(a)*real(b)) - (imag(a)*imag(b)));
}

double dot(Vec2D a, Vec2D b){
    return ( real(a)*real(b) + imag(a)*imag(b));
}

double length( Vec2D a){
    return sqrt(norm(a));
}

Vec2D rotate90(Vec2D a){
    return Vec2D(-imag(a), real(a));
}

Vec2D closestPoint(Vec2D x, Vec2D a, Vec2D b){ // atomic closest
    Vec2D u = b - a;
    double t = std::clamp( dot(x - a, u)/dot(u,u), 0.0, 1.0);
    return(1.0 - t)*a+ t*b;
}


double calculateSphereArea( double r){
    return  M_PI*(r*r);
}

Vec2D randPointOnSphere(Vec2D &x, double r){
    double theta = uniformDis(rnd)*2*M_PI;
    return Vec2D{real(x) + r*cos(theta), imag(x) + r*sin(theta)};
}

Vec2D randPointInSphere(Vec2D &x, double r){
    double radious = r*sqrt(uniformDis(rnd));
    double theta  = uniformDis(rnd)*2*M_PI;
    
    return Vec2D{x.real() + radious + cos(theta), x.imag() + radious + sin(theta)};
}

double G(Vec2D &x, Vec2D &y, double sphereR){

    double r = abs((std::real(y) - x.real()) + (y.imag() - x.imag()));
    return ((1/2*(M_PI))*log(sphereR/r));
}

double gradientG(Vec2D &x, Vec2D &y, double sphereR){
    double r = abs((std::real(y) - std::real(x)) + (std::imag(y) - std::imag(x)));
    double n = (std::real(y) - std::real(x)) + (std::imag(y) - std::imag(x));
    return ((n/2*M_PI)*((1/(r*r)) - (1/(sphereR*sphereR))));
}

double distancePolylines(Vec2D x, const std::vector<Polyline>& P){

    double d = infinity;

    for(int i = 0; i < P.size(); i++){
        for(int j = 0; j < P[i].size() - 1; j++){
            Vec2D y = closestPoint(x, P[i][j], P[i][j+1]);
            d = std::min(d, length(x - y));
        }
    }
    return d;
}

bool isSilhouette(Vec2D x, Vec2D a, Vec2D b, Vec2D c){
    return cross(b-a, x-a) * cross(c-b, x-b) < 0;
}

double silhouetteDistancePolylines(Vec2D x, const std::vector<Polyline>& P){
    double d =  infinity;
    for(int i = 0; i < P.size(); i++){
        for(int j = 0; j < P[i].size() - 1; j++){
            if(isSilhouette(x, P[i][j - 1], P[i][j], P[i][j + 1])){
                d = std::min(d, length(x - P[i][j]));
            }
        }
    }
    return d;
}

double rayIntersection(Vec2D x, Vec2D v, Vec2D a, Vec2D b){
    Vec2D u = b - a;
    Vec2D w = x - a;

    double d = cross(v, u);
    double s = cross(v, w) / d;
    double t = cross(u, w) / d;

    if(t > 0. && 0. <= s && s <= 1.)
        return t;
    return infinity;
}

Vec2D intersectPolylines(Vec2D x, Vec2D v, double r, const std::vector<Polyline> &P, Vec2D& n, bool &onBoundary){

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


double WalkOnStars(Vec2D x0, std::vector<Polyline>boundaryDirichilet, std::vector<Polyline> boundaryNeumann, std::function<double(Vec2D)> g){

    const double epsilon = 0.0001;
    const double rmin = 0.0001; // limits how small the steps will shrink near the silhouette
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

            Dirichilet = distancePolylines(x,boundaryDirichilet);
            Silhouette = silhouetteDistancePolylines(x, boundaryNeumann);
            r = std::max(rmin, std::min(Dirichilet, Silhouette));

            double theta = randomG(-M_PI, M_PI);
            
            if(onBoundary){
                theta = theta/2. + angleOf(n);
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

double WalkOnSphere(Vec2D x0, std::vector<Polyline> boudaryDirichlet, std::function<double(Vec2D)> g){
    
    const double epsilon = 0.0001;
    const int nWalks = 65523;
    const int maxSteps = 65523;
    const double rmin =  0.0001;

    double sum = 0.0;

    for(int i = 0; i < nWalks; i++){

        Vec2D x = x0;
        bool onBoundary = false;
        double Dirichilet;
        double Silhouette;
        int steps = 0;

        do{
            Dirichilet = distancePolylines(x, boudaryDirichlet);    
            
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

double signedAngle(Vec2D x, const std::vector<Polyline> &P){
    
    double theta = 0;
    for(int i = 0; i < P.size(); i++){
        for(int j = 0; j < P[i].size(); j++){
            theta += arg((P[i][j+1] - x)/(P[i][j] - x));
        }
    }

    return theta;
}

bool insideDomain(Vec2D x, const std::vector<Polyline> &boundaryDirichlet, const std::vector<Polyline>& boundaryNeumann){

    double theta = signedAngle(x,boundaryDirichlet) + signedAngle(x, boundaryNeumann);

    const double delta = 0.0001;

    //std::cerr << "Point (" << x.real() << ", " << x.imag() << ") is " << (inside ? "inside" : "outside") << " the domain" << std::endl;

    return std::abs(theta - 2 * M_PI) < delta;;
}

double lines(Vec2D x){
    const double s = 8.0;

    return fmod(floor(real(x))*s, 2.0);
}


std::vector<Polyline> boundaryDirichlet = {
   {{ Vec2D(0.2, 0.2), Vec2D(0.6, 0.0), Vec2D(1.0, 0.2) }},
   {{ Vec2D(1.0, 1.0), Vec2D(0.6, 0.8), Vec2D(0.2, 1.0) }}
};
std::vector<Polyline> boundaryNeumann = {
   {{ Vec2D(1.0, 0.2), Vec2D(0.8, 0.6), Vec2D(1.0, 1.0) }},
   {{ Vec2D(0.2, 1.0), Vec2D(0.0, 0.6), Vec2D(0.2, 0.2) }}
};


int main(){

   srand(time(NULL));

   std::ofstream out( "result.csv" );
 
   int s = 128; // image size
   for( int j = 0; j < s; j++ )
   {
      std::cerr << "row " << j << " of " << s << std::endl;
      for( int i = 0; i < s; i++ )
      {
         Vec2D x0( ((double)i+.5)/((double)s),
                   ((double)j+.5)/((double)s) );
         double u = 0.;
         if( insideDomain(x0, boundaryDirichlet, boundaryNeumann) )
            u = WalkOnStars( x0, boundaryDirichlet, boundaryNeumann, lines);
         out << u;
         if( i < s-1 ) out << ",";
      }
      out << std::endl;
   }
   return 0;

}