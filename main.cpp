#include <geometry.h>

double lines(Geometry::Vec2D x){
    const double s = 8.0;

    return std::fmod(std::floor(x.real())*s, 2.0);
}

std::vector<Geometry::Polyline> boundaryDirichlet = {
   {{ Geometry::Vec2D(0.2, 0.2), Geometry::Vec2D(0.6, 0.0), Geometry::Vec2D(1.0, 0.2) }},
   {{ Geometry::Vec2D(1.0, 1.0), Geometry::Vec2D(0.6, 0.8), Geometry::Vec2D(0.2, 1.0) }}
};
std::vector<Geometry::Polyline> boundaryNeumann = {
   {{ Geometry::Vec2D(1.0, 0.2), Geometry::Vec2D(0.8, 0.6), Geometry::Vec2D(1.0, 1.0) }},
   {{ Geometry::Vec2D(0.2, 1.0), Geometry::Vec2D(0.0, 0.6), Geometry::Vec2D(0.2, 0.2) }}
};


int main(){

    srand(time(NULL));

    Geometry geo;

    std::ofstream out( "out.csv" );

   int s = 128; // image size
   for( int j = 0; j < s; j++ )
   {
      std::cerr << "row " << j << " of " << s << std::endl;
      for( int i = 0; i < s; i++ )
      {
         Geometry::Vec2D x0( ((double)i+.5)/((double)s),
                   ((double)j+.5)/((double)s) );
         double u = 0.;
         if( geo.insideDomain(x0, boundaryDirichlet, boundaryNeumann) )
            u = geo.WalkOnStars( x0, boundaryDirichlet, boundaryNeumann, lines);
         out << u;
         if( i < s-1 ) out << ",";
      }
      out << std::endl;
   }

}