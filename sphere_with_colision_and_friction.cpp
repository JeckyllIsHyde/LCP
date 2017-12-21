#include <iostream>

#include "LittleMath.h"
#include "LCP.h"

const double dt = 0.0025;
const double tmax = 1.5;

const double MASS = 1.0, J = 0.002;
const double TH0 = 30.*M_PI/180;
const double A_GRAVITY = 10.0;
const double K_FLOOR = 1000;
const double KD_FLOOR = 50;

struct RigidBody2D {
  double mMass, mInertia;
  SpatialVector2D v, a, f;
  SpatialTransform2D X;

  RigidBody2D(): mMass(1.0), mInertia(1.0) {}
  RigidBody2D(double m, double I): mMass(m), mInertia(I) {}
  void forwardDynamics();
  void step( double dt );
};

struct PhysicsEngine {
  RigidBody2D body;
  
  void update( double dt );
  void printData( double t );

  void viscoelasticCollision();
  void lcpCollision( double dt );
};

void init_engine_simple_sphere( PhysicsEngine& engine );

int main() {

  PhysicsEngine engine;

  // init engine
  init_engine_simple_sphere( engine );

  // simulate
  double t;
  for ( t=0; t<=tmax+dt; t+=dt ) {
    engine.printData( t );
    engine.update( dt );
  }
  
  return 0;
}

void init_engine_simple_sphere( PhysicsEngine& engine ) {

  engine.body = RigidBody2D(MASS,J);

  engine.body.X.setTheta(TH0); // angular position 
  engine.body.X.r(0) = 0.0; // x-axis position 
  engine.body.X.r(1) = 1.0; // y-axis position

  engine.body.v(0) = 0.0; // angular velocity 
  engine.body.v(1) = 0.5; // x-axis linear velocity 
  engine.body.v(2) = 0.0; // y-axis linear velocity
}

void PhysicsEngine::update( double dt ) {
  body.f = SpatialVector2D();
  SpatialVector2D g_force(0.0,0.0,-MASS*A_GRAVITY);
  body.f+=g_force;
  viscoelasticCollision();
  lcpCollision( dt );
  body.forwardDynamics();
  body.step( dt );
}

void PhysicsEngine::viscoelasticCollision() {
  if ( body.X.r(1)>0 )
    body.f+= SpatialVector2D();
  else
    body.f+= SpatialVector2D(0.0,
			     0.0,
			     -K_FLOOR*body.X.r(1)
			     -KD_FLOOR*body.v(2));
}

void PhysicsEngine::lcpCollision( double dt ) {
  MatrixT<double,3,3> M_;
  M_(0,0) = 1./body.mInertia; M_(1,1) = M_(2,2) = 1./body.mMass;
  MatrixT<double,3,2> D;
  D(0,0) = D(0,1) = 0.0;
  D(1,0) = 1.0; D(1,1) =-1.0;
  D(2,0) = D(2,1) = 0.0;
  MatrixT<double,3,1> N;
  N(0,0) = 0.0;
  N(1,0) = 1.0;
  N(2,0) = 0.0;
  MatrixT<double,2,1> e;
  e(0,0) = e(1,0) = 1.0;
  MatrixT<double,1,1> z;
  z(0,0) = 0.0;
  MatrixT<double,1,1> mu;
  mu(0,0) = 0.6;
  MatrixT<double,3,1> V(body.v);
  MatrixT<double,3,1> q;
  q(0,0)=body.X.getTheta();q(1,0)=body.X.r(0);q(2,0)=body.X.r(1);

  MatrixT<double,4,4> M( ( (D.transpose()*M_*D).
			   cath((D.transpose()*M_*N)).cath(e)).catv
			 ( (N.transpose()*M_*D).
			   cath(N.transpose()*M_*N).cath(z)).catv
			 ( (e.transpose()*(-1.0)).cath(mu).cath(z)));
  MatrixT<double,4,1> b( ( D.transpose()*V ).catv
			 ( N.transpose()*(q/dt+V)).catv(z));

  LCP contact_problem;
  contact_problem.read_LCP( b,M );
  contact_problem.is_verbose = true;
  contact_problem.lemke_algorithm();
  /*
  std::cout << "M_:\n" << M_ << std::endl;
  std::cout << "D:\n" << D << std::endl;
  std::cout << "N:\n" << N << std::endl;
  std::cout << "e:\n" << e << std::endl;
  std::cout << "[D.T*M_*D | D.T*M_*N | e]\n"
	    << "[N.T*M_*D | N.T*M_*N | 0]:\n"
	    << "[-e.T     | mu       | 0]\n"
	    << M
	    << std::endl
	    << "[D.T*V   ]\n[N.T*(q/dt+V)]:\n[0]\n"
	    << b
	    << std::endl;
  */
}

void PhysicsEngine::printData( double t ) {
  const char* separator = ", ";
  std::cout << t << separator
	    << body.X.getTheta()*180./M_PI << separator
	    << body.X.r(0) << separator
	    << body.X.r(1)
	    << std::endl;
}

void RigidBody2D::forwardDynamics() {
  // a = FD( q,v,f )
  a(0) = f(0)/mInertia;
  a(1) = f(1)/mMass;
  a(2) = f(2)/mMass;
}

void RigidBody2D::step( double dt ) {
  //  q += dt*v;
  double qy = X.getTheta();
  qy += dt*v(0); // angular thz
  X.setTheta( qy );
  X.r(0) += dt*v(1); // position x
  X.r(1) += dt*v(2); // position y
  //  v += dt*a;
  v(0) += dt*a(0); // angular velocity wz
  v(1) += dt*a(1); // linear velocity vx
  v(2) += dt*a(2); // linear velocity vy
}
