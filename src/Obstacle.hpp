#ifndef OBSTACLE_HPP
#define OBSTACLE_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"

void AddFlatWallNewBC(particles &vd,
					  const int k0,
					  const int kmax,
					  const Point<DIM, double> Corner,
					  const Point<DIM, double> UnitOffset,
					  const double dx,
					  const Point<DIM, double> obstacle_centre,
					  const Point<DIM, double> obstacle_velocity,
					  const Parameters &arg_p,
					  const double obstacle_omega = 0.0);

void AddFlatWallModNewBC(particles &vd,
						 const int k0,
						 const int kmax,
						 const Point<DIM, double> Corner,
						 const Point<DIM, double> UnitOffset,
						 const double dx,
						 const Point<DIM, double> obstacle_centre,
						 const Point<DIM, double> obstacle_velocity,
						 const Point<DIM, double> given_normal,
						 const double obstacle_omega);

bool isAvobeLine(Point<DIM, double> P, Point<DIM, double> Q, Point<DIM, double> EvalPoint, double dp);
bool isBelowLine(Point<DIM, double> P, Point<DIM, double> Q, Point<DIM, double> EvalPoint, double dp);

class Obstacle
{
public:
	Point<DIM, double> Centre_;
	const Parameters &params_;
	Point<DIM, double> LinearVelocity_;
	double AngularVelocity_;
	double refine_factor;

public:
	Obstacle(Point<DIM, double> Centre,
			 const Parameters &p,
			 Point<DIM, double> vel = {0.0, 0.0},
			 double omega = 0.0,
			 double rf = 1.0);

	virtual ~Obstacle() = default;
	Obstacle(const Parameters &p);
	virtual bool isInside(Point<DIM, double> P) = 0;
	virtual void AddObstacle(particles &vd) = 0;
};

class EmptyObstacle : public Obstacle
{
public:
	EmptyObstacle(const Parameters &p);
	~EmptyObstacle() override = default;
	bool isInside(Point<DIM, double> P) override;
	void AddObstacle(particles &vd) override;
};

class CylinderObstacle : public Obstacle
{
private:
	double Radius_;
	Sphere<DIM, double> Cylinder_;

public:
	CylinderObstacle(double Radius,
					 Point<DIM, double> centre,
					 const Parameters &p,
					 Point<DIM, double> vel = {0.0, 0.0},
					 double omega = 0.0,
					 double rf = 1.0);
	~CylinderObstacle() override = default;
	bool isInside(Point<DIM, double> P) override;
	bool isOutside(Point<DIM, double> P); // for outer cylinder in taylor couette
	void AddObstacle(particles &vd);
};

class EllipticObstacle : public Obstacle
{
private:
	double Major_;
	double Minor_;
	double tilt_;

public:
	EllipticObstacle(double Major_,
					 double Minor_,
					 double tilt_,
					 Point<DIM, double> centre,
					 const Parameters &p,
					 Point<DIM, double> vel = {0.0, 0.0},
					 double omega = 0.0,
					 double rf = 1.0);
	~EllipticObstacle() override = default;
	bool isInside(Point<DIM, double> P) override;
	void AddObstacle(particles &vd);
	Point<DIM, double> parametricEllipse(double theta);
};

class RectangleObstacle : public Obstacle
{
private:
	const double BaseLength_;
	const double HeigthLength_;
	Box<DIM, double> Rectangle_;
	Point<DIM, double> LowerLeft_;
	// Point<DIM, double> LowerRight_;
	// Point<DIM, double> UpperLeft_;
	Point<DIM, double> UpperRight_;

public:
	RectangleObstacle(Point<DIM, double> centre,
					  const Parameters &p,
					  double BaseLength,
					  double HeigthLength,
					  Point<DIM, double> vel = {0.0, 0.0},
					  double omega = 0.0,
					  double rf = 1.0);
	~RectangleObstacle() override = default;
	bool isInside(Point<DIM, double> P) override;
	void AddObstacle(particles &vd);
};

class TriangleObstacle : public Obstacle
{
private:
	const double BaseLength_;
	const double HeigthLength_;
	Box<DIM, double> ContainingRectangle_;
	Point<DIM, double> LowerLeft_;
	Point<DIM, double> LowerRight_;
	Point<DIM, double> UpperRight_;

public:
	TriangleObstacle(Point<DIM, double> centre,
					 const Parameters &p,
					 double BaseLength,
					 double HeigthLength,
					 Point<DIM, double> vel = {0.0, 0.0},
					 double omega = 0.0,
					 double rf = 1.0);
	~TriangleObstacle() override = default;
	bool isInside(Point<DIM, double> P) override;
	void AddObstacle(particles &vd);
};

class TriangleEqui : public Obstacle
{
private:
	const double SideLength_;
	Box<DIM, double> ContainingRectangle_;
	Point<DIM, double> LowerLeft_;
	Point<DIM, double> UpperRight_;
	Point<DIM, double> LowerRight_;
	Point<DIM, double> TriangleTip_;

public:
	TriangleEqui(Point<DIM, double> centre,
				 const Parameters &p,
				 double sidelength,
				 Point<DIM, double> vel = {0.0, 0.0},
				 double omega = 0.0,
				 double rf = 1.0);
	~TriangleEqui() override = default;
	bool isInside(Point<DIM, double> P) override;
	void AddObstacle(particles &vd);
};

#endif // OBSTACLE_HPP
