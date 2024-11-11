#ifndef OBSTACLE_HPP
#define OBSTACLE_HPP

#include "Definitions.hpp"
#include "VectorUtilities.hpp"

void AddFlatWallNewBC(particles &vd,
					  const int k0,
					  const int kmax,
					  const Point<DIM, real_number> Corner,
					  const Point<DIM, real_number> UnitOffset,
					  const real_number dx,
					  const Point<DIM, real_number> obstacle_centre,
					  const Point<DIM, real_number> obstacle_velocity,
					  const Parameters &arg_p,
					  const size_t particle_type,
					  const real_number obstacle_omega = 0.0);

void AddFlatWallModNewBC(particles &vd,
						 const int k0,
						 const int kmax,
						 const Point<DIM, real_number> Corner,
						 const Point<DIM, real_number> UnitOffset,
						 const real_number dx,
						 const Point<DIM, real_number> obstacle_centre,
						 const Point<DIM, real_number> obstacle_velocity,
						 const Parameters &arg_p,
						 const size_t particle_type,
						 const Point<DIM, real_number> given_normal,
						 const real_number given_kappa = 0.0,
						 const real_number obstacle_omega = 0.0);

bool isAvobeLine(Point<DIM, real_number> P, Point<DIM, real_number> Q, Point<DIM, real_number> EvalPoint, real_number dp);
bool isBelowLine(Point<DIM, real_number> P, Point<DIM, real_number> Q, Point<DIM, real_number> EvalPoint, real_number dp);

class Obstacle
{
public:
	Point<DIM, real_number> Centre_;
	const Parameters &params_;
	Point<DIM, real_number> LinearVelocity_;
	real_number AngularVelocity_;
	real_number refine_factor;

public:
	Obstacle(Point<DIM, real_number> Centre,
			 const Parameters &p,
			 Point<DIM, real_number> vel = {0.0, 0.0},
			 real_number omega = 0.0,
			 real_number rf = 1.0);

	virtual ~Obstacle() = default;
	Obstacle(const Parameters &p);
	virtual bool isInside(Point<DIM, real_number> P) = 0;
	virtual void AddObstacle(particles &vd) = 0;
	virtual bool isInsidePlusTol(Point<DIM, real_number> P) = 0;
};

class EmptyObstacle : public Obstacle
{
public:
	EmptyObstacle(const Parameters &p);
	~EmptyObstacle() override = default;
	bool isInside(Point<DIM, real_number> P) override;
	bool isInsidePlusTol(Point<DIM, real_number> P);
	void AddObstacle(particles &vd) override;
};

class CylinderObstacle : public Obstacle
{
private:
	real_number Radius_;

public:
	CylinderObstacle(real_number Radius,
					 Point<DIM, real_number> centre,
					 const Parameters &p,
					 Point<DIM, real_number> vel = {0.0, 0.0},
					 real_number omega = 0.0,
					 real_number rf = 1.0);
	~CylinderObstacle() override = default;
	bool isInside(Point<DIM, real_number> P) override;
	bool isInsidePlusTol(Point<DIM, real_number> P) override;
	bool isInside_minEps(Point<DIM, real_number> P); // for outer cylinder in taylor couette
	void AddObstacle(particles &vd);
};

class SphereObstacle : public Obstacle
{
private:
	real_number Radius_;
	bool autoNormals_;

public:
	SphereObstacle(real_number Radius,
				   Point<DIM, real_number> centre,
				   const Parameters &p,
				   Point<DIM, real_number> vel = {0.0, 0.0, 0.0},
				   real_number omega = 0.0,
				   real_number rf = 1.0,
				   bool autoNormals = false);
	~SphereObstacle() override = default;
	bool isInside(Point<DIM, real_number> P) override;
	bool isInsidePlusTol(Point<DIM, real_number> P);
	void AddObstacle(particles &vd);
};

class EllipticObstacle : public Obstacle
{
private:
	real_number Major_;
	real_number Minor_;
	real_number tilt_;

public:
	EllipticObstacle(real_number Major_,
					 real_number Minor_,
					 real_number tilt_,
					 Point<DIM, real_number> centre,
					 const Parameters &p,
					 Point<DIM, real_number> vel = {0.0, 0.0},
					 real_number omega = 0.0,
					 real_number rf = 1.0);
	~EllipticObstacle() override = default;
	bool isInside(Point<DIM, real_number> P) override;
	bool isInsidePlusTol(Point<DIM, real_number> P);
	void AddObstacle(particles &vd);
	Point<DIM, real_number> parametricEllipse(real_number theta);
};

class RectangleObstacle : public Obstacle
{
private:
	const real_number BaseLength_;
	const real_number HeigthLength_;
	Box<DIM, real_number> Rectangle_;
	Point<DIM, real_number> LowerLeft_;
	// Point<DIM, real_number> LowerRight_;
	// Point<DIM, real_number> UpperLeft_;
	Point<DIM, real_number> UpperRight_;

public:
	RectangleObstacle(Point<DIM, real_number> centre,
					  const Parameters &p,
					  real_number BaseLength,
					  real_number HeigthLength,
					  Point<DIM, real_number> vel = {0.0, 0.0},
					  real_number omega = 0.0,
					  real_number rf = 1.0);
	~RectangleObstacle() override = default;
	bool isInside(Point<DIM, real_number> P) override;
	bool isInsidePlusTol(Point<DIM, real_number> P);
	void AddObstacle(particles &vd);
};

class TriangleObstacle : public Obstacle
{
private:
	const real_number BaseLength_;
	const real_number HeigthLength_;
	Box<DIM, real_number> ContainingRectangle_;
	Point<DIM, real_number> LowerLeft_;
	Point<DIM, real_number> LowerRight_;
	Point<DIM, real_number> UpperRight_;

public:
	TriangleObstacle(Point<DIM, real_number> centre,
					 const Parameters &p,
					 real_number BaseLength,
					 real_number HeigthLength,
					 Point<DIM, real_number> vel = {0.0, 0.0},
					 real_number omega = 0.0,
					 real_number rf = 1.0);
	~TriangleObstacle() override = default;
	bool isInside(Point<DIM, real_number> P) override;
	bool isInsidePlusTol(Point<DIM, real_number> P);
	void AddObstacle(particles &vd);
};

class TriangleEqui : public Obstacle
{
private:
	const real_number SideLength_;
	Box<DIM, real_number> ContainingRectangle_;
	Point<DIM, real_number> LowerLeft_;
	Point<DIM, real_number> UpperRight_;
	Point<DIM, real_number> LowerRight_;
	Point<DIM, real_number> TriangleTip_;

public:
	TriangleEqui(Point<DIM, real_number> centre,
				 const Parameters &p,
				 real_number sidelength,
				 Point<DIM, real_number> vel = {0.0, 0.0},
				 real_number omega = 0.0,
				 real_number rf = 1.0);
	~TriangleEqui() override = default;
	bool isInside(Point<DIM, real_number> P) override;
	bool isInsidePlusTol(Point<DIM, real_number> P);
	void AddObstacle(particles &vd);
};

class TriangleTestObstacle : public Obstacle
{
private:
	const real_number BaseLength_;
	const real_number HeigthLength_;
	Box<DIM, real_number> ContainingRectangle_;
	Point<DIM, real_number> LowerLeft_;
	Point<DIM, real_number> LowerRight_;
	Point<DIM, real_number> UpperRight_;

public:
	TriangleTestObstacle(Point<DIM, real_number> centre,
						 const Parameters &p,
						 real_number BaseLength,
						 real_number HeigthLength,
						 Point<DIM, real_number> vel = {0.0, 0.0},
						 real_number omega = 0.0,
						 real_number rf = 1.0);
	~TriangleTestObstacle() override = default;
	bool isInside(Point<DIM, real_number> P) override;
	bool isInsidePlusTol(Point<DIM, real_number> P);
	void AddObstacle(particles &vd);
};

#endif // OBSTACLE_HPP
