#include "Vector/vector_dist.hpp"
#include "Draw/DrawParticles.hpp"
#include <math.h>
#include <sys/stat.h>
#include <cmath>
#include <tinyxml2.h>

// Type of particles
#define BOUNDARY 0
#define FLUID 1

// Type of boundary condition
#define NO_SLIP 0
#define NEW_NO_SLIP 1

// TYPE OF SCENARIO
#define POISEUILLE 0
#define COUETTE 1
#define HYDROSTATIC 2
#define CYLINDER_ARRAY 3
#define CYLINDER_LATTICE 4
#define SQUARE 5
#define TRIANGLE 6
#define TRIANGLE_EQUILATERAL 7
#define CAVITY 8
#define STEP 9
#define TAYLOR_COUETTE 10
#define MOVING_OBSTACLE 11

// TYPE OF KERNERL
// #define CUBIC 0
// #define QUINTIC 1

// TYPE OF DENSITY COMPUTATION
#define DENSITY_SUMMATION 0
#define DENSITY_DIFFERENTIAL 1

// DIMENSIONALITY, changes the normalizations of the kernels
#define DIM 2

// Initial particle spacing
double dp;

class Parameters
{
public:
	int SCENARIO;
	int BC_TYPE;
	int DENSITY_TYPE;
	int WRITER;				// VTK_WRITER or CSV_WRITER
	int PROBES_ENABLED = 0; // 0 for disabled, 1 for enabled

	//////// DECLARATION OF GLOBAL PARAMETERS /////////////////////////////////////////////////////////////
	// Output file name
	std::string filename;
	std::vector<std::string> probe_filenames;
	// Physical size of the fluid domain, it goes from (0,0,0) to (length[0],length[1],length[2])
	// First fluid particle will be placed at (dp/2,dp/2,dp/2) and the last particle will be placed at (length[0]-dp/2,length[1]-dp/2,length[2]-dp/2)
	double length[DIM];

	// Boundary conditions
	size_t bc[DIM];
	// Number of boundary particles in each direction
	size_t Nboundary[DIM];
	// Number of fluid particles in each direction
	size_t Nfluid[DIM];

	// problem specific length scale to compute the Reynolds number
	double LengthScale;
	double rf;
	// Factor relating H (smoothing length) and dp (particle spacing)
	const double Hconst = 1.0;
	// Smoothing length
	double H;
	// Radius of the kernel support
	double r_threshold;
	// Normalization constant for the kernels
	// double Kcubic;
	double Kquintic;
	// Reynolds number
	double Re;
	// maximum velocity
	double umax;
	// Reference density
	double rho_zero;
	// Gamma in eq of state
	const double gamma_ = 7.0;
	// Constant used for the sound speed, number of times the max velocity
	const double coeff_sound = 10.0;
	// Sound speed
	double cbar;
	// Eq of state constant ( p0 )
	double B;
	// background pressure in eq of state
	const double xi = 0.0;
	// Gravity vector and magnitude
	Point<DIM, double> gravity_vector = {0.0};
	double gravity = 0.0;
	// wall velocity
	Point<DIM, double> vw_top = {0.0};
	Point<DIM, double> vw_bottom = {0.0};
	// Mass of the fluid and boundary particles
	double MassFluid;
	double MassBound;
	// Kinematic viscosity
	double nu;
	// Dynamic viscosity
	double eta;
	// Factor relating Pbackground and B, Pbackground = Bfactor*B
	double Bfactor;
	// Background pressure in the transport force
	double Pbackground;
	// End simulation time
	double t_end;
	// Constant used to define time integration
	const double CFLnumber = 0.1;
	// Controls otput file frequency, 1 means 1 file per time unit, 10 means 10 files per time unit, etc.
	double write_const;
	// iteration counter
	size_t cnt = 0;

	Point<DIM, double> ObstacleCenter;
	Point<DIM, double> ObstacleVelocity;
	double ObstacleOmega;

	// Taylor Couette
	double Rin;
	double Rout;
	double Win;
	double Wout;

	// custom string
	std::string custom_string = "";
};

// int SCENARIO;
// int BC_TYPE;
// //  int KERNEL = QUINTIC;
// int DENSITY_TYPE;
// int WRITER; // VTK_WRITER or CSV_WRITER
// // const int INTERACTION_LIMITER = NO_LIMITER;
// int PROBES_ENABLED = 0; // 0 for disabled, 1 for enabled

// //////// DECLARATION OF GLOBAL PARAMETERS /////////////////////////////////////////////////////////////
// // Output file name
// std::string filename;
// std::vector<std::string> probe_filenames;
// // Initial particle spacing
// double dp;
// // Physical size of the fluid domain, it goes from (0,0,0) to (length[0],length[1],length[2])
// // First fluid particle will be placed at (dp/2,dp/2,dp/2) and the last particle will be placed at (length[0]-dp/2,length[1]-dp/2,length[2]-dp/2)
// double length[DIM];
// // problem specific length scale to compute the Reynolds number
// double LengthScale;
// // Factor relating H (smoothing length) and dp (particle spacing)
// const double Hconst = 1.0;
// // Smoothing length
// double H;
// // Radius of the kernel support
// double r_threshold;
// // Normalization constant for the kernels
// // double Kcubic;
// double Kquintic;
// // Reynolds number
// double Re;
// // maximum velocity
// double umax;
// // Reference density
// double rho_zero;
// // Gamma in eq of state
// const double gamma_ = 1.0;
// // Constant used for the sound speed, number of times the max velocity
// double coeff_sound = 10.0;
// // Sound speed
// double cbar;
// // Eq of state constant ( p0 )
// double B;
// // background pressure in eq of state
// const double xi = 0.0;
// // Gravity vector and magnitude
// Point<DIM, double> gravity_vector = {0.0};
// double gravity = 0.0;
// // wall velocity
// Point<DIM, double> vw_top = {0.0};
// Point<DIM, double> vw_bottom = {0.0};
// // Mass of the fluid and boundary particles
// double MassFluid;
// double MassBound;
// // Kinematic viscosity
// double nu;
// // Dynamic viscosity
// double eta;
// // Factor relating Pbackground and B, Pbackground = Bfactor*B
// double Bfactor;
// // Background pressure in the transport force
// double Pbackground;
// // End simulation time
// double t_end;
// // Constant used to define time integration
// const double CFLnumber = 0.1;
// // Controls otput file frequency, 1 means 1 file per time unit, 10 means 10 files per time unit, etc.
// double write_const;
// // iteration counter
// size_t cnt = 0;

Parameters GlobalParameters;

//////// ALIAS FOR THE PARTICLE PROPERTIES //////////////////////////////////////////
// FLUID or BOUNDARY
const size_t type = 0;
// Density
const int rho = 1;
// pressure
const int pressure = 2;
// Delta rho calculated in the force calculation
const int drho = 3;
// calculated force
const int force = 4;
// velocity
const int velocity = 5;
// background pressure force
const int force_transport = 6;
// transport velocity
const int v_transport = 7;
// normal vector
const int normal_vector = 8;
// curvature
const int curvature_boundary = 9;
// arc length wall
const int arc_length = 10;
// volume
const int vd_volume = 11;
// angular velocity
const int vd_omega = 12;

typedef vector_dist<DIM, double, aggregate<size_t, double, double, double, double[DIM], double[DIM], double[DIM], double[DIM], double[DIM], double, double, double[3], double>> particles;
//                                          |         |     |        |        |            |           |		     |           |	        |       |          | 		   |
//                                          |         |     |        |        |            |           |		     |	     	 |	        |       |	       |	       |
//                                       type        rho  pressure delta   force          vel         F_t           vel_t    normal curvature   arc_length     vol       omega
//                                                                 density

typedef vector_dist<DIM, double, aggregate<double>> probe_particles;

struct ModelCustom
{
	template <typename Decomposition, typename vector>
	inline void addComputation(Decomposition &dec, vector &vd, size_t v, size_t p)
	{
		if (vd.template getProp<type>(p) == FLUID)
			dec.addComputationCost(v, 4);
		else
			dec.addComputationCost(v, 4);
	}

	template <typename Decomposition>
	inline void applyModel(Decomposition &dec, size_t v)
	{
		dec.setSubSubDomainComputationCost(v, dec.getSubSubDomainComputationCost(v) * dec.getSubSubDomainComputationCost(v));
	}

	double distributionTol()
	{
		return 1.01;
	}
};

// Eq state, compute pressure given density, for all fluid particles in vd
void EqState(particles &vd)
{
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto a = it.get();

		if (vd.getProp<type>(a) == FLUID)
		{
			double rho_a = vd.template getProp<rho>(a);
			double rho_frac = rho_a / GlobalParameters.rho_zero;

			vd.template getProp<pressure>(a) = GlobalParameters.B * (std::pow(rho_frac, GlobalParameters.gamma_) - 1.0) + GlobalParameters.xi;
		}

		++it;
	}
}

// Inverted equation of state, compute density given pressure, particle wise
double InvEqState_particle(const double p)
{
	return GlobalParameters.rho_zero * std::pow(((p - GlobalParameters.xi) / GlobalParameters.B + 1.0), 1.0 / GlobalParameters.gamma_);
}

double EqState_particle(const double rho)
{
	return GlobalParameters.B * (std::pow(rho / GlobalParameters.rho_zero, GlobalParameters.gamma_) - 1.0) + GlobalParameters.xi;
}

// Vector utilities
double dotProduct(const Point<DIM, double> &v, const Point<DIM, double> &w)
{
	if (DIM == 2)
		return v.get(0) * w.get(0) + v.get(1) * w.get(1);
	else if (DIM == 3)
		return v.get(0) * w.get(0) + v.get(1) * w.get(1) + v.get(2) * w.get(2);
}

std::array<Point<DIM, double>, DIM> dyadicProduct(const Point<DIM, double> &v, const Point<DIM, double> &w)
{

	std::array<Point<DIM, double>, DIM> dyad;
	dyad[0] = v.get(0) * w;
	dyad[1] = v.get(1) * w;
	if (DIM == 3)
		dyad[2] = v.get(2) * w;
	return dyad;
}

Point<DIM, double> matVec(const std::array<Point<DIM, double>, DIM> &m, const Point<DIM, double> &v)
{
	Point<DIM, double> res;
	res.get(0) = dotProduct(m[0], v);
	res.get(1) = dotProduct(m[1], v);
	if (DIM == 3)
		res.get(2) = dotProduct(m[2], v);
	return res;
}

double getVectorNorm(const Point<DIM, double> &v)
{
	return sqrt(norm2(v));
}

void normalizeVector(Point<DIM, double> &v)
{
	if (DIM == 2 && (v.get(0) == 0.0 && v.get(1) == 0.0))
	{
		v.get(0) = 0.0;
		v.get(1) = 0.0;
		return;
	}
	else if (DIM == 3 && (v.get(0) == 0.0 && v.get(1) == 0.0 && v.get(2) == 0.0))
	{
		v.get(0) = 0.0;
		v.get(1) = 0.0;
		v.get(2) = 0.0;
		return;
	}

	const double norm = getVectorNorm(v);
	v = v / norm;
}
Point<DIM, double> getPerpendicularUnit(const Point<DIM, double> &v)
{
	Point<DIM, double> perp;
	perp.get(0) = -v.get(1);
	perp.get(1) = v.get(0);
	if (DIM == 3)
		perp.get(2) = v.get(2);
	normalizeVector(perp);
	return perp;
}

// Kernel functions
// double Cubic_W(double r, const double sl)
// {
// 	r /= sl;
// 	if (r >= 0.0 && r < 1.0)
// 		return Kcubic * (1.0 - 1.5 * r * r + 0.75 * r * r * r);
// 	else if (r >= 1.0 && r < 2.0)
// 		return Kcubic * (0.25 * (2.0 - r) * (2.0 - r) * (2.0 - r));
// 	else
// 		return 0.0;
// }

// Point<DIM, double> Grad_Cubic_W(const Point<DIM, double> &dx, const double r, const double sl)
// {
// 	Point<DIM, double> DW;
// 	const double q = r / sl;

// 	const double c1 = Kcubic * (-3.0 / H);
// 	const double d1 = Kcubic * (9.0 / 4.0 / H);
// 	const double c2 = Kcubic * (-3.0 / 4.0 / H);

// 	double factor;
// 	if (q > 0.0 && q < 1.0)
// 		factor = c1 * q + d1 * q * q / r;
// 	else if (q >= 1.0 && q < 2.0)
// 		factor = c2 * (2.0 - q) * (2.0 - q) / r;
// 	else
// 		factor = 0.0;

// 	DW = factor * dx;

// 	return DW;
// }
// double Wab(double r, const double sl)
// {
// 	if (KERNEL == CUBIC)
// 		return Cubic_W(r, sl);
// 	else if (KERNEL == QUINTIC)
// 		return Quintic_W(r, sl);
// }
// Point<DIM, double> DWab(const Point<DIM, double> &dx, const double r)
// {
// 	if (KERNEL == CUBIC)
// 		return Grad_Cubic_W(dx, r);
// 	else if (KERNEL == QUINTIC)
// 		return Grad_Quintic_W(dx, r);
// }
double Wab(double r)
{
	r /= GlobalParameters.H;
	const double tmp3 = (3.0 - r) * (3.0 - r) * (3.0 - r) * (3.0 - r) * (3.0 - r);
	const double tmp2 = (2.0 - r) * (2.0 - r) * (2.0 - r) * (2.0 - r) * (2.0 - r);
	const double tmp1 = (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
	if (r >= 0.0 && r < 1.0)
		return GlobalParameters.Kquintic * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1);
	else if (r >= 1.0 && r < 2.0)
		return GlobalParameters.Kquintic * (tmp3 - 6.0 * tmp2);
	else if (r >= 2.0 && r < 3.0)
		return GlobalParameters.Kquintic * tmp3;
	else
		return 0.0;
}

Point<DIM, double> DWab(const Point<DIM, double> &dx, const double r)
{
	Point<DIM, double> DW;
	const double q = r / GlobalParameters.H;

	const double tmp3 = (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q);
	const double tmp2 = (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q);
	const double tmp1 = (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q);

	const double c1 = GlobalParameters.Kquintic * (-5.0 / GlobalParameters.H);
	double factor;

	if (q > 0.0 && q < 1.0)
		factor = c1 * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1) / r;
	else if (q >= 1.0 && q < 2.0)
		factor = c1 * (tmp3 - 6.0 * tmp2) / r;
	else if (q >= 2.0 && q < 3.0)
		factor = c1 * tmp3 / r;
	else
		factor = 0.0;

	DW = factor * dx;

	return DW;
}

Point<DIM, double> Pi_physical(const Point<DIM, double> &dr, const double &r, const Point<DIM, double> &dv, const Point<DIM, double> &dW)
{
	return GlobalParameters.eta * (dv * dotProduct(dr, dW)) / (r * r);
}

double PressureForce(const double &rhoa, const double &rhob, const double &prsa, const double &prsb)
{
	return -1.0 * (rhob * prsa + rhoa * prsb) / (rhoa + rhob);
}

std::array<Point<DIM, double>, 3> GetBoundaryPositions(const Point<DIM, double> &r, const Point<DIM, double> &normal)
{
	// aplies offset to a vector and returns the vectors pointing at the virtual wall particles
	// r needs to be pointing to the wall particle

	Point<DIM, double> offset_1 = -0.5 * normal * dp; // offset from wall to first particle
	Point<DIM, double> offset_2 = -1.0 * normal * dp; // offset from first to second particle, and from second to third
	std::array<Point<DIM, double>, 3> r_virtual;
	r_virtual[0] = r + offset_1;
	r_virtual[1] = r_virtual[0] + offset_2;
	r_virtual[2] = r_virtual[1] + offset_2;

	return r_virtual;
}

template <typename CellList>
void calcFluidVec(particles &vd, CellList &NN)
{
	// This function computes the vector pointing to the average position of the fluid particles with respect to the boundary particle

	auto part = vd.getDomainIterator();

	// Update the cell-list
	vd.updateCellList(NN);

	// For each particle ...
	while (part.isNext())
	{
		// Key of the particle a
		vect_dist_key_dx a = part.get();

		// if particle BOUNDARY
		if (vd.getProp<type>(a) == BOUNDARY)
		{
			// Get the position xa of the particle a
			Point<DIM, double> xa = vd.getPos(a);

			// initialize sum
			Point<DIM, double> r_fluid_sum = {0.0, 0.0};

			// neighborhood particles
			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

			// iterate the neighborhood particles
			while (Np.isNext() == true)
			{
				// Key of b particle
				const unsigned long b = Np.get();

				if (vd.getProp<type>(b) == FLUID)
				{
					// Get the position xb of the particle b
					const Point<DIM, double> xb = vd.getPos(b);

					// Get the vector pointing at fluid from boundary
					Point<DIM, double> dr = xb - xa;

					// take the norm squared of this vector
					const double r2 = norm2(dr);

					// If the particles interact ...
					if (r2 < GlobalParameters.r_threshold * GlobalParameters.r_threshold)
					{
						normalizeVector(dr);
						r_fluid_sum += dr;
					}
				}

				++Np;
			}
			normalizeVector(r_fluid_sum);

			for (int xyz = 0; xyz < DIM; ++xyz)
			{
				vd.template getProp<normal_vector>(a)[xyz] = r_fluid_sum.get(xyz);
			}
			// store in normal temporally, this vector is only needed for computing normals oriented "outside"
		}

		++part;
	}
}

template <typename CellList>
void calcNormalVec(particles &vd, CellList &NN)
{
	// This function computes the normal vector for a boundary particle based on the other boundary particles inside its kernel.
	// it computes the average of the perpendicular vectors to the vectors pointing to the other boundary particles ( oriented to the fluid )

	auto part = vd.getDomainIterator();

	// Update the cell-list
	vd.updateCellList(NN);

	// For each particle ...
	while (part.isNext())
	{
		// Key of the particle a
		vect_dist_key_dx a = part.get();

		// if particle BOUNDARY
		if (vd.getProp<type>(a) == BOUNDARY)
		{
			// Get the position xa of the particle a
			Point<DIM, double> xa = vd.getPos(a);

			// get vector that points at fluid
			Point<DIM, double> Rfluid = vd.getProp<normal_vector>(a);

			// initialize sum
			Point<DIM, double> n_sum = {0.0};

			// neighborhood particles
			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

			// iterate the neighborhood particles
			while (Np.isNext() == true)
			{
				// Key of b particle
				const unsigned long b = Np.get();
				if (a.getKey() == b)
				{
					++Np;
					continue;
				}

				if (vd.getProp<type>(b) == BOUNDARY)
				{
					// Get the position xb of the particle b
					const Point<DIM, double> xb = vd.getPos(b);

					// Get the vector pointing at a from b
					const Point<DIM, double> dr = xa - xb;

					// take the norm squared of this vector
					const double r2 = norm2(dr);

					// If the particles interact ...
					if (r2 < GlobalParameters.r_threshold * GlobalParameters.r_threshold)
					{
						// get perpendicular vector to dr
						Point<DIM, double> perp = getPerpendicularUnit(dr); // this is normalized

						// get scalar product of perp and Rfluid
						const double perp_dot_fluid = dotProduct(perp, Rfluid);

						// we want perp to point towards the fluid
						if (perp_dot_fluid < 0.0)
							perp = -1.0 * perp;

						// evaluate kernel
						double W = Wab(sqrt(r2));
						n_sum += perp * W;
					}
				}

				++Np;
			}
			// normalize the summed vector
			normalizeVector(n_sum);

			// store in normal vector
			for (int xyz = 0; xyz < DIM; ++xyz)
			{
				vd.template getProp<normal_vector>(a)[xyz] = n_sum.get(xyz);
			}
		}

		++part;
	}
}

inline unsigned int clockDifference(const int a, const int b, const int maxval)
{
	const unsigned int diff = std::abs(a - b);
	return std::min(diff, maxval - diff);
}

template <typename CellList>
void calcCurvature(particles &vd, CellList &NN, Vcluster<> &v_cl)
{
	// This function computes the curvature of the boundary particles from the divergence of the normal vector

	// // previous step, find the maxcounter in the particles. (stored  in pressure of marker particles)
	// auto part1 = vd.getDomainIterator();

	// double max_counter = -1e10;
	// // For each particle ...
	// while (part1.isNext())
	// {
	// 	// Key of the particle a
	// 	vect_dist_key_dx a = part1.get();

	// 	// if particle BOUNDARY
	// 	if (vd.getProp<type>(a) == BOUNDARY)
	// 	{
	// 		if (vd.getProp<pressure>(a) > max_counter)
	// 		{
	// 			max_counter = vd.getProp<pressure>(a);
	// 		}
	// 	}
	// 	++part1;
	// }
	// unsigned int max_counter_int = (int)std::round(max_counter);

	// v_cl.max(max_counter_int);
	// v_cl.execute();

	// std::cout << "max counter: " << max_counter_int << std::endl;
	// Now the real calculation

	auto part = vd.getDomainIterator();

	// Update the cell-list
	vd.updateCellList(NN);

	// max curvature is determined form the derivation of the volume formula, for curvature higher than this the volume of the third particle
	// is no longer between the two circles. Actually for 1/(2.5*dp) it becomes 0 and later negative
	const double max_curvature = 1.0 / (3.0 * dp);

	// For each particle ...
	while (part.isNext())
	{
		// Key of the particle a
		vect_dist_key_dx a = part.get();

		// if particle BOUNDARY
		if (vd.getProp<type>(a) == BOUNDARY)
		{
			// Get the position xa of the particle a
			Point<DIM, double> xa = vd.getPos(a);

			// get normal of a
			Point<DIM, double> normal_a = vd.getProp<normal_vector>(a);

			// initialize sums
			double K_sum = 0.0;
			double w_sum = 0.0;

			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

			// iterate the neighborhood particles
			while (Np.isNext() == true)
			{
				// Key of b particle
				const unsigned long b = Np.get();

				if (vd.getProp<type>(b) == BOUNDARY)
				{
					// Get the position xb of the particle b
					const Point<DIM, double> xb = vd.getPos(b);

					// Get the vector pointing at a from b
					const Point<DIM, double> dr = xa - xb;

					// take the norm squared of this vector
					const double r2 = norm2(dr);

					// If the particles interact ...
					if (r2 < GlobalParameters.r_threshold * GlobalParameters.r_threshold)
					{
						// OLD CALCULATION
						// get normal vector of b
						// Point<DIM, double> normal_b = vd.getProp<normal_vector>(b);

						// // evaluate kernel
						// double W = Wab(sqrt(r2));
						// // evaluate kernel gradient
						// Point<DIM, double> dW = DWab(dr, sqrt(r2));

						// K_sum += dotProduct(normal_b - normal_a, dW); // divergence of the normal vector
						// w_sum += W;

						// // NEW CALCULATION
						if (a.getKey() != b)
						{
							Point<DIM, double> normal_b = vd.getProp<normal_vector>(b);
							double r = sqrt(r2);
							double W = Wab(r);
							Point<DIM, double> eab = -1.0 * dr / r;
							double local_k = dotProduct(normal_b - normal_a, eab) / r;

							if (local_k > max_curvature)
							{
								local_k = max_curvature;
							}

							K_sum += local_k * W;
							w_sum += W;
						}
					}
				}

				++Np;
			}

			K_sum = K_sum / w_sum;
			// store in curvature
			vd.template getProp<curvature_boundary>(a) = K_sum;
		}
		++part;
	}
}

void calcVolume(particles &vd)
{
	// This function computes the curvature of the boundary particles from the divergence of the normal vector

	auto part = vd.getDomainIterator();

	// For each particle ...
	while (part.isNext())
	{
		// Key of the particle a
		vect_dist_key_dx a = part.get();

		// if particle BOUNDARY
		if (vd.getProp<type>(a) == BOUNDARY)
		{
			double kappa = vd.getProp<curvature_boundary>(a);
			double dxwall = vd.getProp<arc_length>(a);

			for (int i = 0; i < 3; i++)
			{
				vd.template getProp<vd_volume>(a)[i] = 0.5 * (2.0 * dp + dp * dp * kappa - 2.0 * (i + 1.0) * dp * dp * kappa) * dxwall;
			}
		}
		++part;
	}
}

template <typename CellList>
void calc_density(particles &vd, CellList &NN)
{
	// This function computes the density of particles from the summation of the kernel

	auto part = vd.getDomainIterator();

	// Update the cell-list
	vd.updateCellList(NN);

	// For each particle ...
	while (part.isNext())
	{
		// Key of the particle a
		vect_dist_key_dx a = part.get();

		// if particle FLUID
		if (vd.getProp<type>(a) == FLUID)
		{
			double massa = GlobalParameters.MassFluid;

			// Get the position xb of the particle a
			Point<DIM, double> xa = vd.getPos(a);

			// initialize density sum
			double rho_sum = 0.0;
			// double w_sum = 0.0;
			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

			// iterate the neighborhood particles
			while (Np.isNext() == true)
			{
				// Key of b particle
				unsigned long b = Np.get();

				// // if (a == b) skip this particle
				// if (a.getKey() == b)
				// {
				// 	++Np;
				// 	continue;
				// };

				// Get the position xb of the particle b
				const Point<DIM, double> xb = vd.getPos(b);

				// Get the vector pointing at xa from xb
				const Point<DIM, double> dr = xa - xb;

				// take the norm squared of this vector
				const double r2 = norm2(dr);
				const double r = sqrt(r2);

				// If the particles interact ...
				if (r < GlobalParameters.r_threshold)
				{

					if (vd.getProp<type>(b) == FLUID)
					{

						// evaluate kernel
						const double w = Wab(r);
						// w_sum += w * dp * dp;

						rho_sum += w * GlobalParameters.MassFluid;
					}
					else
					{
						if (GlobalParameters.BC_TYPE == NO_SLIP)
						{
							const double r = sqrt(r2);

							// evaluate kernel
							const double w = Wab(r);
							// w_sum += w * dp * dp;
							rho_sum += w * GlobalParameters.MassBound;
						}
						else if (GlobalParameters.BC_TYPE == NEW_NO_SLIP) // need to evaluate kernel at dummy particles
						{
							// get normal vector of b
							const Point<DIM, double> normal = vd.getProp<normal_vector>(b);
							// get volumes, and curvature of dummy particles
							const Point<3, double> vol = vd.template getProp<vd_volume>(b);
							// const double kappa = vd.template getProp<curvature_boundary>(b);
							// Apply offsets to dr to get 3 vectrors pointing to dummy particles
							const std::array<Point<DIM, double>, 3> R_dummy = GetBoundaryPositions(-1.0 * dr, normal);

							// distance to the marker particle
							// const double dist2marker = sqrt(r2);

							// double k = vd.template getProp<vd_volume>(b)[0];
							for (int i = 0; i < 3; i++)
							{
								// double rmax = sqrt(3.0 * 3.0 - (0.5 + (double)i) * (0.5 + (double)i)) * dp;
								// double rmin = (3.0 - (0.5 + (double)i)) * dp;
								// double kappa_max = 1.0 / (3.0 * dp);

								// kappa 0 gets rmax, kappa = kappa_max gets rmin
								// double r_interp = (rmin - rmax) / kappa_max * kappa + rmax;
								////
								// if (dist2marker < r_interp)
								// {

								const double W = Wab(getVectorNorm(R_dummy[i]));
								rho_sum += W * vol[i] * GlobalParameters.rho_zero; // W*mass
																				   // }
							}
						}
					}
				}

				++Np;
			}
			if (rho_sum != 0.0)
			{
				vd.template getProp<rho>(a) = rho_sum;
				// std::cout << "rho: " << rho_sum << std::endl;
				// std::cout << "w_sum: " << w_sum << std::endl;
			}
			else
			{
				std::cout << "WARNING: NO particles around, density summation zero" << std::endl;
				vd.template getProp<rho>(a) = GlobalParameters.rho_zero;
			}
		}

		++part;
	}
}

template <typename CellList>
void ExtrapolateVelocity(particles &vd, CellList &NN)
{
	// This function fills the value of v_transport for the boundary particles, with the velocity for no slip BC

	auto part = vd.getDomainIterator();

	// Update the cell-list
	vd.updateCellList(NN);

	// For each particle ...
	while (part.isNext())
	{
		// Key of the b particle
		vect_dist_key_dx b = part.get();

		// if particle boundary
		if (vd.getProp<type>(b) != FLUID)
		{

			// Get the position xb of the boundary particle
			Point<DIM, double> xb = vd.getPos(b);

			// Get the Velocity of the boundary particle
			Point<DIM, double> vb = vd.getProp<velocity>(b);

			// initialize sums
			Point<DIM, double> sum_vW = (DIM == 2) ? Point<DIM, double>{0.0, 0.0} : Point<DIM, double>{0.0, 0.0, 0.0};

			double sum_pW = 0.0;
			double sum_W = 0.0;

			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(b)));
			// iterate the neighborhood fluid particles
			while (Np.isNext() == true)
			{
				// Key of fluid particle
				unsigned long f = Np.get();

				// if (b == f) skip this particle
				if (b.getKey() == f)
				{
					++Np;
					continue;
				};
				// Skip other boundary particles
				if (vd.getProp<type>(f) != FLUID)
				{
					++Np;
					continue;
				}

				// Get the position xf of the fluid particle
				Point<DIM, double> xf = vd.getPos(f);

				// Get the velocity of the fluid particle
				Point<DIM, double> vf = vd.getProp<velocity>(f);

				// Get the density of the fluid particle
				double rhof = vd.getProp<rho>(f);

				// Get the pressure of the fluid particle
				double Pf = vd.getProp<pressure>(f);

				// Get the vector pointing at xb from xf rwf
				Point<DIM, double> dr = xb - xf;

				// take the norm squared of this vector
				double r2 = norm2(dr);

				// If the particles interact ...
				if (r2 < GlobalParameters.r_threshold * GlobalParameters.r_threshold)
				{
					// calculate distance
					double r = sqrt(r2);

					// evaluate kernel
					double w = Wab(r);
					// compute v*W
					sum_vW += w * vf;

					// compute Pf +rhof*(g-a) dot rwf
					// at the moment a = 0
					const double dot = dotProduct(dr, GlobalParameters.gravity_vector);
					sum_pW += w * (Pf + rhof * dot);
					sum_W += w;
				}

				++Np;
			}
			if (sum_W != 0.0)
			{
				// Set the velocity of the boundary particle b ( to be used in the momentum equation to impose BC)
				// We use the v_transport field because boundary particles dont use this array

				// marker particles store centre of solid body in force_transport since it is unused
				// vector pointing from centre of rotation to marker particle
				const Point<DIM, double> radial_vec = {xb.get(0) - vd.getProp<force_transport>(b)[0], xb.get(1) - vd.getProp<force_transport>(b)[1]};
				const double radial = getVectorNorm(radial_vec);
				// vector tangential to the radial vector, rotation velocity is in this direction
				const Point<DIM, double> rotation_tangential = getPerpendicularUnit(radial_vec);
				Point<DIM, double> vw = vd.template getProp<velocity>(b);
				vw.get(0) += radial * vd.getProp<vd_omega>(b) * rotation_tangential.get(0);
				vw.get(1) += radial * vd.getProp<vd_omega>(b) * rotation_tangential.get(1);

				for (int xyz = 0; xyz < DIM; ++xyz)
				{
					vd.template getProp<v_transport>(b)[xyz] = 2.0 * vw[xyz] - sum_vW.get(xyz) / sum_W;
				}

				// Set the pressure of the boundary particle b
				vd.template getProp<pressure>(b) = sum_pW / sum_W;
				// Compute density from inverted Eq of state
				vd.template getProp<rho>(b) = InvEqState_particle(vd.template getProp<pressure>(b));
			}
			else
			{

				for (int xyz = 0; xyz < DIM; ++xyz)
				{
					vd.template getProp<v_transport>(b)[xyz] = 2.0 * vd.template getProp<velocity>(b)[xyz];
				}

				vd.template getProp<pressure>(b) = 0.0;
				vd.template getProp<rho>(b) = GlobalParameters.rho_zero;
			}
		}

		++part;
	}
}

// function to know the type of a variable
// used to know types of auto variables
template <class T>
std::string type_name()
{
	typedef typename std::remove_reference<T>::type TR;
	std::unique_ptr<char, void (*)(void *)> own(
#ifndef _MSC_VER
		abi::__cxa_demangle(typeid(TR).name(), nullptr,
							nullptr, nullptr),
#else
		nullptr,
#endif
		std::free);
	std::string r = own != nullptr ? own.get() : typeid(TR).name();
	if (std::is_const<TR>::value)
		r += " const";
	if (std::is_volatile<TR>::value)
		r += " volatile";
	if (std::is_lvalue_reference<T>::value)
		r += "&";
	else if (std::is_rvalue_reference<T>::value)
		r += "&&";
	return r;
}

void interact_fluid_boundary_new(particles &vd,
								 vect_dist_key_dx fluid_key,
								 const double &massf,
								 const double &rhof,
								 const double &Pf,
								 const Point<DIM, double> &xf,
								 const Point<DIM, double> &vf,
								 const Point<DIM, double> &xw,
								 const Point<DIM, double> &r_wall_to_fluid,
								 unsigned long &boundary_key,
								 bool accumulate_force,
								 double &cylinder_force)
{

	// Points from fluid to wall
	Point<DIM, double> r_fluid_to_wall = -1.0 * r_wall_to_fluid;
	double dist2marker = getVectorNorm(r_fluid_to_wall);
	Point<DIM, double> vw = vd.getProp<velocity>(boundary_key);

	double ang_vel = vd.getProp<vd_omega>(boundary_key);

	if (ang_vel != 0.0) // if the solid is rotating, we need to add the tangential velocity of the rotation to vw
	{
		// marker particles store centre of solid body in force_transport since it is unused
		// vector pointing from centre of rotation to marker particle
		const Point<DIM, double> radial_vec = {xw.get(0) - vd.getProp<force_transport>(boundary_key)[0],
											   xw.get(1) - vd.getProp<force_transport>(boundary_key)[1]};
		const double radius = getVectorNorm(radial_vec);

		// get vector tangential to the radial vector, rotation velocity is in this direction
		const Point<DIM, double> tangential_rotation = getPerpendicularUnit(radial_vec);

		// Wall velocity is linear velocity + w*R*tangential
		vw.get(0) += radius * ang_vel * tangential_rotation.get(0);
		vw.get(1) += radius * ang_vel * tangential_rotation.get(1);
	}

	// Get normal and tangential vectors for velocity mirroring
	const Point<DIM, double> normal = vd.getProp<normal_vector>(boundary_key);
	const Point<DIM, double> tangential = getPerpendicularUnit(normal);

	// wall acceleration
	const Point<DIM, double> aw = vd.getProp<force>(boundary_key);

	// Difference between fluid transport and momentum velocity
	const Point<DIM, double> vtf = vd.getProp<v_transport>(fluid_key);
	const Point<DIM, double> vdiff_f = vtf - vf;

	// get curvature, and arc length
	double kappa = vd.getProp<curvature_boundary>(boundary_key);

	// Project vf and vw on tangential and normal directions
	double vt = dotProduct(vf, tangential);
	double vn = dotProduct(vf, normal);
	double vwt = dotProduct(vw, tangential);

	// vertical distance from fluid particle to wall
	double lf = dotProduct(r_fluid_to_wall, normal);
	lf = (lf < 0.0 ? -1.0 * lf : lf); // absolute value

	// Get array of vectors from fluid to 3 boundary particles and its norms
	std::array<Point<DIM, double>, 3> r_boundary = GetBoundaryPositions(r_fluid_to_wall, normal);
	std::array<double, 3> r_boundary_norm = {getVectorNorm(r_boundary[0]), getVectorNorm(r_boundary[1]), getVectorNorm(r_boundary[2])};

	// const double dist2third = r_boundary_norm[2];
	// distance from 3 boundary particles to marker
	std::array<double, 3> lwall = {0.5 * dp, 1.5 * dp, 2.5 * dp};

	// project forces on normal direction
	double g_normal = dotProduct(GlobalParameters.gravity_vector, normal);
	double a_normal = dotProduct(aw, normal);

	// For gradient of Af tensor
	const std::array<Point<DIM, double>, DIM> Af = dyadicProduct(rhof * vf, vdiff_f);
	// to avoid division by zero
	lf = std::max(lf, 0.25 * dp);

	const Point<3, double> vol = vd.getProp<vd_volume>(boundary_key);

	for (int i = 0; i < 3; i++) // for the 3 boundary particles
	{
		// double rmax = sqrt(3.0 * 3.0 - (0.5 + (double)i) * (0.5 + (double)i)) * dp;
		// double rmin = (3.0 - (0.5 + (double)i)) * dp;
		// double kappa_max = 1.0 / (3.0 * dp);

		// kappa 0 gets rmax, kappa = kappa_max gets rmin
		// double r_interp = (rmin - rmax) / kappa_max * kappa + rmax;
		// double r_interp = rmin;
		// if (dist2marker < r_interp)
		// {
		const double Mass_boundary = vol[i] * GlobalParameters.rho_zero;
		const Point<DIM, double> v_boundary = ((vwt - vt) * (lwall[i] / lf) + vwt) * tangential + vn * normal;
		const double p_boundary = Pf - rhof * (g_normal - a_normal) * (lf + lwall[i]); //  dot(r_boundary,normal) = -(lf+lw)
		const double rho_boundary = InvEqState_particle(p_boundary);

		// flip sign of r_boundary to get vector pointing from boundary to fluid (Force routines use the vector pointing from b to a)
		r_boundary[i] = -1.0 * r_boundary[i];

		// Evaluate kernel gradient
		const Point<DIM, double> DW = DWab(r_boundary[i], r_boundary_norm[i]);

		// Compute forces
		const Point<DIM, double> v_rel = vf - v_boundary;
		// const double Va2 = (massf / rhof) * (massf / rhof);
		// const double Vb2 = Volumes[i] * Volumes[i];
		// const double Vb2 = (Mass_boundary / rho_boundary) * (Mass_boundary / rho_boundary); // this is Vol_boundary^2*(rhob/rho0)^2 therefore considers density changes,
		const double Vb = (Mass_boundary / rho_boundary);
		// But I found no big difference with the commented out expression

		const Point<DIM, double> ViscosityTerm = Pi_physical(r_boundary[i], r_boundary_norm[i], v_rel, DW);
		const double PressureTerm = PressureForce(rhof, rho_boundary, Pf, p_boundary); //-p_boundary - Pf;
		const Point<DIM, double> GradATerm = 0.5 * matVec(Af, DW);

		// std::cout << "PressureTerm_x: " << PressureTerm * DW.get(0) << "PressureTerm_y: " << PressureTerm * DW.get(1) << std::endl;
		// std::cout << "ViscosityTerm_x: " << ViscosityTerm.get(0) << "ViscosityTerm_y: " << ViscosityTerm.get(1) << std::endl;
		// std::cout << "GradATerm_x: " << GradATerm.get(0) << "GradATerm_y: " << GradATerm.get(1) << std::endl;
		// if (accumulate_force) // we accumulate x force on cylinder ( this is just to compute drag coefficient)
		// {
		// 	cylinder_force += -1.0 * (Va2 + Vb2) * (PressureTerm * DW.get(0) + ViscosityTerm.get(0) + GradATerm.get(0)) / massf;
		// }
		for (int xyz = 0; xyz < DIM; ++xyz)
		{
			// write to particles
			vd.getProp<force>(fluid_key)[xyz] += 2.0 * (Vb / rhof) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + GradATerm.get(xyz));
			vd.getProp<force_transport>(fluid_key)[xyz] += -2.0 * (Vb / rhof) * (GlobalParameters.Pbackground) * DW.get(xyz);
			// vd.getProp<force>(fluid_key)[xyz] += (Va2 + Vb2) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + GradATerm.get(xyz)) / massf;
			// vd.getProp<force_transport>(fluid_key)[xyz] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(xyz) / massf;
		}

		if (GlobalParameters.DENSITY_TYPE == DENSITY_DIFFERENTIAL) // this doesnt work well I havent touched in a long time and I have made may changes
		{
			vd.getProp<drho>(fluid_key) += Mass_boundary * dotProduct(v_rel, DW);
		}
		// }
	}
}

void interact_fluid_boundary_old(particles &vd,
								 vect_dist_key_dx &fluid_key,
								 const double &massf,
								 const double &rhof,
								 const double &Pf,
								 const Point<DIM, double> &xf,
								 const Point<DIM, double> &vf,
								 const Point<DIM, double> &xb,
								 const Point<DIM, double> &r_ab,
								 const double &r2,
								 unsigned long &boundary_key,
								 bool accumulate_force,
								 double &cylinder_force)
{
	const double massb = GlobalParameters.MassBound;
	const double rhob = vd.getProp<rho>(boundary_key);
	const double Pb = vd.getProp<pressure>(boundary_key);
	const Point<DIM, double> vb = vd.getProp<velocity>(boundary_key);
	const Point<DIM, double> vb_noslip = vd.getProp<v_transport>(boundary_key); // here we store the extrapolated velocity for no slip BC

	const Point<DIM, double> r_fb = xf - xb;
	const double r = sqrt(r2);

	const Point<DIM, double> v_rel = vf - vb;
	const Point<DIM, double> v_rel_aux = vf - vb_noslip;

	const Point<DIM, double> DW = DWab(r_ab, r);

	const double Va2 = (massf / rhof) * (massf / rhof);
	const double Vb2 = (massb / rhob) * (massb / rhob);

	const Point<DIM, double> ViscosityTerm = Pi_physical(r_ab, r, v_rel_aux, DW);
	const double PressureTerm = PressureForce(rhof, rhob, Pf, Pb);

	const Point<DIM, double> vtf = vd.getProp<v_transport>(fluid_key);
	const Point<DIM, double> vdiff_f = vtf - vf;
	// Boundary particles have no transport velocity difference
	std::array<Point<DIM, double>, DIM> Af = dyadicProduct(rhof * vf, vdiff_f);
	Point<DIM, double> GradATerm = 0.5 * matVec(Af, DW);

	if (accumulate_force) // we accumulate x force on cylinder
	{
		cylinder_force += -1.0 * (Va2 + Vb2) * (PressureTerm * DW.get(0) + ViscosityTerm.get(0) + GradATerm.get(0)) / massf;
		// cylinder_force += -1.0 * (Va2 + Vb2) * Pbackground * DW.get(0) / massf;
	}
	for (int xyz = 0; xyz < DIM; ++xyz)
	{
		vd.getProp<force>(fluid_key)[xyz] += (Va2 + Vb2) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + GradATerm.get(xyz)) / massf;
		vd.getProp<force_transport>(fluid_key)[xyz] += -1.0 * (Va2 + Vb2) * GlobalParameters.Pbackground * DW.get(xyz) / massf;
	}

	if (GlobalParameters.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
	{
		vd.getProp<drho>(fluid_key) += massb * dotProduct(v_rel, DW);
	}
}

void interact_fluid_fluid(particles &vd,
						  const vect_dist_key_dx &a_key,
						  const double &massa,
						  const double &rhoa,
						  const double &Pa,
						  const Point<DIM, double> &xa,
						  const Point<DIM, double> &va,
						  const Point<DIM, double> &xb,
						  const Point<DIM, double> &r_ab,
						  const double &r2,
						  const unsigned long &b_key)
{

	const double massb = GlobalParameters.MassFluid;
	const double rhob = vd.getProp<rho>(b_key);
	const double Pb = vd.getProp<pressure>(b_key);
	const Point<DIM, double> vb = vd.getProp<velocity>(b_key);

	const double r = sqrt(r2);

	const Point<DIM, double> v_rel = va - vb;

	const Point<DIM, double> DW = DWab(r_ab, r);

	const double Va2 = (massa / rhoa) * (massa / rhoa);
	const double Vb2 = (massb / rhob) * (massb / rhob);

	const Point<DIM, double> ViscosityTerm = Pi_physical(r_ab, r, v_rel, DW);
	const double PressureTerm = PressureForce(rhoa, rhob, Pa, Pb);

	const Point<DIM, double> vta = vd.getProp<v_transport>(a_key);
	const Point<DIM, double> vtb = vd.getProp<v_transport>(b_key);
	const Point<DIM, double> vdiff_a = vta - va;
	const Point<DIM, double> vdiff_b = vtb - vb;

	// double GradATerm = 0.5 * (rhoa * dotProduct(va, vdiff_a) + rhob * dotProduct(vb, vdiff_b));
	std::array<Point<DIM, double>, DIM> Aa = dyadicProduct(rhoa * va, vdiff_a);
	std::array<Point<DIM, double>, DIM> Ab = dyadicProduct(rhob * vb, vdiff_b);
	std::array<Point<DIM, double>, DIM> SumA;
	SumA[0] = Aa[0] + Ab[0];
	SumA[1] = Aa[1] + Ab[1];
	if (DIM == 3)
		SumA[2] = Aa[2] + Ab[2];

	Point<DIM, double> GradATerm = 0.5 * matVec(SumA, DW);

	for (int xyz = 0; xyz < DIM; ++xyz)
	{
		vd.getProp<force>(a_key)[xyz] += (Va2 + Vb2) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + GradATerm.get(xyz)) / massa;
		vd.getProp<force_transport>(a_key)[xyz] += -1.0 * (Va2 + Vb2) * GlobalParameters.Pbackground * DW.get(xyz) / massa;
	}
	if (GlobalParameters.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
	{
		vd.getProp<drho>(a_key) += massb * dotProduct(v_rel, DW);
	}
}

template <typename CellList>
void calc_forces(particles &vd, CellList &NN, double &cylinder_force, bool calc_drag)
{
	vector_dist_iterator part = vd.getDomainIterator();

	const double max_curvature = 1.0 / (3.0 * dp);

	// Update the cell-list
	vd.updateCellList(NN);

	// For each particle ...
	while (part.isNext())
	{
		// get key of particle a

		auto a = part.get();

		// We threat FLUID particle differently from BOUNDARY PARTICLES
		// Boundary particles dont need any force computation

		if (vd.getProp<type>(a) == FLUID) // INTERACTION OF A FLUID PARTICLE WITH ITS NEIGHBORHOOD
		{

			// Get the position xp of the particle
			Point<DIM, double> xa = vd.getPos(a);

			// Take the mass of the particle dependently if it is FLUID or BOUNDARY
			double massa = GlobalParameters.MassFluid;

			// Get the density and pressure of the particle a
			double rhoa = vd.getProp<rho>(a);
			double Pa = vd.getProp<pressure>(a);

			// Get the Velocity of the particle a
			Point<DIM, double> va = vd.getProp<velocity>(a);

			// Reset the force counter (0 + gravity)
			for (int xyz = 0; xyz < DIM; xyz++)
			{
				vd.template getProp<force>(a)[xyz] = GlobalParameters.gravity_vector.get(xyz);
				vd.template getProp<force_transport>(a)[xyz] = 0.0;
			}
			vd.template getProp<drho>(a) = 0.0;

			// Get an iterator over the neighborhood particles of p
			auto Np = NN.getNNIterator(NN.getCell(xa));

			// For each neighborhood particle
			while (Np.isNext() == true)
			{
				// get particle b
				unsigned long b = Np.get();

				// if (a == b) skip this particle
				if (a.getKey() == b)
				{
					++Np;
					continue;
				};

				// Get the position xb of the particle
				Point<DIM, double> xb = vd.getPos(b);

				// Get the distance between a and b
				// in fluid - boundary its xf-xb i.e. vector pointing at fluid from boundary
				Point<DIM, double> dr = xa - xb;

				// take the norm (squared) of this vector
				double r2 = norm2(dr);
				bool accumulate_force = true;
				// if they interact
				if (r2 < GlobalParameters.r_threshold * GlobalParameters.r_threshold)
				{
					if (vd.getProp<type>(b) == BOUNDARY)
					{
						if (calc_drag && (xb.get(1) > 2.0 * dp || xb.get(1) < GlobalParameters.length[1] - 2.0 * dp)) // to exclude solid walls
						{
							accumulate_force = true;
						}
						else
						{
							accumulate_force = false;
						}

						if (GlobalParameters.BC_TYPE == NO_SLIP)
							interact_fluid_boundary_old(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b, accumulate_force, cylinder_force);
						else if (GlobalParameters.BC_TYPE == NEW_NO_SLIP)
						{
							const Point<DIM, double> normal_b = vd.getProp<normal_vector>(b);
							const double kappa = vd.getProp<curvature_boundary>(b);
							// if (dotProduct(normal_b, dr) > (max_cos / max_curvature) * kappa)
							// if (dotProduct(normal_b, dr) > 0.0)
							interact_fluid_boundary_new(vd, a, massa, rhoa, Pa, xa, va, xb, dr, b, accumulate_force, cylinder_force);
						}
					}
					else
					{
						interact_fluid_fluid(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b);
					}
				}
				++Np;
			}
		}

		++part;
	}
}

void max_velocity(particles &vd, Vcluster<> &v_cl, double &max_vel)
{
	// Calculate the maximum acceleration
	auto part = vd.getDomainIterator();

	while (part.isNext())
	{
		auto a = part.get();

		Point<DIM, double> vel(vd.getProp<velocity>(a));
		double vel2 = norm2(vel);

		if (vel2 >= max_vel)
			max_vel = vel2;

		++part;
	}
	max_vel = std::sqrt(max_vel);

	// Vcluster<> &v_cl = create_vcluster();
	v_cl.max(max_vel);
	v_cl.execute();
}

double calc_deltaT(particles &vd, Vcluster<> &v_cl)
{
	double Maxvel = 0.0;
	max_velocity(vd, v_cl, Maxvel);

	double dt_u = 0.25 * GlobalParameters.H / (GlobalParameters.cbar + abs(Maxvel));
	double dt_visc = 0.25 * GlobalParameters.H * GlobalParameters.H / (GlobalParameters.nu);
	double dt_g = 0.25 * sqrt(GlobalParameters.H / getVectorNorm(GlobalParameters.gravity_vector));
	double dt = GlobalParameters.CFLnumber * std::min({dt_u, dt_visc, dt_g});
	// if (dt < DtMin)
	// 	dt = DtMin;

	return dt;
}

void ApplyRotation(Point<DIM, double> &x, const double theta, const Point<DIM, double> centre)
{

	Point<DIM, double> x_rotated;
	x_rotated.get(0) = cos(theta) * (x.get(0) - centre.get(0)) - sin(theta) * (x.get(1) - centre.get(1)) + centre.get(0);
	x_rotated.get(1) = sin(theta) * (x.get(0) - centre.get(0)) + cos(theta) * (x.get(1) - centre.get(1)) + centre.get(1);

	x = x_rotated;
}

Point<DIM, double> SolidBodyAcceleration(double t)
{
	Point<DIM, double> a;
	if (GlobalParameters.SCENARIO == MOVING_OBSTACLE)
	{
		double period = 10.0;
		double amplitude = 0.25;
		a.get(0) = amplitude * (4.0 * M_PI * M_PI / (period * period)) * sin(2.0 * M_PI * t / period);
		a.get(1) = 0.0;
	}
	else
	{
		a = {0.0, 0.0};
	}

	return a;
}

template <typename CellList>
void kick_drift_int(particles &vd, CellList &NN, const double dt, Vcluster<> &v_cl, double &cylinder_force, bool calc_drag, double t)
{
	// particle iterator
	auto part = vd.getDomainIterator();

	const double dt_2 = dt * 0.5;
	vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

	Point<DIM, double> aSolid_minus = SolidBodyAcceleration(t - dt_2);
	Point<DIM, double> aSolid_plus = SolidBodyAcceleration(t + dt_2);

	// For each particle ...
	while (part.isNext())
	{
		// get particle a key
		vect_dist_key_dx a = part.get();

		if (vd.template getProp<type>(a) == BOUNDARY)
		{
			// move boundary particles with linear velocity and acceleration
			for (int xyz = 0; xyz < DIM; xyz++)
			{
				// trick to avoid accelerating the solid walls in the moving obstacle scenario
				if (vd.template getProp<vd_omega>(a) != 0.0)
					vd.template getProp<force>(a)[xyz] = aSolid_minus.get(xyz);
				else
					vd.template getProp<force>(a)[xyz] = 0.0;

				vd.template getProp<velocity>(a)[xyz] += dt_2 * vd.template getProp<force>(a)[xyz];
				vd.getPos(a)[xyz] += dt * vd.template getProp<velocity>(a)[xyz];

				// centre also moves
				vd.template getProp<force_transport>(a)[xyz] += dt * vd.template getProp<velocity>(a)[xyz];
			}

			// also rotate
			if (vd.template getProp<vd_omega>(a) != 0.0) // if rotating
			{
				double theta = vd.template getProp<vd_omega>(a) * dt;
				Point<DIM, double> normal = vd.template getProp<normal_vector>(a);

				// apply rotation
				Point<DIM, double> centre = vd.getProp<force_transport>(a); // stored in force_transport
				Point<DIM, double> x = vd.getPos(a);
				ApplyRotation(x, theta, centre);
				vd.getPos(a)[0] = x.get(0);
				vd.getPos(a)[1] = x.get(1);
				// update normals ( rotated wrt origin )
				ApplyRotation(normal, theta, {0.0, 0.0});
				vd.template getProp<normal_vector>(a)[0] = normal.get(0);
				vd.template getProp<normal_vector>(a)[1] = normal.get(1);
			}
		}
		else
		{
			for (int xyz = 0; xyz < DIM; xyz++)
			{
				vd.template getProp<velocity>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force>(a)[xyz];
				vd.template getProp<v_transport>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force_transport>(a)[xyz];
				vd.getPos(a)[xyz] += dt * vd.template getProp<v_transport>(a)[xyz];
			}

			if (GlobalParameters.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
				vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
		}
		++part;
	}
	// map particles if they have gone outside the domain
	vd.map();
	vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

	// in density summation we need to update the density after moving the particles
	if (GlobalParameters.DENSITY_TYPE == DENSITY_SUMMATION)
	{
		calc_density(vd, NN);
		vd.ghost_get<rho>(KEEP_PROPERTIES);
	}

	// Calculate pressure from the density
	EqState(vd);
	vd.ghost_get<pressure>(KEEP_PROPERTIES);

	if (GlobalParameters.BC_TYPE == NO_SLIP)
	{
		ExtrapolateVelocity(vd, NN);
		vd.ghost_get<rho, pressure, v_transport>(KEEP_PROPERTIES);
	}

	cylinder_force = 0.0;
	calc_forces(vd, NN, cylinder_force, calc_drag);
	if (calc_drag)
	{
		cylinder_force = cylinder_force * GlobalParameters.MassFluid;
		v_cl.sum(cylinder_force);
		v_cl.execute();
	}

	vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

	// particle iterator
	auto part2 = vd.getDomainIterator();

	// For each particle ...
	while (part2.isNext())
	{
		// get particle a key
		vect_dist_key_dx a = part2.get();

		if (vd.template getProp<type>(a) == BOUNDARY)
		{
			for (int xyz = 0; xyz < DIM; xyz++)
			{
				// trick to avoid accelerating the solid walls in the moving obstacle scenario
				if (vd.template getProp<vd_omega>(a) != 0.0)
					vd.template getProp<force>(a)[xyz] = aSolid_plus.get(xyz);
				else
					vd.template getProp<force>(a)[xyz] = 0.0;

				vd.template getProp<velocity>(a)[xyz] += dt_2 * vd.template getProp<force>(a)[xyz];
			}
		}
		else
		{
			for (int xyz = 0; xyz < DIM; xyz++)
			{
				vd.template getProp<velocity>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force>(a)[xyz];
				vd.template getProp<v_transport>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force_transport>(a)[xyz];
			}
		}

		++part2;
	}

	// increment the iteration counter
	GlobalParameters.cnt++;
}

void SetFilename(std::string &filename, const long int Nprc)
{
	filename = "";
	// Scenario name
	if (GlobalParameters.SCENARIO == POISEUILLE)
		filename = "Poiseuille";
	else if (GlobalParameters.SCENARIO == COUETTE)
		filename = "Couette";
	else if (GlobalParameters.SCENARIO == HYDROSTATIC)
		filename = "Hydrostatic";
	else if (GlobalParameters.SCENARIO == CYLINDER_ARRAY)
		filename = "CylinderArray";
	else if (GlobalParameters.SCENARIO == CYLINDER_LATTICE)
		filename = "CylinderLattice";
	else if (GlobalParameters.SCENARIO == SQUARE)
		filename = "Square";
	else if (GlobalParameters.SCENARIO == TRIANGLE)
		filename = "Triangle";
	else if (GlobalParameters.SCENARIO == CAVITY)
		filename = "Cavity";
	else if (GlobalParameters.SCENARIO == STEP)
		filename = "Step";
	else if (GlobalParameters.SCENARIO == TAYLOR_COUETTE)
		filename = "TaylorCouette";
	else if (GlobalParameters.SCENARIO == TRIANGLE_EQUILATERAL)
		filename = "TriangleEquilateral";
	else if (GlobalParameters.SCENARIO == MOVING_OBSTACLE)
		filename = "MovingObstacle";

	// BC name
	if (GlobalParameters.BC_TYPE == NO_SLIP)
		filename += "_OLD_BC";
	else if (GlobalParameters.BC_TYPE == NEW_NO_SLIP)
		filename += "_NEW_BC";

	// Density name
	if (GlobalParameters.DENSITY_TYPE == DENSITY_SUMMATION)
		filename += "_Summation";
	else if (GlobalParameters.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
		filename += "_Differential";

	// Add the size of the simulation, the number of processors and the refinement factor
	std::string size_proc_name = std::to_string(GlobalParameters.Nfluid[0]) + "_" + std::to_string(GlobalParameters.Nfluid[1]) + "_" + std::to_string(Nprc) + "prc";
	if (GlobalParameters.BC_TYPE == NEW_NO_SLIP)
		size_proc_name += "_" + std::to_string((int)GlobalParameters.rf) + "rf";
	filename += "_" + size_proc_name;
	filename += GlobalParameters.custom_string;
}

void WriteParameters(Vcluster<> &v_cl)
{
	const double Lx = GlobalParameters.length[0]; // channel height
	const double Ly = GlobalParameters.length[1]; // channel width
	if (DIM == 3)
		const double Lz = GlobalParameters.length[2]; // channel length ( +0.5dp due to periodicity)

	const long int Nprc = v_cl.getProcessingUnits();
	SetFilename(GlobalParameters.filename, Nprc);

	std::string constants_filename = GlobalParameters.filename + "_PARAMETERS" + ".txt";
	std::ofstream file(constants_filename);

	std::string scenario_str = "";
	if (GlobalParameters.SCENARIO == POISEUILLE)
		scenario_str = "Poiseuille";
	else if (GlobalParameters.SCENARIO == COUETTE)
		scenario_str = "Couette";
	else if (GlobalParameters.SCENARIO == HYDROSTATIC)
		scenario_str = "Hydrostatic";
	else if (GlobalParameters.SCENARIO == CYLINDER_ARRAY)
		scenario_str = "CylinderArray";
	else if (GlobalParameters.SCENARIO == CYLINDER_LATTICE)
		scenario_str = "CylinderLattice";
	else if (GlobalParameters.SCENARIO == SQUARE)
		scenario_str = "Square";
	else if (GlobalParameters.SCENARIO == TRIANGLE)
		scenario_str = "Triangle";
	else if (GlobalParameters.SCENARIO == CAVITY)
		scenario_str = "Cavity";
	else if (GlobalParameters.SCENARIO == STEP)
		scenario_str = "Step";
	else if (GlobalParameters.SCENARIO == TAYLOR_COUETTE)
		scenario_str = "TaylorCouette";
	else if (GlobalParameters.SCENARIO == TRIANGLE_EQUILATERAL)
		scenario_str = "TriangleEquilateral";
	else if (GlobalParameters.SCENARIO == MOVING_OBSTACLE)
		scenario_str = "MovingObstacle";

	std::string BC_str = "";
	if (GlobalParameters.BC_TYPE == NO_SLIP)
		BC_str = "OldBC";
	else if (GlobalParameters.BC_TYPE == NEW_NO_SLIP)
		BC_str = "NewBC";

	std::string density_str = "";
	if (GlobalParameters.DENSITY_TYPE == DENSITY_SUMMATION)
		density_str = "Summation";
	else if (GlobalParameters.DENSITY_TYPE == DENSITY_DIFFERENTIAL)
		density_str = "Differential";

	file << "gravity = {" << GlobalParameters.gravity_vector.get(0) << ", " << GlobalParameters.gravity_vector.get(1) << "}" << std::endl;
	file << "v_topwall = {" << GlobalParameters.vw_top.get(0) << ", " << GlobalParameters.vw_top.get(1) << "}" << std::endl;
	file << "v_bottomwall = {" << GlobalParameters.vw_bottom.get(0) << ", " << GlobalParameters.vw_bottom.get(1) << "}" << std::endl;
	file << "length = {" << Lx << ", " << Ly << "}" << std::endl;
	file << "Nfluid = {" << GlobalParameters.Nfluid[0] << ", " << GlobalParameters.Nfluid[1] << "}" << std::endl;
	file << "dp = " << dp << std::endl;
	file << "H = " << GlobalParameters.H << std::endl;
	file << "Mass = " << GlobalParameters.MassFluid << std::endl;
	file << "umax = " << GlobalParameters.umax << std::endl;
	file << "nu = " << GlobalParameters.nu << std::endl;
	file << "Re =  " << GlobalParameters.Re << std::endl;
	file << "rho_0 = " << GlobalParameters.rho_zero << std::endl;
	file << "gamma = " << GlobalParameters.gamma_ << std::endl;
	file << "xi = " << GlobalParameters.xi << std::endl;
	file << "c = " << GlobalParameters.cbar << std::endl;
	file << "B = " << GlobalParameters.B << std::endl;
	file << "Pbackground = " << GlobalParameters.Pbackground << std::endl;
	file << "CFLnumber = " << GlobalParameters.CFLnumber << std::endl;
	file << "SCENARIO = " << scenario_str << std::endl;
	file << "BC_TYPE = " << BC_str << std::endl;
	file << "DENSITY_TYPE = " << density_str << std::endl;
}

void AddFlatWallNewBC(particles &vd,
					  const int k0,
					  const int kmax,
					  const Point<DIM, double> Corner,
					  const Point<DIM, double> UnitOffset,
					  const double dx,
					  const Point<DIM, double> obstacle_centre,
					  const Point<DIM, double> obstacle_velocity,
					  const double obstacle_omega)
{

	for (int k = k0; k < kmax; k++)
	{
		Point<DIM, double> WallPos = Corner + k * UnitOffset;

		vd.add();
		vd.template getLastProp<type>() = BOUNDARY;
		for (int xyz = 0; xyz < DIM; xyz++)
		{
			vd.getLastPos()[xyz] = WallPos.get(xyz); //+ ((double)rand() / RAND_MAX - 0.5) * dp;
			vd.template getLastProp<velocity>()[xyz] = obstacle_velocity.get(xyz);
			vd.template getLastProp<force>()[xyz] = 0.0;
			vd.template getLastProp<force_transport>()[xyz] = obstacle_centre.get(xyz);
			vd.template getLastProp<v_transport>()[xyz] = 0.0;
			vd.template getLastProp<normal_vector>()[xyz] = 0.0;
		}

		vd.template getLastProp<pressure>() = 0.0;
		vd.template getLastProp<rho>() = 0.0;
		vd.template getLastProp<drho>() = 0.0;
		vd.template getLastProp<curvature_boundary>() = 0.0;
		vd.template getLastProp<arc_length>() = dx;
		vd.template getLastProp<vd_omega>() = obstacle_omega;
	}
}

bool isAvobeLine(Point<DIM, double> P, Point<DIM, double> Q, Point<DIM, double> EvalPoint)
{
	// P,Q two points forming a line, EvalPoint the point to check if it is above the line or not
	// with a small epsilon to avoid particles to be very close to the line
	double slope = (Q.get(1) - P.get(1)) / (Q.get(0) - P.get(0));
	// we compare y- y1 + m(x-x1)
	// =0 means the point is in the line
	// >0 means the point is above the line
	// <0 means the point is below the line
	double yline = P.get(1) + slope * (EvalPoint.get(0) - P.get(0));
	double epsilon = 0.10 * dp;
	if (EvalPoint.get(1) > yline + epsilon)
		return true;
	else
		return false;
}
bool isBelowLine(Point<DIM, double> P, Point<DIM, double> Q, Point<DIM, double> EvalPoint)
{
	// P,Q two points forming a line, EvalPoint the point to check if it is above the line or not
	// with a small epsilon to avoid particles to be very close to the line
	double slope = (Q.get(1) - P.get(1)) / (Q.get(0) - P.get(0));
	// we compare y- y1 + m(x-x1)
	// =0 means the point is in the line
	// >0 means the point is above the line
	// <0 means the point is below the line
	double yline = P.get(1) + slope * (EvalPoint.get(0) - P.get(0));
	double epsilon = 0.10 * dp;
	if (EvalPoint.get(1) < yline - epsilon)
		return true;
	else
		return false;
}
void PlaceProbes(probe_particles &probes, const int k0, const int kmax,
				 const Point<DIM, double> Corner, const Point<DIM, double> UnitOffset)
{
	for (int k = k0; k < kmax; k++)
	{
		Point<DIM, double> probe_position = Corner + k * UnitOffset;
		probes.add();
		for (int xyz = 0; xyz < DIM; xyz++)
		{
			probes.getLastPos()[xyz] = probe_position.get(xyz);
		}
		probes.template getLastProp<0>() = 0.0;
	}
}
class Obstacle
{
public:
	Point<DIM, double> Centre_;
	Point<DIM, double> LinearVelocity_;
	double AngularVelocity_;
	double refine_factor;

public:
	Obstacle(Point<DIM, double> Centre, Point<DIM, double> vel = {0.0, 0.0}, double omega = 0.0, double rf = 1.0) : Centre_(Centre), LinearVelocity_(vel), AngularVelocity_(omega), refine_factor(rf) {}
	virtual bool isInside(Point<DIM, double> P) = 0;
	virtual void AddObstacle(particles &vd) = 0;
};

class EmptyObstacle : public Obstacle
{
public:
	EmptyObstacle() : Obstacle({0.0, 0.0}, {0.0, 0.0}, 0.0, 1.0) {}
	bool isInside(Point<DIM, double> P) override
	{
		return false;
	}
	void AddObstacle(particles &vd) override
	{
	}
};

class CylinderObstacle : public Obstacle
{
private:
	double Radius_;
	Sphere<DIM, double> Cylinder_;

public:
	CylinderObstacle(double Radius, Point<DIM, double> centre, Point<DIM, double> vel = {0.0, 0.0}, double omega = 0.0, double rf = 1.0) : Obstacle(centre, vel, omega, rf), Radius_(Radius), Cylinder_(centre, Radius) {}

	bool isInside(Point<DIM, double> P) override
	{
		double radius_aux = Radius_;

		if (GlobalParameters.BC_TYPE == NEW_NO_SLIP)
			radius_aux = Radius_ + 0.1 * dp;

		Sphere<DIM, double> Cylinderaux(Centre_, radius_aux);
		return Cylinderaux.isInside(P);
	}
	bool isOutside(Point<DIM, double> P) // for outer cylinder in taylor couette
	{
		double radius_aux = Radius_;
		if (GlobalParameters.BC_TYPE == NEW_NO_SLIP)
			radius_aux = Radius_ - 0.1 * dp;
		Sphere<DIM, double> Cylinderaux(Centre_, radius_aux);
		return Cylinderaux.isInside(P);
	}

	void AddObstacle(particles &vd)
	{
		const double dx_self = dp / refine_factor;
		const double perimeter = 2.0 * M_PI * Radius_;
		const int Np_cylinder = ceil(perimeter / dx_self);
		const double dtheta = 2.0 * M_PI / Np_cylinder;
		const double dxwall = dtheta * Radius_;
		double theta = 0.0;

		Point<DIM, double> Cylinder_particle;

		for (int k = 0; k < Np_cylinder; k++)
		{

			Cylinder_particle[0] = Centre_.get(0) + Radius_ * cos(theta);
			Cylinder_particle[1] = Centre_.get(1) + Radius_ * sin(theta);

			if (DIM == 3)
			{
				Cylinder_particle[2] = Centre_.get(2);
			}

			vd.add();
			vd.template getLastProp<type>() = BOUNDARY;
			for (int xyz = 0; xyz < DIM; xyz++)
			{
				vd.getLastPos()[xyz] = Cylinder_particle.get(xyz); // + ((double)rand() / RAND_MAX - 0.5) * dp;
				vd.template getLastProp<force>()[xyz] = 0.0;
				vd.template getLastProp<force_transport>()[xyz] = Centre_.get(xyz);
				vd.template getLastProp<v_transport>()[xyz] = 0.0;
				vd.template getLastProp<normal_vector>()[xyz] = 0.0;
				vd.template getLastProp<velocity>()[xyz] = LinearVelocity_.get(xyz);
			}

			vd.template getLastProp<pressure>() = 0.0;
			vd.template getLastProp<rho>() = 0.0;
			vd.template getLastProp<drho>() = 0.0;
			vd.template getLastProp<curvature_boundary>() = 0.0; // 1.0 / Radius_;
			vd.template getLastProp<arc_length>() = dxwall;
			vd.template getLastProp<vd_omega>() = AngularVelocity_;
			theta += dtheta;
		}
	}
};

class RectangleObstacle : public Obstacle
{
private:
	const unsigned int BaseLength_;
	const unsigned int HeigthLength_;
	Box<DIM, double> Rectangle_;

	Point<DIM, double> LowerLeft_;
	Point<DIM, double> LowerRight_;
	Point<DIM, double> UpperLeft_;
	Point<DIM, double> UpperRight_;

public:
	RectangleObstacle(Point<DIM, double> centre, unsigned int BaseLength, unsigned int HeigthLength, Point<DIM, double> vel = {0.0, 0.0}, double omega = 0.0, double rf = 1.0)
		: Obstacle(centre, vel, omega, rf), BaseLength_(BaseLength), HeigthLength_(HeigthLength)
	{
		LowerLeft_ = Point<DIM, double>{Centre_.get(0) - (((double)BaseLength_ - 1.0) * dp) / 2.0,
										Centre_.get(1) - (((double)HeigthLength_ - 1.0) * dp) / 2.0};
		LowerRight_ = Point<DIM, double>{Centre_.get(0) + (((double)BaseLength_ - 1.0) * dp) / 2.0,
										 Centre_.get(1) - (((double)HeigthLength_ - 1.0) * dp) / 2.0};
		UpperLeft_ = Point<DIM, double>{Centre_.get(0) - (((double)BaseLength_ - 1.0) * dp) / 2.0,
										Centre_.get(1) + (((double)HeigthLength_ - 1.0) * dp) / 2.0};
		UpperRight_ = Point<DIM, double>{Centre_.get(0) + (((double)BaseLength_ - 1.0) * dp) / 2.0,
										 Centre_.get(1) + (((double)HeigthLength_ - 1.0) * dp) / 2.0};

		Rectangle_ = Box<DIM, double>(LowerLeft_, UpperRight_);
	}
	bool isInside(Point<DIM, double> P) override
	{
		return Rectangle_.isInside(P);
	}

	void AddObstacle(particles &vd)
	{
		const double baseLength = (BaseLength_ - 1) * dp;
		const double heigthLength = (HeigthLength_ - 1) * dp;

		Point<DIM, double> Xoffset = {dp, 0.0};
		Point<DIM, double> Yoffset = {0.0, dp};
		// Lower wall
		AddFlatWallNewBC(vd, 0, BaseLength_, LowerLeft_, Xoffset, dp, Centre_, LinearVelocity_, AngularVelocity_);
		// Right wall
		AddFlatWallNewBC(vd, 1, HeigthLength_ - 1, LowerRight_, Yoffset, dp, Centre_, LinearVelocity_, AngularVelocity_);
		// Upper wall
		AddFlatWallNewBC(vd, 0, BaseLength_, UpperRight_, -1.0 * Xoffset, dp, Centre_, LinearVelocity_, AngularVelocity_);
		// Left wall
		AddFlatWallNewBC(vd, 1, HeigthLength_ - 1, UpperLeft_, -1.0 * Yoffset, dp, Centre_, LinearVelocity_, AngularVelocity_);
	}
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
	TriangleObstacle(Point<DIM, double> centre, double BaseLength, double HeigthLength, Point<DIM, double> vel = {0.0, 0.0}, double omega = 0.0, double rf = 1.0) : Obstacle(centre, vel, omega, rf), BaseLength_(BaseLength), HeigthLength_(HeigthLength)
	{
		LowerLeft_ = Point<DIM, double>{Centre_.get(0) - 2.0 * BaseLength_ / 3.0,
										Centre_.get(1) - HeigthLength_ / 3.0};

		LowerRight_ = Point<DIM, double>{Centre_.get(0) + BaseLength_ / 3.0,
										 Centre_.get(1) - HeigthLength_ / 3.0};

		UpperRight_ = Point<DIM, double>{Centre_.get(0) + BaseLength_ / 3.0,
										 Centre_.get(1) + 2.0 * HeigthLength_ / 3.0};

		ContainingRectangle_ = Box<DIM, double>(LowerLeft_, UpperRight_);
	}

	bool isInside(Point<DIM, double> P) override
	{
		return (ContainingRectangle_.isInside(P) && !isAvobeLine(LowerLeft_, UpperRight_, P));
	}

	void AddObstacle(particles &vd)
	{
		const double dx_self = dp / refine_factor;
		// Lower wall
		const int N_bottom = ceil(BaseLength_ / dx_self);
		const double dxwall_bottom = BaseLength_ / N_bottom;
		const Point<DIM, double> Xoffset = {dxwall_bottom, 0.0};
		AddFlatWallNewBC(vd, 0, N_bottom + 1, LowerLeft_, Xoffset, dxwall_bottom, Centre_, LinearVelocity_, AngularVelocity_);

		// Right wall
		const int N_right = ceil(HeigthLength_ / dx_self);
		const double dxwall_right = HeigthLength_ / N_right;
		const Point<DIM, double> Yoffset = {0.0, dxwall_right};
		AddFlatWallNewBC(vd, 1, N_right + 1, LowerRight_, Yoffset, dxwall_right, Centre_, LinearVelocity_, AngularVelocity_);

		//  Hypothenuse wall
		// We want particles spaced roughly by dp
		const double HypothenuseLength = sqrt(BaseLength_ * BaseLength_ + HeigthLength_ * HeigthLength_);
		const int Ndiag = ceil(HypothenuseLength / dx_self);  // integer number of particles that can fit in the diagonal
		const double dxwall_diag = HypothenuseLength / Ndiag; // actual spacing between particles ( close to dp but not exactly)
		const double sin_theta = HeigthLength_ / HypothenuseLength;
		const double cos_theta = BaseLength_ / HypothenuseLength;
		const Point<DIM, double> Diagoffset{dxwall_diag * cos_theta, dxwall_diag * sin_theta};

		AddFlatWallNewBC(vd, 1, Ndiag, UpperRight_, -1.0 * Diagoffset, dxwall_diag, Centre_, LinearVelocity_, AngularVelocity_);
	}
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
	TriangleEqui(Point<DIM, double> centre, double sidelength, Point<DIM, double> vel = {0.0, 0.0}, double omega = 0.0, double rf = 1.0) : Obstacle(centre, vel, omega, rf), SideLength_(sidelength)
	{

		UpperRight_ = Point<DIM, double>{Centre_.get(0) + sqrt(3) * SideLength_ / 6.0,
										 Centre_.get(1) + SideLength_ / 2.0};
		LowerRight_ = Point<DIM, double>{Centre_.get(0) + sqrt(3) * SideLength_ / 6.0,
										 Centre_.get(1) - SideLength_ / 2.0};
		TriangleTip_ = Point<DIM, double>{Centre_.get(0) - sqrt(3.0) * SideLength_ / 3.0,
										  Centre_.get(1)};

		LowerLeft_ = Point<DIM, double>{TriangleTip_.get(0),
										LowerRight_.get(1)};

		ContainingRectangle_ = Box<DIM, double>(LowerLeft_, UpperRight_);
	}
	bool isInside(Point<DIM, double> P) override
	{
		return (ContainingRectangle_.isInside(P) && !isAvobeLine(TriangleTip_, UpperRight_, P) && !isBelowLine(TriangleTip_, LowerRight_, P));
	}

	void AddObstacle(particles &vd)
	{
		double dx_self = dp / refine_factor;
		const int N_wall = ceil(SideLength_ / dx_self);
		const double dxwall = SideLength_ / N_wall;
		Point<DIM, double> Yoffset = {0.0, dxwall};

		// Right wall
		AddFlatWallNewBC(vd, 0, N_wall + 1, LowerRight_, Yoffset, dxwall, Centre_, LinearVelocity_, AngularVelocity_);

		// Inclined walls

		const double cos_theta = sqrt(3.0) / 2.0;
		const double sin_theta = 1 / 2.0;
		const Point<DIM, double> Diagoffset{dxwall * cos_theta, dxwall * sin_theta};
		Point<DIM, double> DiagoffsetNW = Diagoffset;
		DiagoffsetNW.get(0) = -1.0 * DiagoffsetNW.get(0);
		//  Hypothenuse upper wall
		AddFlatWallNewBC(vd, 1, N_wall + 1, UpperRight_, -1.0 * Diagoffset, dxwall, Centre_, LinearVelocity_, AngularVelocity_);
		// Hypothenuse lower wall
		AddFlatWallNewBC(vd, 1, N_wall, LowerRight_, DiagoffsetNW, dxwall, Centre_, LinearVelocity_, AngularVelocity_);
	}
};

void interact_probe_boundary_new(particles &vd,
								 vect_dist_key_dx probekey,
								 const Point<DIM, double> &r_wall_to_probe,
								 unsigned long &boundary_key,
								 const int component,
								 double &W_sum,
								 double &magnitude_tmp)
{
	Point<DIM, double> r_probe_to_wall = -1.0 * r_wall_to_probe;

	Point<DIM, double> normal = vd.getProp<normal_vector>(boundary_key);
	// get curvature, and arc length
	double kappa = vd.getProp<curvature_boundary>(boundary_key);
	double dxwall = vd.getProp<arc_length>(boundary_key);

	std::array<Point<DIM, double>, 3> r_boundary = GetBoundaryPositions(r_probe_to_wall, normal);
	std::array<double, 3> r_boundary_norm = {getVectorNorm(r_boundary[0]), getVectorNorm(r_boundary[1]), getVectorNorm(r_boundary[2])};
	std::array<double, 3> Volume_boundary;
	std::array<double, 3> Mass_boundary;

	for (int i = 0; i < 3; i++)
	{
		if (r_boundary_norm[i] < GlobalParameters.r_threshold)
		{
			// compute volume of boundary particle, this gives density and pressure
			Volume_boundary[i] = 0.5 * (2.0 * dp + dp * dp * kappa - 2.0 * (i + 1.0) * dp * dp * kappa) * dxwall;
			double ker = Wab(r_boundary_norm[i]) * Volume_boundary[i];
			W_sum += ker;
			magnitude_tmp += vd.getProp<velocity>(boundary_key)[component] * ker;
		}
	}
}

template <typename CellList>
inline void sensor_velocity_comp(particles &vd,
								 probe_particles &probes,
								 Vcluster<> &v_cl,
								 CellList &NN,
								 const int component,
								 Obstacle *obstacle_ptr)
{

	auto it = probes.getDomainIterator();
	while (it.isNext()) // iterate all probes
	{

		auto probekey = it.get();
		Point<DIM, double> xp = probes.getPos(probekey);

		double magnitude_tmp = 0.0;
		double W_sum = 0.0;

		if (obstacle_ptr->isInside(xp)) // if probe is inside the obstacle it gets a 0 value
		{
			probes.getProp<0>(probekey) = 0.0;
			++it;
			continue;
		}
		else
		{
			// get the iterator over the neighbohood particles of the probes position
			auto itg = NN.getNNIterator(NN.getCell(xp));
			while (itg.isNext())
			{
				// get key of the particle
				auto q = itg.get();

				// Get the position of the neighborhood particle q
				Point<DIM, double> xq = vd.getPos(q);
				Point<DIM, double> dr = xp - xq;
				// Calculate the distance
				double r = getVectorNorm(dr);

				if (vd.template getProp<type>(q) == BOUNDARY)
				{
					if (r < 0.01 * dp) // if probe is placed on top of a wall particle
					{
						W_sum = 1.0;
						magnitude_tmp = vd.template getProp<velocity>(q)[component];
						break;
					}

					// if (BC_TYPE == NEW_NO_SLIP)
					// {
					// 	// interact_probe_boundary_new(vd, probekey, dr, q, component, W_sum, magnitude_tmp);
					// }
					else if (GlobalParameters.BC_TYPE == NO_SLIP)
					{
						double ker = Wab(r) * (GlobalParameters.MassBound / GlobalParameters.rho_zero);
						W_sum += ker;
						magnitude_tmp += vd.template getProp<v_transport>(q)[component] * ker;
					}
				}
				else if (vd.template getProp<type>(q) == FLUID)
				{
					double rhoq = vd.template getProp<rho>(q);
					double ker = Wab(r) * (GlobalParameters.MassFluid / rhoq);

					// Also keep track of the calculation of the summed kernel
					W_sum += ker;

					// Add the total magnitude contribution
					magnitude_tmp += vd.template getProp<velocity>(q)[component] * ker;
				}
				// next neighborhood particle
				++itg;
			}

			// We calculate the magnitude normalizing with W_sum
			if (W_sum == 0.0)
				magnitude_tmp = 0.0;
			else
				magnitude_tmp = 1.0 / W_sum * magnitude_tmp;

			probes.getProp<0>(probekey) = magnitude_tmp;
			++it;
		}
	}
}

void computeAverageVelocity(particles &vd, Vcluster<> &v_cl, double t, std::ofstream &avgvelstream, double cylinder_force)
{
	auto it = vd.getDomainIterator();
	double vx{0.0}, vy{0.0}, v{0.0};
	int counter = 0;
	while (it.isNext())
	{
		auto key = it.get();

		if (vd.template getProp<type>(key) == FLUID)
		{
			vx += vd.template getProp<velocity>(key)[0];
			vy += vd.template getProp<velocity>(key)[1];
			v += getVectorNorm(vd.template getProp<velocity>(key));
			counter++;
		}

		++it;
	}

	v_cl.sum(vx);
	v_cl.sum(vy);
	v_cl.sum(v);
	v_cl.sum(counter);
	v_cl.execute();
	vx = vx / counter;
	vy = vy / counter;
	v = v / counter;
	cylinder_force = cylinder_force / (GlobalParameters.eta * v);

	if (v_cl.getProcessUnitID() == 0)
	{
		avgvelstream << t << ", " << v << ", " << cylinder_force << std::endl;
	}
}

void ParseXMLFile(const std::string &filename)
{
	tinyxml2::XMLDocument doc;
	tinyxml2::XMLError eResult = doc.LoadFile(filename.c_str());

	if (eResult != tinyxml2::XML_SUCCESS)
	{
		std::cerr << "Error: Unable to load XML file: " << filename << std::endl;
		return;
	}

	// Parse <simulation> element
	tinyxml2::XMLElement *simulationElement = doc.FirstChildElement("configuration")->FirstChildElement("simulation");
	if (simulationElement)
	{
		const char *scenario_str = simulationElement->Attribute("scenario");
		const char *bc_type_str = simulationElement->Attribute("bcType");
		const char *density_type_Str = simulationElement->Attribute("densityType");
		const char *writer_str = simulationElement->Attribute("writerType");

		if (strcmp(scenario_str, "Poiseuille") == 0)
			GlobalParameters.SCENARIO = POISEUILLE;
		else if (strcmp(scenario_str, "Couette") == 0)
			GlobalParameters.SCENARIO = COUETTE;
		else if (strcmp(scenario_str, "Hydrostatic") == 0)
			GlobalParameters.SCENARIO = HYDROSTATIC;
		else if (strcmp(scenario_str, "CylinderArray") == 0)
			GlobalParameters.SCENARIO = CYLINDER_ARRAY;
		else if (strcmp(scenario_str, "CylinderLattice") == 0)
			GlobalParameters.SCENARIO = CYLINDER_LATTICE;
		else if (strcmp(scenario_str, "Square") == 0)
			GlobalParameters.SCENARIO = SQUARE;
		else if (strcmp(scenario_str, "Triangle") == 0)
			GlobalParameters.SCENARIO = TRIANGLE;
		else if (strcmp(scenario_str, "TriangleEquilateral") == 0)
			GlobalParameters.SCENARIO = TRIANGLE_EQUILATERAL;
		else if (strcmp(scenario_str, "Cavity") == 0)
			GlobalParameters.SCENARIO = CAVITY;
		else if (strcmp(scenario_str, "Step") == 0)
			GlobalParameters.SCENARIO = STEP;
		else if (strcmp(scenario_str, "TaylorCouette") == 0)
			GlobalParameters.SCENARIO = TAYLOR_COUETTE;
		else if (strcmp(scenario_str, "MovingObstacle") == 0)
			GlobalParameters.SCENARIO = MOVING_OBSTACLE;
		else
		{
			throw std::runtime_error("Unknown scenario");
		}

		if (strcmp(bc_type_str, "old") == 0)
			GlobalParameters.BC_TYPE = NO_SLIP;
		else if (strcmp(bc_type_str, "new") == 0)
			GlobalParameters.BC_TYPE = NEW_NO_SLIP;

		if (strcmp(density_type_Str, "Summation") == 0)
			GlobalParameters.DENSITY_TYPE = DENSITY_SUMMATION;
		else if (strcmp(density_type_Str, "Differential") == 0)
			GlobalParameters.DENSITY_TYPE = DENSITY_DIFFERENTIAL;

		if (strcmp(writer_str, "VTK") == 0)
			GlobalParameters.WRITER = VTK_WRITER;
		else if (strcmp(writer_str, "CSV") == 0)
			GlobalParameters.WRITER = CSV_WRITER;

		simulationElement->QueryDoubleAttribute("time", &(GlobalParameters.t_end));
		simulationElement->QueryDoubleAttribute("write_const", &(GlobalParameters.write_const));
		const char *custom_str = simulationElement->Attribute("custom_string");

		GlobalParameters.custom_string += custom_str;
	}

	// Parse <geometry> element
	tinyxml2::XMLElement *geometryElement = doc.FirstChildElement("configuration")->FirstChildElement("geometry");
	if (geometryElement)
	{
		geometryElement->QueryDoubleAttribute("lengthScale", &GlobalParameters.LengthScale);
		geometryElement->QueryDoubleAttribute("rf", &GlobalParameters.rf);
		geometryElement->QueryAttribute("sizeX", &(GlobalParameters.Nfluid[0]));
		geometryElement->QueryAttribute("sizeY", &(GlobalParameters.Nfluid[1]));
		if (DIM == 3)
			geometryElement->QueryAttribute("sizeZ", &GlobalParameters.Nfluid[2]);
	}

	// Parse <physics> element
	tinyxml2::XMLElement *physicsElement = doc.FirstChildElement("configuration")->FirstChildElement("physics");
	if (physicsElement)
	{

		physicsElement->QueryDoubleAttribute("rho0", &GlobalParameters.rho_zero);
		physicsElement->QueryDoubleAttribute("nu", &GlobalParameters.nu);
		physicsElement->QueryDoubleAttribute("Bfactor", &GlobalParameters.Bfactor);
		physicsElement->QueryDoubleAttribute("gravityX", &(GlobalParameters.gravity_vector.get(0)));
		physicsElement->QueryDoubleAttribute("gravityY", &(GlobalParameters.gravity_vector.get(1)));
		if (DIM == 3)
			physicsElement->QueryDoubleAttribute("gravityZ", &(GlobalParameters.gravity_vector.get(2)));
		physicsElement->QueryDoubleAttribute("vTopX", &(GlobalParameters.vw_top.get(0)));
		physicsElement->QueryDoubleAttribute("vTopY", &(GlobalParameters.vw_top.get(1)));
		if (DIM == 3)
			physicsElement->QueryDoubleAttribute("vTopZ", &(GlobalParameters.vw_top.get(2)));

		physicsElement->QueryDoubleAttribute("vBotX", &(GlobalParameters.vw_bottom.get(0)));
		physicsElement->QueryDoubleAttribute("vBotY", &(GlobalParameters.vw_bottom.get(1)));
		if (DIM == 3)
			physicsElement->QueryDoubleAttribute("vBotZ", &(GlobalParameters.vw_bottom.get(2)));
	}

	// Parse <obstacle> element
	tinyxml2::XMLElement *obstacleElement = doc.FirstChildElement("configuration")->FirstChildElement("obstacle");
	if (obstacleElement)
	{
		obstacleElement->QueryDoubleAttribute("velX", &(GlobalParameters.ObstacleVelocity.get(0)));
		obstacleElement->QueryDoubleAttribute("velY", &(GlobalParameters.ObstacleVelocity.get(1)));
		if (DIM == 3)
			obstacleElement->QueryDoubleAttribute("velZ", &(GlobalParameters.ObstacleVelocity.get(2)));

		obstacleElement->QueryDoubleAttribute("omega", &GlobalParameters.ObstacleOmega);
	}

	// Parse <TaylorCouette> element
	tinyxml2::XMLElement *taylorCouetteElement = doc.FirstChildElement("configuration")->FirstChildElement("TaylorCouette");
	if (taylorCouetteElement)
	{
		taylorCouetteElement->QueryDoubleAttribute("Rin", &GlobalParameters.Rin);
		taylorCouetteElement->QueryDoubleAttribute("Rout", &GlobalParameters.Rout);
		taylorCouetteElement->QueryDoubleAttribute("Win", &GlobalParameters.Win);
		taylorCouetteElement->QueryDoubleAttribute("Wout", &GlobalParameters.Wout);
	}
}
void InitializeConstants()
{
	// Given the scenario and xml parameters, we set the constants of GlobalParameters

	// Set boundary conditions periodic or non periodic
	size_t Nbound = (GlobalParameters.BC_TYPE == NEW_NO_SLIP) ? 1 : 3;

	if (GlobalParameters.SCENARIO == POISEUILLE ||
		GlobalParameters.SCENARIO == COUETTE ||
		GlobalParameters.SCENARIO == CYLINDER_ARRAY ||
		GlobalParameters.SCENARIO == SQUARE ||
		GlobalParameters.SCENARIO == TRIANGLE ||
		GlobalParameters.SCENARIO == TRIANGLE_EQUILATERAL)
	{
		GlobalParameters.Nboundary[0] = 0;
		GlobalParameters.Nboundary[1] = Nbound;
		GlobalParameters.bc[0] = PERIODIC;
		GlobalParameters.bc[1] = NON_PERIODIC;
	}
	else if (GlobalParameters.SCENARIO == HYDROSTATIC ||
			 GlobalParameters.SCENARIO == CAVITY ||
			 GlobalParameters.SCENARIO == MOVING_OBSTACLE ||
			 GlobalParameters.SCENARIO == TAYLOR_COUETTE)
	{
		GlobalParameters.Nboundary[0] = Nbound;
		GlobalParameters.Nboundary[1] = Nbound;
		GlobalParameters.bc[0] = NON_PERIODIC;
		GlobalParameters.bc[1] = NON_PERIODIC;
	}
	else if (GlobalParameters.SCENARIO == CYLINDER_LATTICE)
	{
		GlobalParameters.Nboundary[0] = 0;
		GlobalParameters.Nboundary[1] = 0;
		GlobalParameters.bc[0] = PERIODIC;
		GlobalParameters.bc[1] = PERIODIC;
	}

	// Set particle spacing, definition depends on scenario
	if (GlobalParameters.SCENARIO != CYLINDER_ARRAY && GlobalParameters.SCENARIO != CYLINDER_LATTICE)
	{
		dp = GlobalParameters.LengthScale / GlobalParameters.Nfluid[1];
	}
	else
	{
		if (GlobalParameters.SCENARIO == CYLINDER_ARRAY) // chanel height is 4 times the cylinder radius
			dp = 4.0 * GlobalParameters.LengthScale / (double)GlobalParameters.Nfluid[1];
		else if (GlobalParameters.SCENARIO == CYLINDER_LATTICE) // chanel height is 5 times the cylinder radius
			dp = 5.0 * GlobalParameters.LengthScale / (double)GlobalParameters.Nfluid[1];
	}

	GlobalParameters.length[0] = dp * (GlobalParameters.Nfluid[0] - GlobalParameters.bc[0]);
	GlobalParameters.length[1] = dp * (GlobalParameters.Nfluid[1] - GlobalParameters.bc[1]);

	// Set maximum velocity
	if (GlobalParameters.SCENARIO == POISEUILLE)
		GlobalParameters.umax = GlobalParameters.gravity_vector.get(0) * GlobalParameters.LengthScale * GlobalParameters.LengthScale / (8.0 * GlobalParameters.nu);
	else if (GlobalParameters.SCENARIO == COUETTE)
		GlobalParameters.umax = GlobalParameters.vw_top.get(0);
	else if (GlobalParameters.SCENARIO == HYDROSTATIC)
		GlobalParameters.umax = 0.1;
	else if (GlobalParameters.SCENARIO == CYLINDER_ARRAY)
		GlobalParameters.umax = 2.5e-3;
	else if (GlobalParameters.SCENARIO == CYLINDER_LATTICE)
		GlobalParameters.umax = 1.2e-4; // umax = 5.77 * 1e-5; // (morris, to get c=5.77*1e-4)
	else if (GlobalParameters.SCENARIO == SQUARE)
		GlobalParameters.umax = 4.1 * 1e-1;
	else if (GlobalParameters.SCENARIO == TRIANGLE)
		GlobalParameters.umax = 3.5 * 1e-1; // from previous simulations for nu = 0.01
	else if (GlobalParameters.SCENARIO == TRIANGLE_EQUILATERAL)
		GlobalParameters.umax = 4e-1; // from previous simulations for nu = 0.01
	else if (GlobalParameters.SCENARIO == CAVITY)
		GlobalParameters.umax = GlobalParameters.vw_top.get(0);
	else if (GlobalParameters.SCENARIO == MOVING_OBSTACLE)
		GlobalParameters.umax = 1.5; // from previous simulations for nu = 0.01
	else if (GlobalParameters.SCENARIO == TAYLOR_COUETTE)
		GlobalParameters.umax = 4.0 * abs(GlobalParameters.Wout * GlobalParameters.Rout);

	// Enable probes in some scenarios
	if (GlobalParameters.SCENARIO == CAVITY || GlobalParameters.SCENARIO == CYLINDER_LATTICE)
	{
		GlobalParameters.PROBES_ENABLED = 1;
	}
	// Set general parameters
	GlobalParameters.H = GlobalParameters.Hconst * dp;
	GlobalParameters.r_threshold = 3.0 * GlobalParameters.H;
	GlobalParameters.Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / GlobalParameters.H / GlobalParameters.H / GlobalParameters.H : 7.0 / 478.0 / M_PI / GlobalParameters.H / GlobalParameters.H;
	GlobalParameters.MassFluid = GlobalParameters.rho_zero * (DIM == 3 ? dp * dp * dp : dp * dp);
	GlobalParameters.MassBound = GlobalParameters.rho_zero * (DIM == 3 ? dp * dp * dp : dp * dp);
	GlobalParameters.cbar = GlobalParameters.coeff_sound * GlobalParameters.umax;
	GlobalParameters.B = GlobalParameters.rho_zero * GlobalParameters.cbar * GlobalParameters.cbar / GlobalParameters.gamma_;
	GlobalParameters.Pbackground = GlobalParameters.Bfactor * GlobalParameters.B;
	GlobalParameters.eta = GlobalParameters.nu * GlobalParameters.rho_zero;
	GlobalParameters.Re = GlobalParameters.umax * GlobalParameters.LengthScale / GlobalParameters.nu;
	GlobalParameters.gravity = getVectorNorm(GlobalParameters.gravity_vector);
	GlobalParameters.ObstacleCenter = {(GlobalParameters.length[0] + GlobalParameters.bc[0] * dp) / 2.0, (GlobalParameters.length[1] + GlobalParameters.bc[1] * dp) / 2.0};
}
void CreateParticleGeometry(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Vcluster<> &v_cl, Obstacle *&obstacle_ptr)
{

	// Non periodic situation grid of 5 fluid particles and 3 boundary particles
	// We need a virtual grid of 5 + 2*(3+1) particles,
	// therefore the domain is discretized with 13 grid points,
	// when we use DrawParticles::DrawBox we will draw only the particles at the grid positons strictly inside the box,
	// the () repesent the recipient box, and the || represent the fluid box, we can see how this distribution places exactly 5 fluid particles inside and 3 boundary particles
	//           D-(-o--o--o-|-x--x--x--x--x--|-o-o-o-)-D
	// D: domain, o: boundary, x: fluid, --: dp distance
	// in a periodic situation we have the following
	// .....--x--x--D-|-x--x--x--x--x--|-D--x--x--......
	// therefore we need a grid of 5 + 2 particles, and the domain is discretized with 7 grid points

	// Size of the virtual cartesian grid that defines where to place the particles
	size_t sz[DIM];

	// In the case of the new bc we need particles exactly at the wall, for this we need sz_aux
	// We want to put one virtual grid point between each pair of the old ones,
	// so that the new spacing is dp/2, and we can put a marker particle exactly at the wall
	size_t sz_aux[DIM];

	obstacle_ptr = new EmptyObstacle();

	// Initialize obstacle in scenarios where needed
	if (GlobalParameters.SCENARIO == CYLINDER_ARRAY)
		obstacle_ptr = new CylinderObstacle(GlobalParameters.LengthScale, GlobalParameters.ObstacleCenter, GlobalParameters.ObstacleVelocity, GlobalParameters.ObstacleOmega, GlobalParameters.rf);
	else if (GlobalParameters.SCENARIO == CYLINDER_LATTICE)
		obstacle_ptr = new CylinderObstacle(GlobalParameters.LengthScale, GlobalParameters.ObstacleCenter, GlobalParameters.ObstacleVelocity, GlobalParameters.ObstacleOmega, GlobalParameters.rf);
	else if (GlobalParameters.SCENARIO == SQUARE)
	{
		int integerBaseLength = 11;
		int integerHeigthLength = 11;
		obstacle_ptr = new RectangleObstacle(GlobalParameters.ObstacleCenter, integerBaseLength, integerHeigthLength, GlobalParameters.ObstacleVelocity, GlobalParameters.ObstacleOmega, GlobalParameters.rf);
	}
	else if (GlobalParameters.SCENARIO == TRIANGLE)
	{
		const double BaseLength = 0.47;
		const double HeigthLength = 0.3;
		obstacle_ptr = new TriangleObstacle(GlobalParameters.ObstacleCenter, BaseLength, HeigthLength, GlobalParameters.ObstacleVelocity, GlobalParameters.ObstacleOmega, GlobalParameters.rf);
	}
	else if (GlobalParameters.SCENARIO == TRIANGLE_EQUILATERAL)
		obstacle_ptr = new TriangleEqui(GlobalParameters.ObstacleCenter, GlobalParameters.LengthScale / 4.0, GlobalParameters.ObstacleVelocity, GlobalParameters.ObstacleOmega, GlobalParameters.rf);
	else if (GlobalParameters.SCENARIO == MOVING_OBSTACLE)
		obstacle_ptr = new TriangleEqui(GlobalParameters.ObstacleCenter, GlobalParameters.LengthScale / 4.0, GlobalParameters.ObstacleVelocity, GlobalParameters.ObstacleOmega, GlobalParameters.rf);

	double refine_factor = GlobalParameters.rf;

	// Now define the iterator boxes
	// We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
	double offset_domain_left[DIM] = {0.0};
	double offset_domain_right[DIM] = {0.0};
	double offset_recipient[DIM] = {0.0};
	double offset_periodic_fluid[DIM] = {0.0};
	double offset_periodic_recipient[DIM] = {0.0};

	for (int xyz = 0; xyz < DIM; xyz++)
	{
		if (GlobalParameters.bc[xyz] == NON_PERIODIC) // non periodic, fluid covered by boundary
		{
			sz[xyz] = GlobalParameters.Nfluid[xyz] + 2 * (GlobalParameters.Nboundary[xyz] + 1);
			offset_domain_left[xyz] = (0.5 + GlobalParameters.Nboundary[xyz]) * dp;
			offset_domain_right[xyz] = (0.5 + GlobalParameters.Nboundary[xyz]) * dp;

			if (GlobalParameters.Nboundary[xyz] != 0)
				sz_aux[xyz] = 2 * GlobalParameters.Nfluid[xyz] - 1 + 2 * (2 * GlobalParameters.Nboundary[xyz] + 1 + 1);
			else // for a direction with no boundary particles we dont need to add anything
				sz_aux[xyz] = sz[xyz];

			if (GlobalParameters.BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
				offset_recipient[xyz] = 0.1 * GlobalParameters.Nboundary[xyz] * dp;
			else if (GlobalParameters.BC_TYPE == NO_SLIP)
				offset_recipient[xyz] = GlobalParameters.Nboundary[xyz] * dp;
		}
		else // periodic, open ended
		{
			GlobalParameters.Nfluid[xyz] -= 1;

			sz[xyz] = GlobalParameters.Nfluid[xyz] + 2;
			sz_aux[xyz] = sz[xyz];

			offset_domain_left[xyz] = 0.0;
			offset_domain_right[xyz] = dp;
			offset_periodic_fluid[xyz] = 0.75 * dp;
			offset_periodic_recipient[xyz] = 0.85 * dp;

			// correct the number of particles in case of periodicity, we substracted 1 before to accomodate the periodic boundary
			GlobalParameters.Nfluid[xyz] += 1;
		}
	}

	// Define the boxes
	Box<DIM, double> domain({-offset_domain_left[0],
							 -offset_domain_left[1]},
							{GlobalParameters.length[0] + offset_domain_right[0],
							 GlobalParameters.length[1] + offset_domain_right[1]});

	Box<DIM, double> fluid_box({0.0,
								0.0},
							   {GlobalParameters.length[0] + offset_periodic_fluid[0],
								GlobalParameters.length[1] + offset_periodic_fluid[1]});

	Box<DIM, double> recipient({-offset_recipient[0],
								-offset_recipient[1]},
							   {GlobalParameters.length[0] + offset_recipient[0] + offset_periodic_recipient[0],
								GlobalParameters.length[1] + offset_recipient[1] + offset_periodic_recipient[1]});

	// Will only be used in the new bc
	Box<DIM, double> recipient_hole({offset_recipient[0],
									 offset_recipient[1]},
									{GlobalParameters.length[0] - offset_recipient[0] + offset_periodic_fluid[0],
									 GlobalParameters.length[1] - offset_recipient[1] + offset_periodic_fluid[1]});

	for (int xyz = 0; xyz < DIM; xyz++) // correct length in periodic case
	{
		if (GlobalParameters.bc[xyz] == PERIODIC)
			GlobalParameters.length[xyz] += dp;
	}

	// extended boundary around the domain, and the processor domain
	Ghost<DIM, double> g(GlobalParameters.r_threshold);

	// create particle object
	particles vd_loc(0, domain, GlobalParameters.bc, g, DEC_GRAN(512));
	// vd is argument passed as reference we want to fill with particles
	vd = vd_loc;

	// Write constants on file
	WriteParameters(v_cl);

	// place probes
	if (GlobalParameters.PROBES_ENABLED)
	{
		// we want to place probes  in a vertical line at this locations
		Point<DIM, double> EndChannel = {GlobalParameters.length[0], 0.0};
		Point<DIM, double> HalfChannel = {GlobalParameters.length[0] / 2.0, 0.0};
		Point<DIM, double> HalfHeight = {0.0, GlobalParameters.length[1] / 2.0};
		Point<DIM, double> VerticalOffset = {0.0, dp};
		Point<DIM, double> HorizontalOffset = {dp, 0.0};
		int k0 = 0;
		int kendHeight = GlobalParameters.Nfluid[1] + 1;
		int kendWidth = GlobalParameters.Nfluid[0] + 1;

		std::vector<Point<DIM, double>> ProbePoints; // start points for the PlaceProbes function
		std::vector<int> ProbeComponents;			 // component to measure 0 for x 1 for y
		std::vector<Point<DIM, double>> Offsets;
		std::vector<int> maxIters;

		if (GlobalParameters.SCENARIO == CAVITY)
		{
			ProbePoints.push_back(HalfChannel);
			ProbePoints.push_back(HalfHeight);

			ProbeComponents.push_back(0); // measure x velocity
			ProbeComponents.push_back(1); // measure y velocity

			Offsets.push_back(VerticalOffset);
			Offsets.push_back(HorizontalOffset);

			maxIters.push_back(kendHeight);
			maxIters.push_back(kendWidth);
		}
		else
		{
			ProbePoints.push_back(HalfChannel);
			ProbePoints.push_back(EndChannel);

			ProbeComponents.push_back(0); // measure x velocity
			ProbeComponents.push_back(0); // measure x velocity

			Offsets.push_back(VerticalOffset);
			Offsets.push_back(VerticalOffset);

			maxIters.push_back(kendHeight);
			maxIters.push_back(kendHeight);
		}

		for (int k = 0; k < ProbePoints.size(); k++)
		{
			// create probe object
			Ghost<DIM, double> gp(0);
			size_t bc_p[DIM] = {NON_PERIODIC, NON_PERIODIC};
			probe_particles vp_loc(0, domain, bc_p, gp, DEC_GRAN(512));
			if (ProbeComponents[k] == 0)
			{
				openfpm::vector<std::string> names_p = {"vx"};
				vp_loc.setPropNames(names_p);
			}
			else if (ProbeComponents[k] == 1)
			{
				openfpm::vector<std::string> names_p = {"vy"};
				vp_loc.setPropNames(names_p);
			}

			if (v_cl.getProcessUnitID() == 0)
			{
				PlaceProbes(vp_loc, k0, maxIters[k], ProbePoints[k], Offsets[k]);
			}
			std::pair<probe_particles, int> tmp = std::make_pair(vp_loc, ProbeComponents[k]);
			vp_vec.push_back(tmp);
			GlobalParameters.probe_filenames.push_back("probes_" + std::to_string(k) + "_" + GlobalParameters.filename);
		}
	}

	// Add the obstacle/walls as marker particles only on processor 0
	if (GlobalParameters.BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
	{
		// Add obstacle
		obstacle_ptr->AddObstacle(vd);

		// Add walls
		if (GlobalParameters.bc[0] == PERIODIC && GlobalParameters.bc[1] == NON_PERIODIC) // Channel like scenario
		{
			double dx_wall = dp / refine_factor;
			int Nwall = ceil(GlobalParameters.length[0] / dx_wall);
			dx_wall = GlobalParameters.length[0] / Nwall;
			Point<DIM, double> X_Offset = {dx_wall, 0.0};

			Point<DIM, double> LL_corner = {0.0, 0.0};
			Point<DIM, double> UL_corner = {0.0, GlobalParameters.length[1]};
			// Top And Bottom Walls
			AddFlatWallNewBC(vd, 0, Nwall, LL_corner, X_Offset, dx_wall, {0.0, 0.0}, GlobalParameters.vw_bottom, 0.0);
			AddFlatWallNewBC(vd, 0, Nwall, UL_corner, X_Offset, dx_wall, {0.0, 0.0}, GlobalParameters.vw_top, 0.0);
		}
		else if (GlobalParameters.bc[0] == NON_PERIODIC && GlobalParameters.bc[1] == NON_PERIODIC) // Box like scenario
		{
			double dx_wall_x = dp / refine_factor;
			int Nwall_x = ceil(GlobalParameters.length[0] / dx_wall_x);
			dx_wall_x = GlobalParameters.length[0] / Nwall_x;
			Point<DIM, double> X_Offset = {dx_wall_x, 0.0};

			double dx_wall_y = dp / refine_factor;
			int Nwall_y = ceil(GlobalParameters.length[1] / dx_wall_y);
			dx_wall_y = GlobalParameters.length[1] / Nwall_y;
			Point<DIM, double> Y_Offset = {0.0, dx_wall_y};

			Point<DIM, double> LL_corner = {0.0, 0.0};
			Point<DIM, double> LR_corner = {GlobalParameters.length[0], 0.0};

			Point<DIM, double> UL_corner = {0.0, GlobalParameters.length[1]};

			// Top And Bottom Walls
			AddFlatWallNewBC(vd, 0, Nwall_x + 1, LL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, GlobalParameters.vw_bottom, 0.0);
			AddFlatWallNewBC(vd, 0, Nwall_x + 1, UL_corner, X_Offset, dx_wall_x, {0.0, 0.0}, GlobalParameters.vw_top, 0.0);

			// Left And Right Walls
			AddFlatWallNewBC(vd, 1, Nwall_y, LL_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, 0.0);
			AddFlatWallNewBC(vd, 1, Nwall_y, LR_corner, Y_Offset, dx_wall_y, {0.0, 0.0}, {0.0, 0.0}, 0.0);
		}
	}

	// return an iterator to the fluid particles to add to vd
	auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

	// for each particle inside the fluid box ...
	while (fluid_it.isNext())
	{

		Point<DIM, double> iterator_position = fluid_it.get();

		if ((*obstacle_ptr).isInside(iterator_position)) // if inside the obstacle region
		{
			if (GlobalParameters.BC_TYPE == NO_SLIP) // add particle but set it as boundary
			{
				// ... add a particle ...
				vd.add();
				vd.template getLastProp<type>() = BOUNDARY;
				vd.template getLastProp<vd_omega>() = (*obstacle_ptr).AngularVelocity_;
				for (int xyz = 0; xyz < DIM; xyz++)
				{
					vd.template getLastProp<velocity>()[xyz] = ((*obstacle_ptr).LinearVelocity_).get(xyz);
					vd.template getLastProp<force_transport>()[xyz] = ((*obstacle_ptr).Centre_).get(xyz);
				}
			}
			else if (GlobalParameters.BC_TYPE == NEW_NO_SLIP) // not add particle because already added
			{
				++fluid_it;
				continue;
			}
		}
		else // if not inside obstacle at all just add fluid particles
		{
			// ... add a particle ...
			vd.add();
			vd.template getLastProp<type>() = FLUID;
			vd.template getLastProp<vd_omega>() = 0.0;
			for (int xyz = 0; xyz < DIM; xyz++)
			{
				vd.template getLastProp<velocity>()[xyz] = 0.0;
				vd.template getLastProp<force_transport>()[xyz] = 0.0;
			}
		}
		// Set properties
		vd.template getLastProp<rho>() = GlobalParameters.rho_zero;
		vd.template getLastProp<pressure>() = 0.0;
		vd.template getLastProp<drho>() = 0.0;

		for (int xyz = 0; xyz < DIM; xyz++)
		{
			vd.getLastPos()[xyz] = iterator_position.get(xyz);
			vd.template getLastProp<v_transport>()[xyz] = 0.0;
			vd.template getLastProp<normal_vector>()[xyz] = 0.0;
		}

		vd.template getLastProp<curvature_boundary>() = 0.0;
		vd.template getLastProp<arc_length>() = dp;

		// next fluid particle
		++fluid_it;
	}

	// Now place solid walls using iterators (only for OLD BC)

	if (GlobalParameters.BC_TYPE == NO_SLIP)
	{

		openfpm::vector<Box<DIM, double>> holes;
		holes.add(fluid_box);
		Box<DIM, double> hole_box = holes.get(0);
		auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);

		if (GlobalParameters.bc[0] != PERIODIC || GlobalParameters.bc[1] != PERIODIC) // no walls in all periodic scenario
		{
			while (bound_box.isNext())
			{
				Point<DIM, double> position = bound_box.get();

				// periodic bc, with no boundary particles in y direction has a bug, it puts 3 extra particles outside in the y direction
				// When running on multiple cores, with this we check if particle is outside the recipient box
				// Another bug places boundary particles in the correct plane, but inside the fluid box;
				// if (bc[0] == PERIODIC && position.get(0) > dp / 2.0 && position.get(0) < length[0] - dp / 2.0)
				// {
				// 	++bound_box;
				// 	continue;
				// }

				if (!recipient.isInside((position)))
				{
					++bound_box;
					continue;
				}
				if (hole_box.isInside(position))
				{
					++bound_box;
					continue;
				}

				// if (GlobalParameters.BC_TYPE == NEW_NO_SLIP && (GlobalParameters.bc[0] == NON_PERIODIC && GlobalParameters.bc[1] == NON_PERIODIC))
				// {
				// 	// Check if x and z coordinates are multiples of dp, keep multiples, discard the rest
				// 	double remx = fmod(position.get(0), dp);
				// 	double remz = fmod(position.get(1), dp);
				// 	double tol = 0.5 * dp * 10e-2;

				// 	if (remx > tol && remx < dp - tol)
				// 	{
				// 		++bound_box;
				// 		continue;
				// 	}
				// 	if (remz > tol && remz < dp - tol)
				// 	{
				// 		++bound_box;
				// 		continue;
				// 	}
				// }
				vd.add();

				vd.template getLastProp<type>() = BOUNDARY;
				vd.template getLastProp<rho>() = 0.0;
				vd.template getLastProp<pressure>() = 0.0;
				vd.template getLastProp<drho>() = 0.0;

				for (int xyz = 0; xyz < DIM; xyz++)
				{
					vd.getLastPos()[xyz] = bound_box.get().get(xyz);
					vd.template getLastProp<force>()[xyz] = 0.0;
					vd.template getLastProp<force_transport>()[xyz] = 0.0;
					vd.template getLastProp<v_transport>()[xyz] = 0.0;
					vd.template getLastProp<normal_vector>()[xyz] = 0.0;
					if (position.get(1) < dp / 4.0) // bottom wall
					{
						vd.template getLastProp<velocity>()[xyz] = GlobalParameters.vw_bottom.get(xyz);
					}
					else if (position.get(1) > GlobalParameters.length[1] - dp / 4.0) // top wall
					{
						vd.template getLastProp<velocity>()[xyz] = GlobalParameters.vw_top.get(xyz);
					}
				}

				vd.template getLastProp<curvature_boundary>() = 0.0;
				vd.template getLastProp<arc_length>() = dp;

				++bound_box;
			}
		}
	}

	openfpm::vector<std::string> names({"type",
										"rho",
										"pressure",
										"drho",
										"force",
										"velocity",
										"force_transport",
										"v_transport",
										"normal",
										"curvature",
										"arc_length",
										"volume",
										"vd_omega"});
	vd.setPropNames(names);
}
void CreateParticleGeometryTaylorCouette(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Vcluster<> &v_cl, Obstacle *&obstacle_ptr)
{
	// Size of the virtual cartesian grid that defines where to place the particles
	size_t sz[DIM];

	// In the case of the new bc we need particles at the wall, for this we need sz_aux
	// We want to put one virtual grid point between each pair of the old ones,
	// so that the new spacing is dp/2, and we can put a fluid particle exactly at the wall
	size_t sz_aux[DIM];

	// We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
	double offset_domain_left[DIM] = {0.0};
	double offset_domain_right[DIM] = {0.0};
	double offset_recipient[DIM] = {0.0};
	double offset_periodic_fluid[DIM] = {0.0};
	double offset_periodic_recipient[DIM] = {0.0};

	// non periodic situation grid of 5 fluid particles and 3 boundary particles
	// We need a virtual grid of 5 + 2*(3+1) particles,
	// therefore the domain is discretized with 13 grid points,
	// when we use DrawParticles::DrawBox we will draw only the particles at the grid positons strictly inside the box,
	// the () repesent the recipient box, and the || represent the fluid box, we can see how this distribution places exactly 5 fluid particles inside and 3 boundary particles
	//           D-(-o--o--o-|-x--x--x--x--x--|-o-o-o-)-D
	// D: domain, o: boundary, x: fluid, --: dp distance
	// in a periodic situation we have the following
	// .....--x--x--D-|-x--x--x--x--x--|-D--x--x--......
	// therefore we need a grid of 5 + 2 particles, and the domain is discretized with 7 grid points

	double Rin = GlobalParameters.Rin;
	double Rout = GlobalParameters.Rout;
	double Win = GlobalParameters.Win;
	double Wout = GlobalParameters.Wout;
	double Vin = Win * Rin;
	double Vout = Wout * Rout;
	double a_tc = -((Rout * Rout * Rin * Rin) / (Rout * Rout - Rin * Rin)) * (Wout - Win);
	double b_tc = (Wout * Rout * Rout - Win * Rin * Rin) / (Rout * Rout - Rin * Rin);

	size_t Nbound = (GlobalParameters.BC_TYPE == NEW_NO_SLIP) ? 1 : 3;

	for (int dim = 0; dim < DIM; dim++)
	{
		GlobalParameters.length[dim] = dp * GlobalParameters.Nfluid[dim];
		sz[dim] = GlobalParameters.Nfluid[dim] + 2 * (Nbound + 1);
		offset_domain_left[dim] = (0.5 + Nbound) * dp;
		offset_domain_right[dim] = (0.5 + Nbound) * dp;
	}

	// Define the boxes
	Box<DIM, double> domain({-GlobalParameters.length[0] / 2.0 - offset_domain_left[0],
							 -GlobalParameters.length[1] / 2.0 - offset_domain_left[1]},
							{GlobalParameters.length[0] / 2.0 + offset_domain_right[0],
							 GlobalParameters.length[1] / 2.0 + offset_domain_right[1]});
	Box<DIM, double> fluid_box({-GlobalParameters.length[0] / 2.0,
								-GlobalParameters.length[1] / 2.0},
							   {GlobalParameters.length[0] / 2.0,
								GlobalParameters.length[1] / 2.0});

	// extended boundary around the domain, and the processor domain
	Ghost<DIM, double> g(GlobalParameters.r_threshold);

	// create particle object
	particles vd_loc(0, domain, GlobalParameters.bc, g, DEC_GRAN(512));
	vd = vd_loc;

	// correct the number of particles in case of periodicity, we substracted 1 before to accomodate the periodic boundary
	for (int dim = 0; dim < DIM; dim++)
	{
		if (GlobalParameters.bc[dim] == PERIODIC)
		{
			GlobalParameters.Nfluid[dim] += 1;
			GlobalParameters.length[dim] += dp;
		}
	}

	// Write constants on file
	double rf = GlobalParameters.rf;

	WriteParameters(v_cl);

	// Set cylindrical object parameters
	Point<DIM, double> CylinderCentre = {0.0, 0.0};

	obstacle_ptr = new EmptyObstacle();

	const Point<DIM, double> vel = {0.0, 0.0};

	CylinderObstacle *obstacle_ptr_out = new CylinderObstacle(Rout, CylinderCentre, vel, Wout, rf);
	CylinderObstacle *obstacle_ptr_in = new CylinderObstacle(Rin, CylinderCentre, vel, Win, rf);

	CylinderObstacle *obstacle_ptr_out_aux = new CylinderObstacle(Rout + 3.0 * dp, CylinderCentre, vel, Wout, rf);
	CylinderObstacle *obstacle_ptr_in_aux = new CylinderObstacle(Rin - 3.0 * dp, CylinderCentre, vel, Win, rf);

	// Add the obstacle as marker particles only on processor 0
	if (GlobalParameters.BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
	{
		obstacle_ptr_out->AddObstacle(vd);
		obstacle_ptr_in->AddObstacle(vd);
	}

	Box<DIM, double> fluid_box_aux({-Rout - 3.5 * dp,
									-Rout - 3.5 * dp},
								   {Rout + 3.5 * dp,
									Rout + 3.5 * dp});
	auto out_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box_aux);
	// for each particle inside the fluid box ...
	if (GlobalParameters.BC_TYPE == NO_SLIP)
	{
		while (out_it.isNext())
		{

			Point<DIM, double> iterator_position = out_it.get();
			if (!(*obstacle_ptr_out).isOutside(iterator_position) && (*obstacle_ptr_out_aux).isInside(iterator_position)) // if outside the outer cylinder and inside outer cylinder aux
			{
				if (GlobalParameters.BC_TYPE == NO_SLIP)
				{
					vd.add();
					// Set properties
					vd.template getLastProp<type>() = BOUNDARY;
					vd.template getLastProp<vd_omega>() = (*obstacle_ptr_out).AngularVelocity_;
					for (int xyz = 0; xyz < DIM; xyz++)
					{
						vd.template getLastProp<velocity>()[xyz] = ((*obstacle_ptr_out).LinearVelocity_).get(xyz);
						vd.template getLastProp<force_transport>()[xyz] = ((*obstacle_ptr_out).Centre_).get(xyz);
					}
					vd.template getLastProp<rho>() = GlobalParameters.rho_zero;
					vd.template getLastProp<pressure>() = 0.0;
					vd.template getLastProp<drho>() = 0.0;

					for (int xyz = 0; xyz < DIM; xyz++)
					{
						vd.getLastPos()[xyz] = iterator_position.get(xyz);
						vd.template getLastProp<v_transport>()[xyz] = 0.0;
						vd.template getLastProp<normal_vector>()[xyz] = 0.0;
					}

					vd.template getLastProp<curvature_boundary>() = 0.0;
					vd.template getLastProp<arc_length>() = dp;

					// next fluid particle
					++out_it;
					continue;
				}
				else
				{
					++out_it;
					continue;
				}
			}
			else // skip fluid particle
			{
				++out_it;
				continue;
			}
		}
	}

	// return an iterator to the fluid particles to add to vd
	auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

	// for each particle inside the fluid box ...
	while (fluid_it.isNext())
	{

		Point<DIM, double> iterator_position = fluid_it.get();
		if ((*obstacle_ptr_out).isOutside(iterator_position)) // if inside the outer cylinder
		{
			if ((*obstacle_ptr_in).isInside(iterator_position)) // if inside the inner cylinder region
			{
				if (!(*obstacle_ptr_in_aux).isOutside(iterator_position))
				{
					if (GlobalParameters.BC_TYPE == NO_SLIP) // add particle but set it as boundary
					{
						// ... add a particle ...
						vd.add();
						vd.template getLastProp<type>() = BOUNDARY;
						vd.template getLastProp<vd_omega>() = (*obstacle_ptr_in).AngularVelocity_;
						for (int xyz = 0; xyz < DIM; xyz++)
						{
							vd.template getLastProp<velocity>()[xyz] = ((*obstacle_ptr_in).LinearVelocity_).get(xyz);
							vd.template getLastProp<force_transport>()[xyz] = ((*obstacle_ptr_in).Centre_).get(xyz);
						}
					}
					else if (GlobalParameters.BC_TYPE == NEW_NO_SLIP) // not add particle because already added
					{
						++fluid_it;
						continue;
					}
				}
				else
				{
					++fluid_it;
					continue;
				}
			}
			else // if no cylinder at all just add fluid particles
			{
				// ... add a particle ...
				vd.add();
				vd.template getLastProp<type>() = FLUID;
				vd.template getLastProp<vd_omega>() = 0.0;
				double r = iterator_position.get(0) * iterator_position.get(0) + iterator_position.get(1) * iterator_position.get(1);
				r = sqrt(r);
				double theta = atan2(iterator_position.get(1), iterator_position.get(0));
				double uth = a_tc / r + b_tc * r;
				double ur = 0.0;

				double ux = uth * (-iterator_position.get(1) / r);
				double uy = uth * (iterator_position.get(0) / r);

				vd.template getLastProp<velocity>()[0] = ux;
				vd.template getLastProp<velocity>()[1] = uy;

				for (int xyz = 0; xyz < DIM; xyz++)
				{
					vd.template getLastProp<force_transport>()[xyz] = 0.0;
				}
			}
		}
		else // skip fluid particle
		{
			++fluid_it;
			continue;
		}
		// Set properties
		vd.template getLastProp<rho>() = GlobalParameters.rho_zero;
		vd.template getLastProp<pressure>() = 0.0;
		vd.template getLastProp<drho>() = 0.0;

		for (int xyz = 0; xyz < DIM; xyz++)
		{
			vd.getLastPos()[xyz] = iterator_position.get(xyz);
			vd.template getLastProp<v_transport>()[xyz] = 0.0;
			vd.template getLastProp<normal_vector>()[xyz] = 0.0;
		}

		vd.template getLastProp<curvature_boundary>() = 0.0;
		vd.template getLastProp<arc_length>() = dp;

		// next fluid particle
		++fluid_it;
	}

	openfpm::vector<std::string> names({"type",
										"rho",
										"pressure",
										"drho",
										"force",
										"velocity",
										"force_transport",
										"v_transport",
										"normal",
										"curvature",
										"arc_length",
										"volume",
										"vd_omega"});
	vd.setPropNames(names);
}
void CreateParticleGeometryStep(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Vcluster<> &v_cl)
{
	// Size of the virtual cartesian grid that defines where to place the particles
	size_t sz[DIM];

	double length_small[DIM];
	double length_big[DIM];
	// In the case of the new bc we need particles at the wall, for this we need sz_aux
	// We want to put one virtual grid point between each pair of the old ones,
	// so that the new spacing is dp/2, and we can put a fluid particle exactly at the wall
	size_t sz_aux[DIM];

	// Boundary conditions
	size_t bc[DIM];

	// Number of boundary particles in each direction
	size_t Nboundary_big[DIM];

	// Number of fluid particles in each direction
	size_t Nfluid_big[DIM];
	size_t Nfluid_small[DIM];

	// We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
	double offset_domain[DIM] = {0.0};
	double offset_recipient[DIM] = {0.0};
	double offset_periodic_fluid[DIM] = {0.0};
	double offset_periodic_recipient[DIM] = {0.0};

	// non periodic situation grid of 5 fluid particles and 3 boundary particles
	// We need a virtual grid of 5 + 2*(3+1) particles,
	// therefore the domain is discretized with 13 grid points,
	// when we use DrawParticles::DrawBox we will draw only the particles at the grid positons strictly inside the box,
	// the () repesent the recipient box, and the || represent the fluid box, we can see how this distribution places exactly 5 fluid particles inside and 3 boundary particles
	//           D-(-o--o--o-|-x--x--x--x--x--|-o-o-o-)-D
	// D: domain, o: boundary, x: fluid, --: dp distance
	// in a periodic situation we have the following
	// .....--x--x--D-|-x--x--x--x--x--|-D--x--x--......
	// therefore we need a grid of 5 + 2 particles, and the domain is discretized with 7 grid points

	double StepHeight = 4.9;
	GlobalParameters.LengthScale = StepHeight;

	Nfluid_big[0] = 275;
	Nfluid_big[1] = 31;

	Nfluid_small[0] = 100;
	Nfluid_small[1] = 16;

	size_t Nbound = (GlobalParameters.BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
	Nboundary_big[0] = 0;
	Nboundary_big[1] = Nbound;
	double Nboundary_small_up = Nbound;
	double Nboundary_small_down = Nbound + Nfluid_big[1] - Nfluid_small[1];

	bc[0] = PERIODIC;
	bc[1] = NON_PERIODIC;
	dp = GlobalParameters.LengthScale / ((double)Nfluid_big[1] - Nfluid_small[1]);
	GlobalParameters.umax = 1.4 * 1e-1;

	GlobalParameters.H = GlobalParameters.Hconst * dp;
	// r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
	GlobalParameters.r_threshold = 3.0 * GlobalParameters.H;
	GlobalParameters.Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / GlobalParameters.H / GlobalParameters.H / GlobalParameters.H : 7.0 / 478.0 / M_PI / GlobalParameters.H / GlobalParameters.H;
	GlobalParameters.MassFluid = GlobalParameters.rho_zero * (DIM == 3 ? dp * dp * dp : dp * dp);
	GlobalParameters.MassBound = GlobalParameters.rho_zero * (DIM == 3 ? dp * dp * dp : dp * dp);
	GlobalParameters.cbar = GlobalParameters.coeff_sound * GlobalParameters.umax;
	GlobalParameters.B = GlobalParameters.rho_zero * GlobalParameters.cbar * GlobalParameters.cbar / GlobalParameters.gamma_;
	GlobalParameters.Pbackground = GlobalParameters.Bfactor * GlobalParameters.B;
	GlobalParameters.eta = GlobalParameters.nu * GlobalParameters.rho_zero;
	GlobalParameters.Re = GlobalParameters.umax * 2.0 * 5.2 / GlobalParameters.nu;

	GlobalParameters.gravity = getVectorNorm(GlobalParameters.gravity_vector);

	for (int dim = 0; dim < DIM; dim++)
	{
		if (bc[dim] == NON_PERIODIC) // non periodic, fluid covered by boundary
		{
			GlobalParameters.length[dim] = dp * (Nfluid_big[dim]);
			length_small[dim] = dp * (Nfluid_small[dim]);
			length_big[dim] = dp * (Nfluid_big[dim]);
			sz[dim] = Nfluid_big[dim] + 2 * (Nboundary_big[dim] + 1);
			offset_domain[dim] = (0.5 + Nboundary_big[dim]) * dp;

			if (Nboundary_big[dim] != 0)
				sz_aux[dim] = 2 * Nfluid_big[dim] - 1 + 2 * (2 * Nboundary_big[dim] + 1 + 1);
			else // for a direction with no boundary particles we dont need to add anything
				sz_aux[dim] = sz[dim];

			if (GlobalParameters.BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
				offset_recipient[dim] = 0.25 * Nboundary_big[dim] * dp;
			else if (GlobalParameters.BC_TYPE == NO_SLIP)
				offset_recipient[dim] = Nboundary_big[dim] * dp;
		}
		else // periodic, open ended
		{
			Nfluid_big[dim] -= 1;
			GlobalParameters.length[dim] = dp * (Nfluid_big[dim] + Nfluid_small[dim]);
			length_small[dim] = dp * (Nfluid_small[dim]);
			length_big[dim] = dp * (Nfluid_big[dim]);

			sz[dim] = (Nfluid_big[dim] + Nfluid_small[dim]) + 2;
			sz_aux[dim] = sz[dim];

			offset_domain[dim] = 0.5 * dp;
			offset_periodic_fluid[dim] = 0.75 * dp;
			offset_periodic_recipient[dim] = 0.85 * dp;
		}
	}

	// Define the boxes
	Box<DIM, double> domain({-offset_domain[0],
							 -offset_domain[1]},
							{GlobalParameters.length[0] + offset_domain[0],
							 GlobalParameters.length[1] + offset_domain[1]});

	Box<DIM, double> fluid_box_small({0.0,
									  (Nfluid_big[1] - Nfluid_small[1]) * dp},
									 {length_small[0],
									  (Nfluid_big[1] - Nfluid_small[1]) * dp + length_small[1]});

	Box<DIM, double> fluid_box_big({length_small[0],
									0.0},
								   {length_small[0] + length_big[0] + offset_periodic_fluid[0],
									length_big[1]});

	Box<DIM, double> recipient({-offset_recipient[0],
								-offset_recipient[1]},
							   {GlobalParameters.length[0] + offset_recipient[0] + offset_periodic_recipient[0],
								GlobalParameters.length[1] + offset_recipient[1] + offset_periodic_recipient[1]});

	// Will only be used in the new bc
	Box<DIM, double> recipient_hole_small({offset_recipient[0],
										   (Nfluid_big[1] - Nfluid_small[1]) * dp + offset_recipient[1]},
										  {length_small[0] - offset_recipient[0],
										   (Nfluid_big[1] - Nfluid_small[1]) * dp + length_small[1] - offset_recipient[1]});

	Box<DIM, double> recipient_hole_big({length_small[0] + offset_recipient[0],
										 offset_recipient[1]},
										{length_small[0] + length_big[0] - offset_recipient[0] + offset_periodic_fluid[0],
										 length_big[1] - offset_recipient[1]});

	Box<DIM, double> CornerHole{{3 * dp, -3 * dp}, {(3 + Nfluid_small[0] - 6) * dp, (Nfluid_big[1] - Nfluid_small[1] - 3) * dp}};
	Box<DIM, double> CornerHole_New{{dp, -1 * dp}, {(1 + Nfluid_small[0] - 2) * dp, (Nfluid_big[1] - Nfluid_small[1]) * dp - 0.5 * dp}};

	// extended boundary around the domain, and the processor domain
	Ghost<DIM, double> g(GlobalParameters.r_threshold);

	// create particle object
	particles vd_loc(0, domain, bc, g, DEC_GRAN(512));
	vd = vd_loc;

	// correct the number of particles in case of periodicity, we substracted 1 before to accomodate the periodic boundary
	for (int dim = 0; dim < DIM; dim++)
	{
		if (bc[dim] == PERIODIC)
		{
			Nfluid_big[dim] += 1;
			GlobalParameters.length[dim] += dp;
		}
	}

	// Write constants on file
	WriteParameters(v_cl);

	// we want to place probes  in a vertical line at this locations
	Point<DIM, double> P1 = {0.75 * Nfluid_small[0] * dp, 0.0};
	Point<DIM, double> P2 = {1.2 * Nfluid_small[0] * dp, 0.0};
	Point<DIM, double> P3 = {1.4 * Nfluid_small[0] * dp, 0.0};
	Point<DIM, double> P4 = {Nfluid_small[0] * dp + Nfluid_big[0] * dp * 0.7, 0.0};

	std::vector<Point<DIM, double>> ProbePoints = {P1, P2, P3, P4};

	Point<DIM, double> VerticalOffset = {0.0, dp};
	int k0 = 0;
	int kendHeight = Nfluid_big[1];

	// place probes
	if (GlobalParameters.PROBES_ENABLED)
	{
		for (int k = 0; k < 4; k++)
		{
			// create probe object
			Ghost<DIM, double> gp(0);
			size_t bc_p[DIM] = {NON_PERIODIC, NON_PERIODIC};
			probe_particles vp_loc(0, domain, bc_p, gp, DEC_GRAN(512));
			openfpm::vector<std::string> names_p({"vx"});
			vp_loc.setPropNames(names_p);

			PlaceProbes(vp_loc, k0, kendHeight, ProbePoints[k], VerticalOffset);
			std::pair<probe_particles, int> tmp = std::make_pair(vp_loc, 0);
			vp_vec.push_back(tmp);
			GlobalParameters.probe_filenames.push_back("probes_" + std::to_string(k) + "_" + GlobalParameters.filename);
		}
	}

	// return an iterator to the fluid particles to add to vd
	auto fluid_it1 = DrawParticles::DrawBox(vd, sz, domain, fluid_box_big);
	auto fluid_it2 = DrawParticles::DrawBox(vd, sz, domain, fluid_box_small);

	// for each particle inside the fluid box ...
	while (fluid_it1.isNext())
	{

		Point<DIM, double> iterator_position = fluid_it1.get();

		// ... add a particle ...
		vd.add();
		vd.template getLastProp<type>() = FLUID;

		// Set properties
		vd.template getLastProp<rho>() = GlobalParameters.rho_zero;
		vd.template getLastProp<pressure>() = 0.0;
		vd.template getLastProp<drho>() = 0.0;

		for (int xyz = 0; xyz < DIM; xyz++)
		{
			vd.getLastPos()[xyz] = iterator_position.get(xyz);
			vd.template getLastProp<force_transport>()[xyz] = 0.0;
			vd.template getLastProp<v_transport>()[xyz] = 0.0;
			vd.template getLastProp<normal_vector>()[xyz] = 0.0;
			vd.template getLastProp<velocity>()[xyz] = 0.0;
		}

		vd.template getLastProp<curvature_boundary>() = 0.0;
		vd.template getLastProp<arc_length>() = dp;

		// next fluid particle
		++fluid_it1;
	}
	// for each particle inside the fluid box ...
	while (fluid_it2.isNext())
	{

		Point<DIM, double> iterator_position = fluid_it2.get();

		// ... add a particle ...
		vd.add();
		vd.template getLastProp<type>() = FLUID;

		// Set properties
		vd.template getLastProp<rho>() = GlobalParameters.rho_zero;
		vd.template getLastProp<pressure>() = 0.0;
		vd.template getLastProp<drho>() = 0.0;

		for (int xyz = 0; xyz < DIM; xyz++)
		{
			vd.getLastPos()[xyz] = iterator_position.get(xyz);
			vd.template getLastProp<force_transport>()[xyz] = 0.0;
			vd.template getLastProp<v_transport>()[xyz] = 0.0;
			vd.template getLastProp<normal_vector>()[xyz] = 0.0;
			vd.template getLastProp<velocity>()[xyz] = 0.0;
		}

		vd.template getLastProp<curvature_boundary>() = 0.0;
		vd.template getLastProp<arc_length>() = dp;

		// next fluid particle
		++fluid_it2;
	}

	// Now place solid walls
	openfpm::vector<Box<DIM, double>> holes;

	if (GlobalParameters.BC_TYPE == NEW_NO_SLIP)
	{
		holes.add(recipient_hole_small);
		holes.add(recipient_hole_big);
		holes.add(CornerHole_New);
		sz[0] = sz_aux[0];
		sz[1] = sz_aux[1];
	}
	else if (GlobalParameters.BC_TYPE == NO_SLIP)
	{
		holes.add(fluid_box_big);
		holes.add(fluid_box_small);
		holes.add(CornerHole);
	}
	Box<DIM, double> hole_get0 = holes.get(0);
	Box<DIM, double> hole_get1 = holes.get(1);
	auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);

	if (bc[0] != PERIODIC || bc[1] != PERIODIC) // no walls in all periodic scenario
	{
		while (bound_box.isNext())
		{
			Point<DIM, double> position = bound_box.get();

			// periodic bc, with no boundary particles in y direction has a bug, it puts 3 extra particles outside in the y direction
			// When running on multiple cores, with this we check if particle is outside the recipient box
			// Another bug places boundary particles in the correct plane, but inside the fluid box;
			// the first bug seems to be fixed
			// if (!recipient.isInside((position)))
			// {
			// 	++bound_box;
			// 	continue;
			// }
			// if (bc[0] == PERIODIC && position.get(0) > dp / 2.0 && position.get(0) < length[0] - dp / 2.0)
			// {
			// 	++bound_box;
			// 	continue;
			// }

			if (hole_get0.isInside(position) || hole_get1.isInside(position))
			{
				++bound_box;
				continue;
			}

			if (GlobalParameters.BC_TYPE == NEW_NO_SLIP)
			{
				// Check if x and y coordinates are multiples of dp, keep multiples, discard the rest
				double remx = fmod(position.get(0), dp);
				double remy = fmod(position.get(1), dp);
				double tol = 0.5 * dp * 10e-2;

				// if (remx > tol && remx < dp - tol)
				// {
				// 	++bound_box;
				// 	continue;
				// }
				if (remy > tol && remy < dp - tol)
				{
					++bound_box;
					continue;
				}
			}
			vd.add();

			vd.template getLastProp<type>() = BOUNDARY;
			vd.template getLastProp<rho>() = GlobalParameters.rho_zero;
			vd.template getLastProp<pressure>() = 0.0;
			vd.template getLastProp<drho>() = 0.0;

			for (int xyz = 0; xyz < DIM; xyz++)
			{
				vd.getLastPos()[xyz] = bound_box.get().get(xyz);
				vd.template getLastProp<force>()[xyz] = 0.0;
				vd.template getLastProp<force_transport>()[xyz] = 0.0;
				vd.template getLastProp<v_transport>()[xyz] = 0.0;
				vd.template getLastProp<normal_vector>()[xyz] = 0.0;
				if (position.get(1) < dp / 4.0) // bottom wall
				{
					vd.template getLastProp<velocity>()[xyz] = GlobalParameters.vw_bottom.get(xyz);
				}
				else if (position.get(1) > GlobalParameters.length[1] - dp / 4.0) // top wall
				{
					vd.template getLastProp<velocity>()[xyz] = GlobalParameters.vw_top.get(xyz);
				}
			}

			vd.template getLastProp<curvature_boundary>() = 0.0;
			vd.template getLastProp<arc_length>() = dp;

			++bound_box;
		}
	}

	openfpm::vector<std::string> names({"type",
										"rho",
										"pressure",
										"drho",
										"force",
										"velocity",
										"force_transport",
										"v_transport",
										"normal",
										"curvature",
										"arc_length",
										"volume",
										"vd_omega"});
	vd.setPropNames(names);
}
int main(int argc, char *argv[])
{

	// initialize the library
	openfpm_init(&argc, &argv);

	std::string tmp = argv[1];

	ParseXMLFile(tmp);
	InitializeConstants();

	// create a Vcluster object ( like MPI communicator )
	Vcluster<> &v_cl = create_vcluster();

	// Create a particle vector
	particles vd;
	std::vector<std::pair<probe_particles, int>> vp;

	Obstacle *obstacle_ptr = nullptr;
	if (GlobalParameters.SCENARIO == STEP)
	{
		CreateParticleGeometryStep(vd, vp, v_cl);
	}
	else if (GlobalParameters.SCENARIO == TAYLOR_COUETTE)
	{
		CreateParticleGeometryTaylorCouette(vd, vp, v_cl, obstacle_ptr);
	}
	else
	{
		CreateParticleGeometry(vd, vp, v_cl, obstacle_ptr);
	}
	vd.map();
	vd.write_frame(GlobalParameters.filename, 0, GlobalParameters.WRITER);

	for (int k = 0; k < vp.size(); k++)
	{
		(vp[k].first).map();
	}
	// Now that we fill the vector with particles, decompose the domain
	// ModelCustom md;
	// vd.addComputationCosts(md);
	// vd.getDecomposition().decompose();
	// vd.map();

	// for (int k = 0; k < vp.size(); k++)
	// {
	// 	(vp[k].first).getDecomposition().decompose();
	// 	(vp[k].first).map();
	// }
	vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

	auto NN = vd.getCellList(GlobalParameters.r_threshold);
	// vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();
	// vd.updateCellList(NN);

	if (GlobalParameters.BC_TYPE == NO_SLIP) // set up boundary particle velocity for the first iteration
	{
		ExtrapolateVelocity(vd, NN);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();
	}
	else if (GlobalParameters.BC_TYPE == NEW_NO_SLIP) // Set up fluid vector and normal vector of the boundary particles
	{
		calcFluidVec(vd, NN);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

		calcNormalVec(vd, NN);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

		calcCurvature(vd, NN, v_cl);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

		calcVolume(vd);
		vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();
	}

	// Evolve
	size_t write = 1;
	size_t it = 0;
	size_t it_reb = 0;
	double t = 0.0;
	std::ofstream avgvelstream("avgvel.csv");
	bool calc_drag = false;
	double cylinder_force = 0.0;

	while (t <= GlobalParameters.t_end)
	{

		//// Do rebalancing every 200 timesteps
		// it_reb++;
		// if (it_reb == 500)
		// {
		// 	vd.map();

		// 	it_reb = 0;
		// 	ModelCustom md;
		// 	vd.addComputationCosts(md);
		// 	vd.getDecomposition().decompose();
		// 	vd.map();

		// 	if (v_cl.getProcessUnitID() == 0)
		// 		std::cout << "REBALANCED " << std::endl;
		// }

		// Calculate dt for time stepping
		double dt = calc_deltaT(vd, v_cl);

		// in general we dont compute drag coeficient every time step
		calc_drag = false;
		if (write < (t + dt) * GlobalParameters.write_const)
		{
			calc_drag = true;
		}

		// Integrate one time step
		kick_drift_int(vd, NN, dt, v_cl, cylinder_force, calc_drag, t);

		// increment time
		t += dt;
		if (write < t * GlobalParameters.write_const)
		{
			// // sensor calculation require ghost and update cell-list
			if (GlobalParameters.PROBES_ENABLED)
			{
				vd.map();
				vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

				vd.updateCellList(NN);
				for (int k = 0; k < vp.size(); k++)
				{
					probe_particles &probe = vp[k].first;
					// probe.map();
					int &component = vp[k].second;
					sensor_velocity_comp(vd, probe, v_cl, NN, component, obstacle_ptr);
					probe.write_frame(GlobalParameters.probe_filenames[k], write, GlobalParameters.WRITER);
				}
			}

			vd.deleteGhost();
			vd.write_frame(GlobalParameters.filename, write, GlobalParameters.WRITER);
			vd.ghost_get<type, rho, pressure, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length, vd_volume, vd_omega>();

			write++;
			computeAverageVelocity(vd, v_cl, t, avgvelstream, cylinder_force);
			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "TIME: " << t << "  write " << GlobalParameters.cnt << std::endl;
			}
		}
		// else
		// {
		// 	if (v_cl.getProcessUnitID() == 0)
		// 	{
		// 		std::cout << "TIME: " << t << "   " << cnt << std::endl;
		// 	}
		// }
	}
	delete obstacle_ptr;
	openfpm_finalize();
}
