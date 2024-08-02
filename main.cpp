#include "Vector/vector_dist.hpp"
#include "Draw/DrawParticles.hpp"
#include <math.h>
#include <sys/stat.h>
#include <cmath>

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
#define CAVITY 7
#define STEP 8
#define TRIANGLE_SYM 9
#define TAYLOR_COUETTE 10

// TYPE OF KERNERL
// #define CUBIC 0
// #define QUINTIC 1

// TYPE OF DENSITY COMPUTATION
#define DENSITY_SUMMATION 0
#define DENSITY_DIFFERENTIAL 1

// DIMENSIONALITY, changes the normalizations of the kernels
#define DIM 2

const int BC_TYPE = NEW_NO_SLIP;
const int SCENARIO = CYLINDER_ARRAY;
// const int KERNEL = QUINTIC;
const int DENSITY_TYPE = DENSITY_SUMMATION;
const int WRITER = VTK_WRITER; // VTK_WRITER or CSV_WRITER
int PROBES_ENABLED = 0;		   // 0 for disabled, 1 for enabled

//////// DECLARATION OF GLOBAL PARAMETERS /////////////////////////////////////////////////////////////
// Output file name
std::string filename;
std::vector<std::string> probe_filenames;
// Initial particle spacing
double dp;
// Physical size of the fluid domain, it goes from (0,0,0) to (length[0],length[1],length[2])
// First particle will always be placed at (dp/2,dp/2,dp/2) and the last particle will be placed at (length[0]-dp/2,length[1]-dp/2,length[2]-dp/2)
double length[DIM];
// problem specific length scale to compute the Reynolds number
double LengthScale;

// Factor relating H (smoothing length) and dp (particle spacing)
const double Hconst = 1.0;
// Smoothing length
double H;
// Radius of the kernel support
double r_threshold;
// Normalization constant for the kernels
double Kcubic;
double Kquintic;
// Reynolds number
double Re;
// maximum velocity
double umax;
// Reference density
double rho_zero;
// Gamma in eq of state
const double gamma_ = 1.0;
// Constant used for the sound speed, number of times the max velocity
double coeff_sound = 10.0;
// Sound speed
double cbar;
// Eq of state constant
double B;
// background pressure in eq of state
const double xi = 0.0;

const double max_angle = 5.0; // degrees

// Gravity vector and magnitude
Point<DIM, double> gravity_vector = {0.0};
double gravity;
// Wall velocity
Point<DIM, double> vw_top = {0.0};
Point<DIM, double> vw_bottom = {0.0};

// Mass of the fluid and boundary particles M=rho*V, V=dp^3 or V = dp^2
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
// Controls otput file frequency, low means less frequent
double write_const;
// iteration counter
size_t cnt = 0;

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

typedef vector_dist<DIM, double, aggregate<size_t, double, double, double, double[DIM], double[DIM], double[DIM], double[DIM], double[DIM], double, double>> particles;
//                                       |         |     |        |        |        |           |		     |             |		|       |
//                                       |         |     |        |        |        |           |		     |			   |		|       |
//                                    type        rho   pressure delta   force     velocity force_transport   v_transport    normal	curvature   arc_length
//                                                              density

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

// Eq state, compute pressure given density, for all particles in vd
void EqState(particles &vd)
{
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto a = it.get();

		if (vd.getProp<type>(a) == FLUID)
		{
			double rho_a = vd.template getProp<rho>(a);
			double rho_frac = rho_a / rho_zero;

			vd.template getProp<pressure>(a) = B * (std::pow(rho_frac, gamma_) - 1.0) + xi;
		}

		++it;
	}
}

// Inverted equation of state, compute density given pressure, particle wise
double InvEqState_particle(const double p)
{
	return rho_zero * std::pow(((p - xi) / B + 1.0), 1.0 / gamma_);
}

double EqState_particle(const double rho)
{
	return B * (std::pow(rho / rho_zero, gamma_) - 1.0) + xi;
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
inline double getVectorNorm(const Point<DIM, double> &v)
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
	r /= H;
	const double tmp3 = (3.0 - r) * (3.0 - r) * (3.0 - r) * (3.0 - r) * (3.0 - r);
	const double tmp2 = (2.0 - r) * (2.0 - r) * (2.0 - r) * (2.0 - r) * (2.0 - r);
	const double tmp1 = (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
	if (r >= 0.0 && r < 1.0)
		return Kquintic * (tmp3 - 6.0 * tmp2 + 15.0 * tmp1);
	else if (r >= 1.0 && r < 2.0)
		return Kquintic * (tmp3 - 6.0 * tmp2);
	else if (r >= 2.0 && r < 3.0)
		return Kquintic * tmp3;
	else
		return 0.0;
}

Point<DIM, double> DWab(const Point<DIM, double> &dx, const double r)
{
	Point<DIM, double> DW;
	const double q = r / H;

	const double tmp3 = (3.0 - q) * (3.0 - q) * (3.0 - q) * (3.0 - q);
	const double tmp2 = (2.0 - q) * (2.0 - q) * (2.0 - q) * (2.0 - q);
	const double tmp1 = (1.0 - q) * (1.0 - q) * (1.0 - q) * (1.0 - q);

	const double c1 = Kquintic * (-5.0 / H);
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
	return eta * (dv * dotProduct(dr, dW)) / (r * r);
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
					if (r2 < r_threshold * r_threshold)
					{
						normalizeVector(dr);
						r_fluid_sum += dr;
					}
				}

				++Np;
			}
			normalizeVector(r_fluid_sum);
			// store in force_transport because we are not using it for boundary particles

			for (int xyz = 0; xyz < DIM; ++xyz)
			{
				vd.template getProp<force_transport>(a)[xyz] = r_fluid_sum.get(xyz);
			}
		}

		++part;
	}
}

template <typename CellList>
void calcNormalVec(particles &vd, CellList &NN)
{
	// This function computes the vector computes the normal vector for a boundary particle based on the other boundary particles inside its kernel.

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
			Point<DIM, double> Rfluid = vd.getProp<force_transport>(a);

			// initialize sum
			Point<DIM, double> n_sum = (DIM == 2) ? Point<DIM, double>{0.0, 0.0} : Point<DIM, double>{0.0, 0.0, 0.0};
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

					// Get the vector pointing at b from a
					const Point<DIM, double> dr = xa - xb;

					// take the norm squared of this vector
					const double r2 = norm2(dr);

					// If the particles interact ...
					if (r2 < r_threshold * r_threshold)
					{
						// get perpendicular vector to dr
						Point<DIM, double> perp = getPerpendicularUnit(dr); // this is normalized

						// get scalar product of perp and Rfluid
						double perp_dot_fluid = dotProduct(perp, Rfluid);

						// we want perp to point towards the fluid
						if (perp_dot_fluid < 0.0)
							perp = -1.0 * perp;

						// evaluate kernel
						double W = Wab(std::sqrt(r2));
						n_sum += perp * W;
					}
				}

				++Np;
			}
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

			// if (vd.getProp<pressure>(a) == 0.0) // flat wall
			// {
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
					if (r2 < r_threshold * r_threshold)
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
			// }
			// else // obstacle
			// {
			// 	// iterate the neighborhood particles
			// 	while (Np.isNext() == true)
			// 	{
			// 		// Key of b particle
			// 		const unsigned long b = Np.get();

			// 		if (vd.getProp<type>(b) == BOUNDARY)
			// 		{
			// 			// Get the position xb of the particle b
			// 			const Point<DIM, double> xb = vd.getPos(b);

			// 			// Get the vector pointing at a from b
			// 			Point<DIM, double> dr = xa - xb;

			// 			// take the norm squared of this vector
			// 			const double r2 = norm2(dr);
			// 			const double r = sqrt(r2);
			// 			// the distance between the two particles is the clock distance
			// 			const unsigned int counter_a = (int)std::round(vd.getProp<pressure>(a));
			// 			const unsigned int counter_b = (int)std::round(vd.getProp<pressure>(b));

			// 			const unsigned int clock_dist = clockDifference(counter_a, counter_b, max_counter_int);
			// 			const double rclk = clock_dist * vd.getProp<arc_length>(a);
			// 			// If the particles interact ...
			// 			if (rclk < r_threshold)
			// 			{
			// 				// OLD CALCULATION
			// 				// Point<DIM, double> normal_b = vd.getProp<normal_vector>(b);

			// 				// double W = Wab(rclk);
			// 				// normalizeVector(dr);
			// 				// dr = rclk * dr;
			// 				// Point<DIM, double> dW = DWab(dr, rclk);

			// 				// K_sum += dotProduct(normal_b - normal_a, dW); // divergence of the normal vector
			// 				// w_sum += W;

			// 				// NEW CALCULATION

			// 				if (a.getKey() != b)
			// 				{
			// 					Point<DIM, double> normal_b = vd.getProp<normal_vector>(b);
			// 					double W = Wab(rclk);
			// 					Point<DIM, double> eab = -1.0 * dr / r;
			// 					double local_k = dotProduct(normal_b - normal_a, eab) / rclk;

			// 					K_sum += local_k * W;
			// 					w_sum += W;
			// 				}
			// 			}
			// 		}

			// 		++Np;
			// 	}
			// }

			K_sum = K_sum / w_sum;
			// store in curvature
			vd.template getProp<curvature_boundary>(a) = K_sum;
		}
		++part;
	}
}

template <typename CellList>
void fix_mass(particles &vd, CellList &NN)
{
	// This function adjusts the mass of particles so that the initial density (determined by summation) is rho_zero

	auto part = vd.getDomainIterator();

	// Update the cell-list
	vd.updateCellList(NN);
	double global_sum = 0.0;
	int count = 0;

	// For each particle ...
	while (part.isNext())
	{
		// Key of the particle a
		vect_dist_key_dx a = part.get();

		// if particle FLUID
		if (vd.getProp<type>(a) == FLUID)
		{
			// Get the position xb of the particle a
			Point<DIM, double> xa = vd.getPos(a);

			// initialize sum
			double W_sum = 0.0;
			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

			// iterate the neighborhood particles
			while (Np.isNext() == true)
			{
				// Key of b particle
				const unsigned long b = Np.get();

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

				// If the particles interact ...
				if (r2 < r_threshold * r_threshold)
				{
					if (vd.getProp<type>(b) == FLUID)
					{
						// calculate distance
						const double r = std::sqrt(r2);

						// evaluate kernel
						const double w = Wab(r);

						W_sum += w;
					}
					else
					{
						if (BC_TYPE == NO_SLIP)
						{
							// calculate distance
							const double r = std::sqrt(r2);

							// evaluate kernel
							const double w = Wab(r);

							W_sum += w;
						}
						else if (BC_TYPE == NEW_NO_SLIP) // need to evaluate kernel at dummy particles
						{
							// xw.get(0) this is the x coordinate of the wall
							Point<DIM, double> normal = vd.getProp<normal_vector>(b);

							// Apply offsets to dr to get 3 vectrors pointing to dummy particles
							std::array<Point<DIM, double>, 3> R_dummy = GetBoundaryPositions(-1.0 * dr, normal);

							const Point<DIM, double> r1 = R_dummy[0];
							const Point<DIM, double> r2 = R_dummy[1];
							const Point<DIM, double> r3 = R_dummy[2];

							const double W1 = Wab(getVectorNorm(r1));
							const double W2 = Wab(getVectorNorm(r2));
							const double W3 = Wab(getVectorNorm(r3));

							W_sum += (W1 + W2 + W3);
						}
					}
				}

				++Np;
			}

			global_sum += W_sum;
			++count;
		}

		++part;
	}
	Vcluster<> &v_cl = create_vcluster();
	v_cl.sum(count);
	v_cl.sum(global_sum);
	v_cl.execute();

	double avg_sum = global_sum / count;
	double V0 = 1.0 / avg_sum;
	double newmass = rho_zero * V0;
	MassFluid = newmass;
	MassBound = newmass;
}

template <typename CellList>
void calc_density(particles &vd, CellList &NN)
{
	// This function computes the density of particles from the summation of the kernel

	const double max_curvature = 1.0 / (3.0 * dp);
	double max_cos = cos(max_angle * M_PI / 180.0 + M_PI / 2.0); // angle is measured from below the plane of the normal vector
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
			double massa = MassFluid;

			// Get the position xb of the particle a
			Point<DIM, double> xa = vd.getPos(a);

			// initialize density sum
			double rho_sum = 0.0;
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

				// If the particles interact ...
				if (r2 < r_threshold * r_threshold)
				{

					if (vd.getProp<type>(b) == FLUID)
					{
						// calculate distance
						const double r = std::sqrt(r2);

						// evaluate kernel
						const double w = Wab(r);

						rho_sum += w * MassFluid;
					}
					else
					{
						if (BC_TYPE == NO_SLIP)
						{
							// calculate distance
							const double r = sqrt(r2);

							// evaluate kernel
							const double w = Wab(r);

							rho_sum += w * MassBound;
						}
						else if (BC_TYPE == NEW_NO_SLIP) // need to evaluate kernel at dummy particles
						{
							const Point<DIM, double> normal = vd.getProp<normal_vector>(b);
							const double kappa = vd.getProp<curvature_boundary>(b);

							double dist2marker = sqrt(r2);
							std::array<double, 3> lwall = {0.5 * dp, 1.5 * dp, 2.5 * dp};

							// Apply offsets to dr to get 3 vectrors pointing to dummy particles
							const std::array<Point<DIM, double>, 3> R_dummy = GetBoundaryPositions(-1.0 * dr, normal);
							const double dxwall = vd.getProp<arc_length>(b);

							for (int i = 0; i < 3; ++i)
							{
								if (dist2marker + lwall[i] < r_threshold)
								{
									const double Vol = 0.5 * (2.0 * dp + dp * dp * kappa - 2.0 * i * dp * dp * kappa) * dxwall;
									const double mass = Vol * rho_zero;

									const Point<DIM, double> ri = R_dummy[i];

									const double W = Wab(getVectorNorm(ri));

									rho_sum += W * mass;
								}
							}
							// if (dotProduct(normal, dr) > (max_cos / max_curvature) * kappa)
							// // if (dotProduct(normal, dr) > 0.0)
							// {

							// }
						}
					}
				}

				++Np;
			}
			if (rho_sum != 0.0)
			{
				vd.template getProp<rho>(a) = rho_sum;
			}
			else
			{
				std::cout << "WARNING: NO particles around, density summation zero" << std::endl;
				vd.template getProp<rho>(a) = rho_zero;
			}
		}

		++part;
	}
}

template <typename CellList>
void calc_boundary(particles &vd, CellList &NN)
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
				if (r2 < r_threshold * r_threshold)
				{
					// calculate distance
					double r = sqrt(r2);

					// evaluate kernel
					double w = Wab(r);
					// compute v*W
					sum_vW += w * vf;

					// compute Pf +rhof*(g-a) dot rwf
					// at the moment a = 0
					const double dot = dotProduct(dr, gravity_vector);
					sum_pW += w * (Pf + rhof * dot);
					sum_W += w;
				}

				++Np;
			}
			if (sum_W != 0.0)
			{
				// Set the velocity of the boundary particle b ( to be used in the momentum equation to impose BC)
				// We use the v_transport field because boundary particles dont use this array

				for (int xyz = 0; xyz < DIM; ++xyz)
				{
					vd.template getProp<v_transport>(b)[xyz] = 2.0 * vd.template getProp<velocity>(b)[xyz] - sum_vW.get(xyz) / sum_W;
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
				vd.template getProp<rho>(b) = rho_zero;
			}
		}

		++part;
	}
}

// function to know the type of variables currently at auto
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
	const double massw = MassBound;
	const double rhow = vd.getProp<rho>(boundary_key);
	const double Pw = vd.getProp<pressure>(boundary_key);
	const Point<DIM, double> vw = vd.getProp<velocity>(boundary_key);
	const Point<DIM, double> vtf = vd.getProp<v_transport>(fluid_key);
	const Point<DIM, double> vdiff_f = vtf - vf;

	// Get normal vector
	Point<DIM, double> normal = vd.getProp<normal_vector>(boundary_key);
	Point<DIM, double> tangential = (DIM == 2) ? Point<DIM, double>{-normal.get(1), normal.get(0)} : Point<DIM, double>{normal.get(1), normal.get(0), 0.0};

	// get curvature, and arc length
	double kappa = vd.getProp<curvature_boundary>(boundary_key);
	double dxwall = vd.getProp<arc_length>(boundary_key);

	// Project va on tangential and normal directions
	double vt = dotProduct(vf, tangential);
	double vn = dotProduct(vf, normal);
	double vwt = dotProduct(vw, tangential);

	// vertical distance from fluid particle to wall
	double lf = dotProduct(r_fluid_to_wall, normal);
	lf = (lf < 0.0 ? -1.0 * lf : lf); // absolute value

	// Get array of vectors from fluid to boundary particles, get norm, and get normal distance to wall
	std::array<Point<DIM, double>, 3> r_boundary = GetBoundaryPositions(r_fluid_to_wall, normal);
	std::array<double, 3> r_boundary_norm = {getVectorNorm(r_boundary[0]), getVectorNorm(r_boundary[1]), getVectorNorm(r_boundary[2])};
	std::array<double, 3> lwall = {0.5 * dp, 1.5 * dp, 2.5 * dp};

	// Initialize arrays for boundary particles
	std::array<double, 3> p_boundary;
	std::array<double, 3> rho_boundary;
	std::array<Point<DIM, double>, 3> v_boundary;
	std::array<double, 3> Volume_boundary;
	std::array<double, 3> Mass_boundary;
	double g_normal = dotProduct(gravity_vector, normal);

	size_t interact_count = 0;
	const std::array<Point<DIM, double>, DIM> Af = dyadicProduct(rhof * vf, vdiff_f);
	lf = std::max(lf, 0.05 * dp);
	for (int i = 0; i < 3; i++)
	{
		if (dist2marker + lwall[i] < r_threshold)
		{

			// compute volume of boundary particle, this gives density and pressure
			Volume_boundary[i] = 0.5 * (2.0 * dp + dp * dp * kappa - 2.0 * (i + 1.0) * dp * dp * kappa) * dxwall;
			Mass_boundary[i] = Volume_boundary[i] * rho_zero;

			// rho_boundary[i] = MassBound / Volume_boundary[i];
			// p_boundary[i] = EqState_particle(rho_boundary[i]);
			// compute velocity of boundary particle
			// v_boundary[i] = 2.0 * vw - vt * (lwall[i] / lf) * tangential - vn * (lwall[i] / lf) * normal; // vn
			// v_boundary[i] = ((vwt - vt) * (lwall[i] / lf) + vwt) * tangential - vn * normal; // _min_vn

			v_boundary[i] = ((vwt - vt) * (lwall[i] / lf) + vwt) * tangential + vn * normal; // no label
			// v_boundary[i] = dotProduct((vw - vf) * (lwall[i] / lf) + vw, tangential) * tangential + dotProduct(vf, normal) * normal;
			// p_boundary[i] = std::max(0.0, Pf + rhof * g_normal * (lf + lwall[i])); // if negative set to 0
			p_boundary[i] = Pf + rhof * g_normal * (lf + lwall[i]); // if negative set to 0

			rho_boundary[i] = InvEqState_particle(p_boundary[i]);

			// flip sing of r_boundary to get vector pointing from boundary to fluid (Force routines use the vector pointing from b to a)
			r_boundary[i] = -1.0 * r_boundary[i];

			// Evaluate kernel gradient
			const Point<DIM, double> DW = DWab(r_boundary[i], r_boundary_norm[i]);

			// Compute forces
			const Point<DIM, double> v_rel = vf - v_boundary[i];
			const double Va2 = (massf / rhof) * (massf / rhof);
			// const double Vb2 = (Mass_boundary[i] / rho_zero) * (Mass_boundary[i] / rho_zero);
			const double Vb2 = (Mass_boundary[i] / rho_boundary[i]) * (Mass_boundary[i] / rho_boundary[i]);
			// const double Vb2 = Volume_boundary[i] * Volume_boundary[i];

			const Point<DIM, double> ViscosityTerm = Pi_physical(r_boundary[i], getVectorNorm(r_boundary[i]), v_rel, DW);
			const double PressureTerm = PressureForce(rhof, rho_boundary[i], Pf, p_boundary[i]);
			Point<DIM, double> GradATerm = 0.5 * matVec(Af, DW);

			if (accumulate_force) // we accumulate x force on cylinder
			{
				cylinder_force += -1.0 * (Va2 + Vb2) * (PressureTerm * DW.get(0) + ViscosityTerm.get(0) + GradATerm.get(0)) / massf;
			}
			for (int xyz = 0; xyz < DIM; ++xyz)
			{
				vd.getProp<force>(fluid_key)[xyz] += (Va2 + Vb2) * (PressureTerm * DW.get(xyz) + ViscosityTerm.get(xyz) + GradATerm.get(xyz)) / massf;
				vd.getProp<force_transport>(fluid_key)[xyz] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(xyz) / massf;
			}

			if (DENSITY_TYPE == DENSITY_DIFFERENTIAL)
			{
				vd.getProp<drho>(fluid_key) += Mass_boundary[i] * dotProduct(v_rel, DW);
			}
		}
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
	const double massb = MassBound;
	const double rhob = vd.getProp<rho>(boundary_key);
	const double Pb = vd.getProp<pressure>(boundary_key);
	const Point<DIM, double> vb = vd.getProp<velocity>(boundary_key);
	const Point<DIM, double> vb_noslip = vd.getProp<v_transport>(boundary_key);

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
		vd.getProp<force_transport>(fluid_key)[xyz] += -1.0 * (Va2 + Vb2) * Pbackground * DW.get(xyz) / massf;
	}

	if (DENSITY_TYPE == DENSITY_DIFFERENTIAL)
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

	const double massb = MassFluid;
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
		vd.getProp<force_transport>(a_key)[xyz] += -1.0 * (Va2 + Vb2) * Pbackground * DW.get(xyz) / massa;
	}
	if (DENSITY_TYPE == DENSITY_DIFFERENTIAL)
	{
		vd.getProp<drho>(a_key) += massb * dotProduct(v_rel, DW);
	}
}

template <typename CellList>
void calc_forces(particles &vd, CellList &NN, double &cylinder_force, bool calc_drag)
{
	vector_dist_iterator part = vd.getDomainIterator();

	const double max_curvature = 1.0 / (3.0 * dp);
	double max_cos = cos(max_angle * M_PI / 180.0 + M_PI / 2.0); // angle is measured from below the plane of the normal vector

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
			double massa = MassFluid;

			// Get the density and pressure of the particle a
			double rhoa = vd.getProp<rho>(a);
			double Pa = vd.getProp<pressure>(a);

			// Get the Velocity of the particle a
			Point<DIM, double> va = vd.getProp<velocity>(a);

			// Reset the force counter (0 + gravity)
			for (int xyz = 0; xyz < DIM; xyz++)
			{
				vd.template getProp<force>(a)[xyz] = gravity_vector.get(xyz);
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
				if (r2 < r_threshold * r_threshold)
				{
					if (vd.getProp<type>(b) == BOUNDARY)
					{
						if (calc_drag && (xb.get(1) > 2.0 * dp || xb.get(1) < length[1] - 2.0 * dp)) // to exclude solid walls
						{
							accumulate_force = true;
						}
						else
						{
							accumulate_force = false;
						}

						if (BC_TYPE == NO_SLIP)
							interact_fluid_boundary_old(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b, accumulate_force, cylinder_force);
						else if (BC_TYPE == NEW_NO_SLIP)
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
		else // for boundary particles we only need to update their acceleration according to the time law
		{
			for (int xyz = 0; xyz < DIM; xyz++)
			{
				vd.template getProp<force>(a)[xyz] = 0.0;
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

	double dt_u = 0.25 * H / (cbar + abs(Maxvel));
	double dt_visc = 0.25 * H * H / (nu);
	double dt_g = 0.25 * sqrt(H / getVectorNorm(gravity_vector));
	double dt = CFLnumber * std::min({dt_u, dt_visc, dt_g});
	// if (dt < DtMin)
	// 	dt = DtMin;

	return dt;
}

template <typename CellList>
void kick_drift_int(particles &vd, CellList &NN, const double dt, Vcluster<> &v_cl, double &cylinder_force, bool calc_drag)
{
	// particle iterator
	auto part = vd.getDomainIterator();

	const double dt_2 = dt * 0.5;
	vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();

	// For each particle ...
	while (part.isNext())
	{
		// get particle a key
		vect_dist_key_dx a = part.get();

		if (vd.template getProp<type>(a) == BOUNDARY)
		{
			for (int xyz = 0; xyz < DIM; xyz++)
			{
				vd.template getProp<velocity>(a)[xyz] += dt_2 * vd.template getProp<force>(a)[xyz];
				vd.getPos(a)[xyz] += dt * vd.template getProp<velocity>(a)[xyz];
			}
		}
		else
		{
			for (int xyz = 0; xyz < DIM; xyz++)
			{
				vd.template getProp<velocity>(a)[xyz] += dt_2 * vd.template getProp<force>(a)[xyz];
				vd.template getProp<v_transport>(a)[xyz] += dt_2 * vd.template getProp<force_transport>(a)[xyz];
				vd.getPos(a)[xyz] += dt * vd.template getProp<v_transport>(a)[xyz];
			}

			if (DENSITY_TYPE == DENSITY_DIFFERENTIAL)
				vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
		}
		++part;
	}
	// map particles if they have gone outside the domain
	vd.map();
	vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();

	// in density summation we need to update the density after moving the particles
	if (DENSITY_TYPE == DENSITY_SUMMATION)
	{
		calc_density(vd, NN);
		vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();
	}

	// Calculate pressure from the density
	EqState(vd);
	vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();

	if (BC_TYPE == NO_SLIP)
	{
		calc_boundary(vd, NN);
		vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();
	}

	cylinder_force = 0.0;
	calc_forces(vd, NN, cylinder_force, calc_drag);
	if (calc_drag)
	{
		cylinder_force = cylinder_force * MassFluid;
		v_cl.sum(cylinder_force);
		v_cl.execute();
	}

	vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();

	// particle iterator
	auto part2 = vd.getDomainIterator();

	// For each particle ...
	while (part2.isNext())
	{
		// get particle a key
		vect_dist_key_dx a = part2.get();

		// if the particle is boundary skip
		if (vd.template getProp<type>(a) == BOUNDARY)
		{
			++part2;
			continue;
		}
		for (int xyz = 0; xyz < DIM; xyz++)
		{
			vd.template getProp<velocity>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force>(a)[xyz];
			vd.template getProp<v_transport>(a)[xyz] = vd.template getProp<velocity>(a)[xyz] + dt_2 * vd.template getProp<force_transport>(a)[xyz];
		}
		++part2;
	}

	// increment the iteration counter
	cnt++;
}

void SetFilename(std::string &filename, const size_t Nfluid[DIM], const long int Nprc, const std::string custom_string = "")
{
	filename = "";
	// Scenario name
	if (SCENARIO == POISEUILLE)
		filename = "Poiseuille";
	else if (SCENARIO == COUETTE)
		filename = "Couette";
	else if (SCENARIO == HYDROSTATIC)
		filename = "Hydrostatic";
	else if (SCENARIO == CYLINDER_ARRAY)
		filename = "CylinderArray";
	else if (SCENARIO == CYLINDER_LATTICE)
		filename = "CylinderLattice";
	else if (SCENARIO == SQUARE)
		filename = "Square";
	else if (SCENARIO == TRIANGLE)
		filename = "Triangle";
	else if (SCENARIO == TRIANGLE_SYM)
		filename = "TriangleSym";
	else if (SCENARIO == CAVITY)
		filename = "Cavity";
	else if (SCENARIO == STEP)
		filename = "Step";
	else if (SCENARIO == TAYLOR_COUETTE)
		filename = "TaylorCouette";
	// BC name
	if (BC_TYPE == NO_SLIP)
		filename += "_OLD_BC";
	else if (BC_TYPE == NEW_NO_SLIP)
		filename += "_NEW_BC";
	// Kernel name
	// if (KERNEL == CUBIC)
	// 	filename += "_Cubic";
	// else if (KERNEL == QUINTIC)
	// 	filename += "_Quintic";
	// Density name
	if (DENSITY_TYPE == DENSITY_SUMMATION)
		filename += "_Summation";
	else if (DENSITY_TYPE == DENSITY_DIFFERENTIAL)
		filename += "_Differential";

	// Add the size of the simulation and the number of processors
	std::string size_proc_name = std::to_string(Nfluid[0]) + "_" + std::to_string(Nfluid[1]) + "_" + std::to_string(Nprc) + "prc";
	filename += "_" + size_proc_name;
	filename += custom_string;
}

void WriteParameters(size_t Nfluid[3], Vcluster<> &v_cl, std::string customString)
{
	const double Lx = length[0]; // channel height
	const double Ly = length[1]; // channel width
	if (DIM == 3)
		const double Lz = length[2]; // channel length ( +0.5dp due to periodicity)

	const long int Nprc = v_cl.getProcessingUnits();
	SetFilename(filename, Nfluid, Nprc, customString);

	std::string constants_filename = filename + "_PARAMETERS" + ".txt";
	std::ofstream file(constants_filename);

	std::string scenario_str = "";
	if (SCENARIO == POISEUILLE)
		scenario_str = "Poiseuille";
	else if (SCENARIO == COUETTE)
		scenario_str = "Couette";
	else if (SCENARIO == HYDROSTATIC)
		scenario_str = "Hydrostatic";
	else if (SCENARIO == CYLINDER_ARRAY)
		scenario_str = "CylinderArray";
	else if (SCENARIO == CYLINDER_LATTICE)
		scenario_str = "CylinderLattice";
	else if (SCENARIO == SQUARE)
		scenario_str = "Square";
	else if (SCENARIO == TRIANGLE)
		scenario_str = "Triangle";
	else if (SCENARIO == TRIANGLE_SYM)
		scenario_str = "TriangleSym";
	else if (SCENARIO == CAVITY)
		scenario_str = "Cavity";
	else if (SCENARIO == STEP)
		scenario_str = "Step";
	else if (SCENARIO == TAYLOR_COUETTE)
		scenario_str = "TaylorCouette";

	std::string BC_str = "";
	if (BC_TYPE == NO_SLIP)
		BC_str = "OldBC";
	else if (BC_TYPE == NEW_NO_SLIP)
		BC_str = "NewBC";

	// std::string kernel_str = "";
	// if (KERNEL == CUBIC)
	// 	kernel_str = "Cubic";
	// else if (KERNEL == QUINTIC)
	// 	kernel_str = "Quintic";

	std::string density_str = "";
	if (DENSITY_TYPE == DENSITY_SUMMATION)
		density_str = "Summation";
	else if (DENSITY_TYPE == DENSITY_DIFFERENTIAL)
		density_str = "Differential";

	file << "gravity = {" << gravity_vector.get(0) << ", " << gravity_vector.get(1) << "}" << std::endl;
	file << "v_topwall = {" << vw_top.get(0) << ", " << vw_top.get(1) << "}" << std::endl;
	file << "v_bottomwall = {" << vw_bottom.get(0) << ", " << vw_bottom.get(1) << "}" << std::endl;
	file << "Lx, Ly, Lz = " << Lx << ", " << Ly << std::endl;
	file << "dp = " << dp << std::endl;
	file << "H = " << H << std::endl;
	file << "Mass = " << MassFluid << std::endl;
	file << "umax = " << umax << std::endl;
	file << "nu = " << nu << std::endl;
	file << "Re =  " << Re << std::endl;
	file << "rho_0 = " << rho_zero << std::endl;
	file << "gamma = " << gamma_ << std::endl;
	file << "xi = " << xi << std::endl;
	file << "c = " << cbar << std::endl;
	file << "B = " << B << std::endl;
	file << "Pbackground = " << Pbackground << std::endl;
	file << "CFLnumber = " << CFLnumber << std::endl;
	file << "SCENARIO = " << scenario_str << std::endl;
	file << "BC_TYPE = " << BC_str << std::endl;
	// file << "KERNEL = " << kernel_str << std::endl;
	file << "DENSITY_TYPE = " << density_str << std::endl;
}

void AddFlatWallNewBC(particles &vd, const int k0, const int kmax, const Point<DIM, double> Corner, const Point<DIM, double> UnitOffset, const double dx, Point<DIM, double> SolidVelocity = {0.0, 0.0})
{

	for (int k = k0; k < kmax; k++)
	{
		Point<DIM, double> WallPos = Corner + k * UnitOffset;

		vd.add();
		vd.template getLastProp<type>() = BOUNDARY;
		for (int xyz = 0; xyz < DIM; xyz++)
		{
			vd.getLastPos()[xyz] = WallPos.get(xyz); //+ ((double)rand() / RAND_MAX - 0.5) * dp;
			vd.template getLastProp<velocity>()[xyz] = SolidVelocity.get(xyz);
			vd.template getLastProp<force>()[xyz] = 0.0;
			vd.template getLastProp<force_transport>()[xyz] = 0.0;
			vd.template getLastProp<v_transport>()[xyz] = 0.0;
			vd.template getLastProp<normal_vector>()[xyz] = 0.0;
		}

		vd.template getLastProp<pressure>() = 0.0;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<drho>() = 0.0;
		vd.template getLastProp<curvature_boundary>() = 0.0;
		vd.template getLastProp<arc_length>() = dx;
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
	double epsilon = 0.45 * dp;
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
	double epsilon = 0.45 * dp;
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
private:
public:
	virtual bool isInside(Point<DIM, double> P) = 0;
	virtual void AddObstacle(particles &vd) = 0;
};

class EmptyObstacle : public Obstacle
{
public:
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
	Point<DIM, double> Centre_;
	double Radius_;
	Sphere<DIM, double> Cylinder_;
	double AnglularVelocity_;
	Point<DIM, double> LinearVelocity_;

public:
	CylinderObstacle(Point<DIM, double> Centre, double Radius, double omega, Point<DIM, double> vel) : Centre_(Centre), Radius_(Radius), AnglularVelocity_(omega), LinearVelocity_(vel), Cylinder_(Centre, Radius) {}

	bool isInside(Point<DIM, double> P) override
	{
		double radius_aux = Radius_;

		if (BC_TYPE == NEW_NO_SLIP)
			radius_aux = Radius_ + 0.1 * dp;

		Sphere<DIM, double> Cylinderaux(Centre_, radius_aux);
		return Cylinderaux.isInside(P);
	}
	bool isOutside(Point<DIM, double> P) // for outer cylinder in taylor couette
	{
		double radius_aux = Radius_;
		if (BC_TYPE == NEW_NO_SLIP)
			radius_aux = Radius_ - 0.1 * dp;
		Sphere<DIM, double> Cylinderaux(Centre_, radius_aux);
		return Cylinderaux.isInside(P);
	}

	void AddObstacle(particles &vd)
	{
		const double perimeter = 2.0 * M_PI * Radius_;
		const int Np_cylinder = ceil(perimeter / dp);
		const double dtheta = 2.0 * M_PI / Np_cylinder;
		const double dxwall = dtheta * Radius_;
		double theta = 0.0;
		const double rotation_speed = AnglularVelocity_ * Radius_;

		Point<DIM, double> Cylinder_particle;

		for (int k = 0; k < Np_cylinder; k++)
		{

			const Point<DIM, double> tangential = {-sin(theta), cos(theta)};

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
				vd.template getLastProp<force_transport>()[xyz] = 0.0;
				vd.template getLastProp<v_transport>()[xyz] = 0.0;
				vd.template getLastProp<normal_vector>()[xyz] = 0.0;
				vd.template getLastProp<velocity>()[xyz] = LinearVelocity_.get(xyz) + rotation_speed * tangential.get(xyz);
			}

			vd.template getLastProp<pressure>() = 0.0;
			vd.template getLastProp<rho>() = rho_zero;
			vd.template getLastProp<drho>() = 0.0;
			vd.template getLastProp<curvature_boundary>() = 0.0; // 1.0 / Radius_;
			vd.template getLastProp<arc_length>() = dxwall;
			theta += dtheta;
		}
	}
};

class RectangleObstacle : public Obstacle
{
private:
	const Point<DIM, double> Centre_;
	const unsigned int BaseLength_;
	const unsigned int HeigthLength_;
	Box<DIM, double> Rectangle_;

	Point<DIM, double> LowerLeft_;
	Point<DIM, double> LowerRight_;
	Point<DIM, double> UpperLeft_;
	Point<DIM, double> UpperRight_;

public:
	RectangleObstacle(Point<DIM, double> Centre, unsigned int BaseLength, unsigned int HeigthLength) : Centre_(Centre), BaseLength_(BaseLength), HeigthLength_(HeigthLength)
	{
		LowerLeft_ = Point<DIM, double>{Centre.get(0) - (((double)BaseLength_ - 1.0) * dp) / 2.0,
										Centre.get(1) - (((double)HeigthLength_ - 1.0) * dp) / 2.0};
		LowerRight_ = Point<DIM, double>{Centre.get(0) + (((double)BaseLength_ - 1.0) * dp) / 2.0,
										 Centre.get(1) - (((double)HeigthLength_ - 1.0) * dp) / 2.0};
		UpperLeft_ = Point<DIM, double>{Centre.get(0) - (((double)BaseLength_ - 1.0) * dp) / 2.0,
										Centre.get(1) + (((double)HeigthLength_ - 1.0) * dp) / 2.0};
		UpperRight_ = Point<DIM, double>{Centre.get(0) + (((double)BaseLength_ - 1.0) * dp) / 2.0,
										 Centre.get(1) + (((double)HeigthLength_ - 1.0) * dp) / 2.0};

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
		AddFlatWallNewBC(vd, 0, BaseLength_, LowerLeft_, Xoffset, dp);
		// Right wall
		AddFlatWallNewBC(vd, 1, HeigthLength_ - 1, LowerRight_, Yoffset, dp);
		// Upper wall
		AddFlatWallNewBC(vd, 0, BaseLength_, UpperRight_, -1.0 * Xoffset, dp);
		// Left wall
		AddFlatWallNewBC(vd, 1, HeigthLength_ - 1, UpperLeft_, -1.0 * Yoffset, dp);
	}
};

class TriangleObstacle : public Obstacle
{
private:
	const Point<DIM, double> Centre_;
	const unsigned int BaseLength_;
	const unsigned int HeigthLength_;
	Box<DIM, double> ContainingRectangle_;
	Point<DIM, double> LowerLeft_;
	Point<DIM, double> LowerRight_;
	Point<DIM, double> UpperRight_;

public:
	TriangleObstacle(Point<DIM, double> Centre, unsigned int BaseLength, unsigned int HeigthLength) : Centre_(Centre), BaseLength_(BaseLength), HeigthLength_(HeigthLength)
	{
		LowerLeft_ = Point<DIM, double>{Centre.get(0) - (((double)BaseLength_ - 1.0) * dp) / 2.0,
										Centre.get(1) - (((double)HeigthLength_ - 1.0) * dp) / 2.0};
		LowerRight_ = Point<DIM, double>{Centre.get(0) + (((double)BaseLength_ - 1.0) * dp) / 2.0,
										 Centre.get(1) - (((double)HeigthLength_ - 1.0) * dp) / 2.0};
		UpperRight_ = Point<DIM, double>{Centre.get(0) + (((double)BaseLength_ - 1.0) * dp) / 2.0,
										 Centre.get(1) + (((double)HeigthLength_ - 1.0) * dp) / 2.0};

		ContainingRectangle_ = Box<DIM, double>(LowerLeft_, UpperRight_);
	}

	bool isInside(Point<DIM, double> P) override
	{
		return (ContainingRectangle_.isInside(P) && !isAvobeLine(LowerLeft_, UpperRight_, P));
	}

	void AddObstacle(particles &vd)
	{
		const double baseLength = (BaseLength_ - 1) * dp;
		const double heigthLength = (HeigthLength_ - 1) * dp;

		const Point<DIM, double> Xoffset = {dp, 0.0};
		const Point<DIM, double> Yoffset = {0.0, dp};

		// Lower wall
		AddFlatWallNewBC(vd, 0, BaseLength_, LowerLeft_, Xoffset, dp);
		// Right wall
		AddFlatWallNewBC(vd, 1, HeigthLength_, LowerRight_, Yoffset, dp);

		// We want particles spaced roughly by dp
		const double hypothenuseLength = sqrt(baseLength * baseLength + heigthLength * heigthLength);
		const int Ndiag = ceil(hypothenuseLength / dp);	 // integer number of particles that can fit in the diagonal
		const double dxwall = hypothenuseLength / Ndiag; // actual spacing between particles ( close to dp but not exactly)
		const double sin_theta = heigthLength / hypothenuseLength;
		const double cos_theta = baseLength / hypothenuseLength;
		const Point<DIM, double> Diagoffset{dxwall * cos_theta, dxwall * sin_theta};

		Point<DIM, double> h_normal = {-1.0 * sin_theta, cos_theta};
		//  Hypothenuse wall
		AddFlatWallNewBC(vd, 1, Ndiag, UpperRight_, -1.0 * Diagoffset, dxwall);
	}
};

class TriangleSymObstacle : public Obstacle
{
private:
	const Point<DIM, double> Centre_;
	const unsigned int BaseLength_;
	const unsigned int HeigthLength_;
	Box<DIM, double> ContainingRectangle_;

	Point<DIM, double> LowerLeft_;
	Point<DIM, double> UpperRight_;
	Point<DIM, double> LowerRight_;
	Point<DIM, double> TriangleTip_;

public:
	TriangleSymObstacle(Point<DIM, double> Centre, unsigned int BaseLength, unsigned int HeigthLength) : Centre_(Centre), BaseLength_(BaseLength), HeigthLength_(HeigthLength)
	{
		LowerLeft_ = Point<DIM, double>{Centre.get(0) - (((double)BaseLength_ - 1.0) * dp) / 2.0,
										Centre.get(1) - (((double)HeigthLength_ - 1.0) * dp)};
		UpperRight_ = Point<DIM, double>{Centre.get(0) + (((double)BaseLength_ - 1.0) * dp) / 2.0,
										 Centre.get(1) + (((double)HeigthLength_ - 1.0) * dp)};
		LowerRight_ = Point<DIM, double>{Centre.get(0) + (((double)BaseLength_ - 1.0) * dp) / 2.0,
										 Centre.get(1) - (((double)HeigthLength_ - 1.0) * dp)};
		TriangleTip_ = Point<DIM, double>{Centre.get(0) - (((double)BaseLength_ - 1.0) * dp) / 2.0,
										  Centre.get(1)};

		ContainingRectangle_ = Box<DIM, double>(LowerLeft_, UpperRight_);
	}
	bool isInside(Point<DIM, double> P) override
	{
		return (ContainingRectangle_.isInside(P) && !isAvobeLine(TriangleTip_, UpperRight_, P) && !isBelowLine(TriangleTip_, LowerRight_, P));
	}

	void AddObstacle(particles &vd)
	{
		const double baseLength = (BaseLength_ - 1) * dp;
		const double heigthLength = (HeigthLength_ - 1) * dp;

		Point<DIM, double> Xoffset = {dp, 0.0};
		Point<DIM, double> Yoffset = {0.0, dp};

		// Lower right wall
		AddFlatWallNewBC(vd, 0, HeigthLength_, LowerRight_, Yoffset, dp);
		// Right wall
		AddFlatWallNewBC(vd, 0, HeigthLength_ - 1, UpperRight_, -1.0 * Yoffset, dp);

		// We want particles spaced roughly by dp
		const double hypothenuseLength = sqrt(baseLength * baseLength + heigthLength * heigthLength);
		const int Ndiag = ceil(hypothenuseLength / dp);	 // integer number of particles that can fit in the diagonal
		const double dxwall = hypothenuseLength / Ndiag; // actual spacing between particles ( close to dp but not exactly)
		const double sin_theta = heigthLength / hypothenuseLength;
		const double cos_theta = baseLength / hypothenuseLength;
		const Point<DIM, double> Diagoffset{dxwall * cos_theta, dxwall * sin_theta};
		Point<DIM, double> DiagoffsetNW = Diagoffset;
		DiagoffsetNW.get(0) = -1.0 * DiagoffsetNW.get(0);
		//  Hypothenuse upper wall
		AddFlatWallNewBC(vd, 1, Ndiag + 1, UpperRight_, -1.0 * Diagoffset, dxwall);
		// Hypothenuse lower wall
		AddFlatWallNewBC(vd, 1, Ndiag, LowerRight_, DiagoffsetNW, dxwall);
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
		if (r_boundary_norm[i] < r_threshold)
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
					else if (BC_TYPE == NO_SLIP)
					{
						double ker = Wab(r) * (MassBound / rho_zero);
						W_sum += ker;
						magnitude_tmp += vd.template getProp<v_transport>(q)[component] * ker;
					}
				}
				else if (vd.template getProp<type>(q) == FLUID)
				{
					double rhoq = vd.template getProp<rho>(q);
					double ker = Wab(r) * (MassFluid / rhoq);

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
	cylinder_force = cylinder_force / (eta * v);

	if (v_cl.getProcessUnitID() == 0)
	{
		avgvelstream << t << ", " << v << ", " << cylinder_force << std::endl;
	}
}

void CreateParticleGeometry(particles &vd, std::vector<std::pair<probe_particles, int>> &vp_vec, Vcluster<> &v_cl, Obstacle *&obstacle_ptr)
{
	// Size of the virtual cartesian grid that defines where to place the particles
	size_t sz[DIM];

	// In the case of the new bc we need particles at the wall, for this we need sz_aux
	// We want to put one virtual grid point between each pair of the old ones,
	// so that the new spacing is dp/2, and we can put a fluid particle exactly at the wall
	size_t sz_aux[DIM];

	// Boundary conditions
	size_t bc[DIM];

	// Number of boundary particles in each direction
	size_t Nboundary[DIM];

	// Number of fluid particles in each direction
	size_t Nfluid[DIM];

	// We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
	double offset_domain_left[DIM] = {0.0};
	double offset_domain_right[DIM] = {0.0};
	double offset_recipient[DIM] = {0.0};
	double offset_periodic_fluid[DIM] = {0.0};
	double offset_periodic_recipient[DIM] = {0.0};

	// for (int xyz = 0; xyz < DIM; xyz++)
	// {
	// 	offset_domain[0] = 0.0;
	// 	offset_recipient[0] = 0.0;
	// 	offset_periodic_fluid[0] = 0.0;
	// 	offset_periodic_recipient[0] = 0.0;
	// }

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

	double CylinderRadius;

	t_end = 100.0;
	write_const = 10;
	if (SCENARIO == POISEUILLE)
	{
		Nfluid[0] = 20;
		Nfluid[1] = 40;
		size_t Nbound = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = 0;
		Nboundary[1] = Nbound;
		bc[0] = PERIODIC;
		bc[1] = NON_PERIODIC;
		LengthScale = 1.0;
		dp = LengthScale / Nfluid[1];
		rho_zero = 1.0;
		nu = 0.01;
		Bfactor = 1.0;
		gravity_vector.get(0) = 0.1;
		umax = gravity_vector.get(0) * LengthScale * LengthScale / (8.0 * nu);
		t_end = 100.0;
		write_const = 10;
		PROBES_ENABLED = 0;
	}
	else if (SCENARIO == COUETTE)
	{
		Nfluid[0] = 20;
		Nfluid[1] = 40;
		size_t Nbound = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = 0;
		Nboundary[1] = Nbound;
		bc[0] = PERIODIC;
		bc[1] = NON_PERIODIC;
		LengthScale = 1.0;
		dp = LengthScale / Nfluid[1];
		rho_zero = 1.0;
		nu = 0.01;
		Bfactor = 1.0;
		// gravity_vector.get(0) = 0.1;
		vw_top.get(0) = 1.25;
		umax = vw_top.get(0);
		t_end = 100.0;
		write_const = 10;
		PROBES_ENABLED = 0;
	}
	else if (SCENARIO == HYDROSTATIC)
	{
		Nfluid[0] = 40;
		Nfluid[1] = 40;
		size_t Nbound = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = Nbound;
		Nboundary[1] = Nbound;
		bc[0] = NON_PERIODIC;
		bc[1] = NON_PERIODIC;
		LengthScale = 1.0;
		dp = LengthScale / Nfluid[0];
		rho_zero = 1.0;
		nu = 0.1;
		Bfactor = 1.0;
		gravity_vector.get(1) = 0.1;
		umax = 1.0;
		t_end = 100.0;
		write_const = 10;
		PROBES_ENABLED = 0;
	}
	else if (SCENARIO == CAVITY)
	{
		Nfluid[0] = 80;
		Nfluid[1] = 80;
		size_t Nbound = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = Nbound;
		Nboundary[1] = Nbound;
		bc[0] = NON_PERIODIC;
		bc[1] = NON_PERIODIC;
		LengthScale = 1.0;
		dp = LengthScale / Nfluid[0];
		rho_zero = 1.0;
		nu = 0.01;
		Bfactor = 3.0;
		vw_top.get(0) = 1.0;
		umax = vw_top.get(0);
		t_end = 50.0;
		write_const = 10;
		PROBES_ENABLED = 1;
	}
	else if (SCENARIO == CYLINDER_ARRAY)
	{
		Nfluid[0] = 72;
		Nfluid[1] = 48;
		size_t Nbound = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = 0;
		Nboundary[1] = Nbound;
		bc[0] = PERIODIC;
		bc[1] = NON_PERIODIC;
		CylinderRadius = 0.02;
		LengthScale = CylinderRadius;
		dp = 4.0 * CylinderRadius / (double)Nfluid[1]; // chanel height is 4 times the cylinder radius
		rho_zero = 1000.0;
		nu = 1e-4;
		Bfactor = 3.0;
		gravity_vector.get(0) = 2.5 * 1e-4;
		// umax = 2.4037 * 1e-1; // from previous simulations for nu = 0.01
		umax = 2.5e-3;
		t_end = 1000.0;
		write_const = 1;
		PROBES_ENABLED = 0;
	}
	else if (SCENARIO == CYLINDER_LATTICE)
	{
		Nfluid[0] = 100;
		Nfluid[1] = 100;
		Nboundary[0] = 0;
		Nboundary[1] = 0;
		bc[0] = PERIODIC;
		bc[1] = PERIODIC;
		CylinderRadius = 0.02;
		LengthScale = CylinderRadius;
		dp = 5.0 * CylinderRadius / (double)Nfluid[1]; // chanel height is 5 times the cylinder radius
		rho_zero = 1000.0;
		nu = 1e-6;
		Bfactor = 3.0;
		gravity_vector.get(0) = 1.5 * 1e-7;
		umax = 1.2e-4;
		// umax = 5.77 * 1e-5; // (morris, to get c=5.77*1e-4)
		t_end = 10000.0;
		write_const = 0.1;
		PROBES_ENABLED = 1;
	}
	else if (SCENARIO == SQUARE)
	{
		Nfluid[0] = 60;
		Nfluid[1] = 40;
		size_t Nbound = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = 0;
		Nboundary[1] = Nbound;
		bc[0] = PERIODIC;
		bc[1] = NON_PERIODIC;
		LengthScale = 1.0;
		dp = LengthScale / Nfluid[1];
		rho_zero = 1.0;
		nu = 0.01;
		Bfactor = 3.0;
		gravity_vector.get(0) = 0.1;
		umax = 4.1 * 1e-1;
		t_end = 40.0;
		write_const = 100;
		PROBES_ENABLED = 0;
	}
	else if (SCENARIO == TRIANGLE)
	{
		Nfluid[0] = 60;
		Nfluid[1] = 40;
		size_t Nbound = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = 0;
		Nboundary[1] = Nbound;
		bc[0] = PERIODIC;
		bc[1] = NON_PERIODIC;
		LengthScale = 1.0;
		dp = LengthScale / (double)Nfluid[1];
		rho_zero = 1.0;
		nu = 0.01;
		Bfactor = 3.0;
		gravity_vector.get(0) = 0.1;
		umax = 4.7 * 1e-1; // from previous simulations for nu = 0.01
		t_end = 100.0;
		write_const = 10;
		PROBES_ENABLED = 0;
	}
	else if (SCENARIO == TRIANGLE_SYM)
	{
		Nfluid[0] = 60;
		Nfluid[1] = 40;
		size_t Nbound = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = 0;
		Nboundary[1] = Nbound;
		bc[0] = PERIODIC;
		bc[1] = NON_PERIODIC;
		LengthScale = 1.0;
		dp = LengthScale / (double)Nfluid[1];
		rho_zero = 1.0;
		nu = 0.01;
		Bfactor = 3.0;
		gravity_vector.get(0) = 0.1;
		umax = 3.7 * 1e-1; // from previous simulations for nu = 0.01
		t_end = 100.0;
		write_const = 10;
		PROBES_ENABLED = 0;
	}

	H = Hconst * dp;
	// r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
	r_threshold = 3.0 * H;
	// Kcubic = (DIM == 3) ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;
	Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / H / H / H : 7.0 / 478.0 / M_PI / H / H;
	MassFluid = rho_zero * (DIM == 3 ? dp * dp * dp : dp * dp);
	MassBound = rho_zero * (DIM == 3 ? dp * dp * dp : dp * dp);
	cbar = coeff_sound * umax;
	B = rho_zero * cbar * cbar / gamma_;
	Pbackground = Bfactor * B;
	eta = nu * rho_zero;
	Re = umax * LengthScale / nu;

	gravity = getVectorNorm(gravity_vector);

	for (int dim = 0; dim < DIM; dim++)
	{
		if (bc[dim] == NON_PERIODIC) // non periodic, fluid covered by boundary
		{
			length[dim] = dp * Nfluid[dim];
			sz[dim] = Nfluid[dim] + 2 * (Nboundary[dim] + 1);
			offset_domain_left[dim] = (0.5 + Nboundary[dim]) * dp;
			offset_domain_right[dim] = (0.5 + Nboundary[dim]) * dp;

			if (Nboundary[dim] != 0)
				sz_aux[dim] = 2 * Nfluid[dim] - 1 + 2 * (2 * Nboundary[dim] + 1 + 1);
			else // for a direction with no boundary particles we dont need to add anything
				sz_aux[dim] = sz[dim];

			if (BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
				offset_recipient[dim] = 0.1 * Nboundary[dim] * dp;
			else if (BC_TYPE == NO_SLIP)
				offset_recipient[dim] = Nboundary[dim] * dp;
		}
		else // periodic, open ended
		{
			Nfluid[dim] -= 1;
			length[dim] = dp * Nfluid[dim];

			sz[dim] = Nfluid[dim] + 2;
			sz_aux[dim] = sz[dim];

			offset_domain_left[dim] = 0.0;
			offset_domain_right[dim] = dp;
			offset_periodic_fluid[dim] = 0.75 * dp;
			offset_periodic_recipient[dim] = 0.85 * dp;
		}
	}

	// Define the boxes
	Box<DIM, double> domain({-offset_domain_left[0],
							 -offset_domain_left[1]},
							{length[0] + offset_domain_right[0],
							 length[1] + offset_domain_right[1]});

	Box<DIM, double> fluid_box({0.0,
								0.0},
							   {length[0] + offset_periodic_fluid[0],
								length[1] + offset_periodic_fluid[1]});

	Box<DIM, double> recipient({-offset_recipient[0],
								-offset_recipient[1]},
							   {length[0] + offset_recipient[0] + offset_periodic_recipient[0],
								length[1] + offset_recipient[1] + offset_periodic_recipient[1]});

	// Will only be used in the new bc
	Box<DIM, double> recipient_hole({offset_recipient[0],
									 offset_recipient[1]},
									{length[0] - offset_recipient[0] + offset_periodic_fluid[0],
									 length[1] - offset_recipient[1] + offset_periodic_fluid[1]});

	// extended boundary around the domain, and the processor domain
	Ghost<DIM, double> g(r_threshold);

	// create particle object
	particles vd_loc(0, domain, bc, g, DEC_GRAN(512));
	vd = vd_loc;

	// correct the number of particles in case of periodicity, we substracted 1 before to accomodate the periodic boundary
	for (int dim = 0; dim < DIM; dim++)
	{
		if (bc[dim] == PERIODIC)
		{
			Nfluid[dim] += 1;
			length[dim] += dp;
		}
	}

	// Write constants on file
	std::string customString = "";
	WriteParameters(Nfluid, v_cl, customString);

	// place probes
	if (PROBES_ENABLED)
	{
		// we want to place probes  in a vertical line at this locations
		Point<DIM, double> EndChannel = {length[0], 0.0};
		Point<DIM, double> HalfChannel = {length[0] / 2.0, 0.0};
		Point<DIM, double> HalfHeight = {0.0, length[1] / 2.0};

		std::vector<Point<DIM, double>> ProbePoints;
		std::vector<int> ProbeComponents;

		Point<DIM, double> VerticalOffset = {0.0, dp};
		Point<DIM, double> HorizontalOffset = {dp, 0.0};
		std::vector<Point<DIM, double>> Offsets;

		int k0 = 0;
		int kendHeight = Nfluid[1] + 1;
		int kendWidth = Nfluid[0] + 1;
		std::vector<int> maxIters;

		if (SCENARIO == CAVITY)
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
				openfpm::vector<std::string> names_p;
				names_p.add("vx");
				vp_loc.setPropNames(names_p);
			}
			else if (ProbeComponents[k] == 1)
			{
				openfpm::vector<std::string> names_p;
				names_p.add("vy");
				vp_loc.setPropNames(names_p);
			}

			if (v_cl.getProcessUnitID() == 0)
			{
				PlaceProbes(vp_loc, k0, maxIters[k], ProbePoints[k], Offsets[k]);
			}
			std::pair<probe_particles, int> tmp = std::make_pair(vp_loc, ProbeComponents[k]);
			vp_vec.push_back(tmp);
			probe_filenames.push_back("probes_" + std::to_string(k) + "_" + filename);
		}
	}

	// Set cylindrical object parameters
	Point<DIM, double> CylinderCentre;
	CylinderCentre.get(0) = length[0] / 2.0;
	CylinderCentre.get(1) = length[1] / 2.0;

	// Set square/triangle obstacle parameters
	const Point<DIM, double> RectangleCentre = CylinderCentre;
	const int integerBaseLength = 15;
	const int integerHeigthLength = 8;

	if (SCENARIO == CYLINDER_ARRAY || SCENARIO == CYLINDER_LATTICE)
		obstacle_ptr = new CylinderObstacle(CylinderCentre, CylinderRadius, 0.0, {-4e-4, 0.0});
	else if (SCENARIO == SQUARE)
		obstacle_ptr = new RectangleObstacle(RectangleCentre, integerBaseLength, integerHeigthLength);
	else if (SCENARIO == TRIANGLE)
		obstacle_ptr = new TriangleObstacle(RectangleCentre, integerBaseLength, integerHeigthLength);
	else if (SCENARIO == TRIANGLE_SYM)
		obstacle_ptr = new TriangleSymObstacle(RectangleCentre, integerBaseLength, integerHeigthLength);
	else
		obstacle_ptr = new EmptyObstacle();

	// Add the obstacle as marker particles only on processor 0
	if (BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
	{
		obstacle_ptr->AddObstacle(vd);
	}

	// return an iterator to the fluid particles to add to vd
	auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

	// for each particle inside the fluid box ...
	while (fluid_it.isNext())
	{

		Point<DIM, double> iterator_position = fluid_it.get();

		if ((*obstacle_ptr).isInside(iterator_position)) // if inside the obstacle region
		{
			if (BC_TYPE == NO_SLIP) // add particle but set it as boundary
			{
				// ... add a particle ...
				vd.add();
				vd.template getLastProp<type>() = BOUNDARY;
			}
			else if (BC_TYPE == NEW_NO_SLIP) // not add particle because already added
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
		}
		// Set properties
		vd.template getLastProp<rho>() = rho_zero;
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
		++fluid_it;
	}

	// Now place solid walls
	openfpm::vector<Box<DIM, double>> holes;

	if (BC_TYPE == NEW_NO_SLIP)
	{
		holes.add(recipient_hole);
		sz[0] = sz_aux[0];
		sz[1] = sz_aux[1];
	}
	else if (BC_TYPE == NO_SLIP)
		holes.add(fluid_box);

	Box<DIM, double> hole_get = holes.get(0);
	auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);

	if (bc[0] != PERIODIC || bc[1] != PERIODIC) // no walls in all periodic scenario
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
			if (hole_get.isInside(position))
			{
				++bound_box;
				continue;
			}

			if (BC_TYPE == NEW_NO_SLIP && (bc[0] == NON_PERIODIC && bc[1] == NON_PERIODIC))
			{
				// Check if x and z coordinates are multiples of dp, keep multiples, discard the rest
				double remx = fmod(position.get(0), dp);
				double remz = fmod(position.get(1), dp);
				double tol = 0.5 * dp * 10e-2;

				if (remx > tol && remx < dp - tol)
				{
					++bound_box;
					continue;
				}
				if (remz > tol && remz < dp - tol)
				{
					++bound_box;
					continue;
				}
			}
			vd.add();

			vd.template getLastProp<type>() = BOUNDARY;
			vd.template getLastProp<rho>() = rho_zero;
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
					vd.template getLastProp<velocity>()[xyz] = vw_bottom.get(xyz);
				}
				else if (position.get(1) > length[1] - dp / 4.0) // top wall
				{
					vd.template getLastProp<velocity>()[xyz] = vw_top.get(xyz);
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
										"arc_length"});
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

	// Boundary conditions
	size_t bc[DIM];

	// Number of boundary particles in each direction
	size_t Nboundary[DIM];

	// Number of fluid particles in each direction
	size_t Nfluid[DIM];

	// We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
	double offset_domain_left[DIM] = {0.0};
	double offset_domain_right[DIM] = {0.0};
	double offset_recipient[DIM] = {0.0};
	double offset_periodic_fluid[DIM] = {0.0};
	double offset_periodic_recipient[DIM] = {0.0};

	// for (int xyz = 0; xyz < DIM; xyz++)
	// {
	// 	offset_domain[0] = 0.0;
	// 	offset_recipient[0] = 0.0;
	// 	offset_periodic_fluid[0] = 0.0;
	// 	offset_periodic_recipient[0] = 0.0;
	// }

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

	double Rin = 1.5;
	double Rout = 2.0;
	double Win = 0.1;
	double Wout = -0.075;
	double Vin = Win * Rin;
	double Vout = Wout * Rout;

	int Nfluid_in = 20;
	Nfluid[0] = 80;
	Nfluid[1] = 80;
	size_t Nbound = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
	Nboundary[0] = 0;
	Nboundary[1] = 0;
	bc[0] = NON_PERIODIC;
	bc[1] = NON_PERIODIC;
	LengthScale = 2.0 * Rout;
	dp = LengthScale / (double)Nfluid[1];
	rho_zero = 1.0;
	nu = 0.1;
	Bfactor = 5.0;
	gravity_vector.get(0) = 0.0;
	umax = 1.0 * abs(Vout);
	t_end = 100.0;
	write_const = 10;
	PROBES_ENABLED = 0;

	H = Hconst * dp;
	// r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
	r_threshold = 3.0 * H;
	// Kcubic = (DIM == 3) ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;
	Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / H / H / H : 7.0 / 478.0 / M_PI / H / H;
	MassFluid = rho_zero * (DIM == 3 ? dp * dp * dp : dp * dp);
	MassBound = rho_zero * (DIM == 3 ? dp * dp * dp : dp * dp);
	cbar = coeff_sound * umax;
	B = rho_zero * cbar * cbar / gamma_;
	Pbackground = Bfactor * B;
	eta = nu * rho_zero;
	Re = umax * LengthScale / nu;

	gravity = getVectorNorm(gravity_vector);

	for (int dim = 0; dim < DIM; dim++)
	{
		if (bc[dim] == NON_PERIODIC) // non periodic, fluid covered by boundary
		{
			length[dim] = dp * Nfluid[dim];
			sz[dim] = Nfluid[dim] + 2 * (Nboundary[dim] + 1);
			offset_domain_left[dim] = (0.5 + Nboundary[dim]) * dp;
			offset_domain_right[dim] = (0.5 + Nboundary[dim]) * dp;

			if (Nboundary[dim] != 0)
				sz_aux[dim] = 2 * Nfluid[dim] - 1 + 2 * (2 * Nboundary[dim] + 1 + 1);
			else // for a direction with no boundary particles we dont need to add anything
				sz_aux[dim] = sz[dim];

			if (BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
				offset_recipient[dim] = 0.1 * Nboundary[dim] * dp;
			else if (BC_TYPE == NO_SLIP)
				offset_recipient[dim] = Nboundary[dim] * dp;
		}
		else // periodic, open ended
		{
			Nfluid[dim] -= 1;
			length[dim] = dp * Nfluid[dim];

			sz[dim] = Nfluid[dim] + 2;
			sz_aux[dim] = sz[dim];

			offset_domain_left[dim] = 0.0;
			offset_domain_right[dim] = dp;
			offset_periodic_fluid[dim] = 0.75 * dp;
			offset_periodic_recipient[dim] = 0.85 * dp;
		}
	}

	// Define the boxes
	Box<DIM, double> domain({-Rout - dp,
							 -Rout - dp},
							{Rout + dp,
							 Rout + dp});

	Box<DIM, double> fluid_box({-Rout + 0.1 * dp,
								-Rout + 0.1 * dp},
							   {Rout + 0.1 * dp,
								Rout + 0.1 * dp});

	// extended boundary around the domain, and the processor domain
	Ghost<DIM, double> g(r_threshold);

	// create particle object
	particles vd_loc(0, domain, bc, g, DEC_GRAN(512));
	vd = vd_loc;

	// correct the number of particles in case of periodicity, we substracted 1 before to accomodate the periodic boundary
	for (int dim = 0; dim < DIM; dim++)
	{
		if (bc[dim] == PERIODIC)
		{
			Nfluid[dim] += 1;
			length[dim] += dp;
		}
	}

	// Write constants on file
	std::string customString = "";
	WriteParameters(Nfluid, v_cl, customString);

	// place probes
	if (PROBES_ENABLED)
	{
		// we want to place probes  in a vertical line at this locations
		Point<DIM, double> EndChannel = {length[0], 0.0};
		Point<DIM, double> HalfChannel = {length[0] / 2.0, 0.0};
		Point<DIM, double> HalfHeight = {0.0, length[1] / 2.0};

		std::vector<Point<DIM, double>> ProbePoints;
		std::vector<int> ProbeComponents;

		Point<DIM, double> VerticalOffset = {0.0, dp};
		Point<DIM, double> HorizontalOffset = {dp, 0.0};
		std::vector<Point<DIM, double>> Offsets;

		int k0 = 0;
		int kendHeight = Nfluid[1] + 1;
		int kendWidth = Nfluid[0] + 1;
		std::vector<int> maxIters;

		if (SCENARIO == CAVITY)
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
				openfpm::vector<std::string> names_p;
				names_p.add("vx");
				vp_loc.setPropNames(names_p);
			}
			else if (ProbeComponents[k] == 1)
			{
				openfpm::vector<std::string> names_p;
				names_p.add("vy");
				vp_loc.setPropNames(names_p);
			}

			if (v_cl.getProcessUnitID() == 0)
			{
				PlaceProbes(vp_loc, k0, maxIters[k], ProbePoints[k], Offsets[k]);
			}
			std::pair<probe_particles, int> tmp = std::make_pair(vp_loc, ProbeComponents[k]);
			vp_vec.push_back(tmp);
			probe_filenames.push_back("probes_" + std::to_string(k) + "_" + filename);
		}
	}

	// Set cylindrical object parameters
	Point<DIM, double> CylinderCentre;
	CylinderCentre.get(0) = 0.0;
	CylinderCentre.get(1) = 0.0;
	obstacle_ptr = new EmptyObstacle();
	CylinderObstacle *obstacle_ptr_out = new CylinderObstacle(CylinderCentre, Rout, Wout, 0.0);
	CylinderObstacle *obstacle_ptr_in = new CylinderObstacle(CylinderCentre, Rin, Win, 0.0);

	// Add the obstacle as marker particles only on processor 0
	if (BC_TYPE == NEW_NO_SLIP && v_cl.getProcessUnitID() == 0)
	{
		obstacle_ptr_out->AddObstacle(vd);
		obstacle_ptr_in->AddObstacle(vd);
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
				if (BC_TYPE == NO_SLIP) // add particle but set it as boundary
				{
					// ... add a particle ...
					vd.add();
					vd.template getLastProp<type>() = BOUNDARY;
				}
				else if (BC_TYPE == NEW_NO_SLIP) // not add particle because already added
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
			}
		}
		else // skip fluid particle
		{
			++fluid_it;
			continue;
		}
		// Set properties
		vd.template getLastProp<rho>() = rho_zero;
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
		++fluid_it;
	}

	// // Now place solid walls
	// openfpm::vector<Box<DIM, double>> holes;

	// if (BC_TYPE == NEW_NO_SLIP)
	// {
	// 	holes.add(recipient_hole);
	// 	sz[0] = sz_aux[0];
	// 	sz[1] = sz_aux[1];
	// }
	// else if (BC_TYPE == NO_SLIP)
	// 	holes.add(fluid_box);

	// Box<DIM, double> hole_get = holes.get(0);
	// auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);

	// if (bc[0] != PERIODIC || bc[1] != PERIODIC) // no walls in all periodic scenario
	// {
	// 	while (bound_box.isNext())
	// 	{
	// 		Point<DIM, double> position = bound_box.get();

	// 		// periodic bc, with no boundary particles in y direction has a bug, it puts 3 extra particles outside in the y direction
	// 		// When running on multiple cores, with this we check if particle is outside the recipient box
	// 		// Another bug places boundary particles in the correct plane, but inside the fluid box;
	// 		// if (bc[0] == PERIODIC && position.get(0) > dp / 2.0 && position.get(0) < length[0] - dp / 2.0)
	// 		// {
	// 		// 	++bound_box;
	// 		// 	continue;
	// 		// }

	// 		if (!recipient.isInside((position)))
	// 		{
	// 			++bound_box;
	// 			continue;
	// 		}
	// 		if (hole_get.isInside(position))
	// 		{
	// 			++bound_box;
	// 			continue;
	// 		}

	// 		if (BC_TYPE == NEW_NO_SLIP && (bc[0] == NON_PERIODIC && bc[1] == NON_PERIODIC))
	// 		{
	// 			// Check if x and z coordinates are multiples of dp, keep multiples, discard the rest
	// 			double remx = fmod(position.get(0), dp);
	// 			double remz = fmod(position.get(1), dp);
	// 			double tol = 0.5 * dp * 10e-2;

	// 			if (remx > tol && remx < dp - tol)
	// 			{
	// 				++bound_box;
	// 				continue;
	// 			}
	// 			if (remz > tol && remz < dp - tol)
	// 			{
	// 				++bound_box;
	// 				continue;
	// 			}
	// 		}
	// 		vd.add();

	// 		vd.template getLastProp<type>() = BOUNDARY;
	// 		vd.template getLastProp<rho>() = rho_zero;
	// 		vd.template getLastProp<pressure>() = 0.0;
	// 		vd.template getLastProp<drho>() = 0.0;

	// 		for (int xyz = 0; xyz < DIM; xyz++)
	// 		{
	// 			vd.getLastPos()[xyz] = bound_box.get().get(xyz);
	// 			vd.template getLastProp<force>()[xyz] = 0.0;
	// 			vd.template getLastProp<force_transport>()[xyz] = 0.0;
	// 			vd.template getLastProp<v_transport>()[xyz] = 0.0;
	// 			vd.template getLastProp<normal_vector>()[xyz] = 0.0;
	// 			if (position.get(1) < dp / 4.0) // bottom wall
	// 			{
	// 				vd.template getLastProp<velocity>()[xyz] = vw_bottom.get(xyz);
	// 			}
	// 			else if (position.get(1) > length[1] - dp / 4.0) // top wall
	// 			{
	// 				vd.template getLastProp<velocity>()[xyz] = vw_top.get(xyz);
	// 			}
	// 		}

	// 		vd.template getLastProp<curvature_boundary>() = 0.0;
	// 		vd.template getLastProp<arc_length>() = dp;

	// 		++bound_box;
	// 	}
	// }
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
										"arc_length"});
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
	LengthScale = StepHeight;

	Nfluid_big[0] = 275;
	Nfluid_big[1] = 31;

	Nfluid_small[0] = 100;
	Nfluid_small[1] = 16;

	size_t Nbound = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
	Nboundary_big[0] = 0;
	Nboundary_big[1] = Nbound;
	double Nboundary_small_up = Nbound;
	double Nboundary_small_down = Nbound + Nfluid_big[1] - Nfluid_small[1];

	bc[0] = PERIODIC;
	bc[1] = NON_PERIODIC;
	dp = LengthScale / ((double)Nfluid_big[1] - Nfluid_small[1]);
	rho_zero = 1.0;
	nu = 1.456 * 1e-2;
	Bfactor = 5.0;
	gravity_vector.get(0) = 0.000603077 / 2.0;
	umax = 1.4 * 1e-1;
	t_end = 2000.0;
	write_const = 1;

	H = Hconst * dp;
	// r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
	r_threshold = 3.0 * H;
	Kcubic = (DIM == 3) ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;
	Kquintic = (DIM == 3) ? 1.0 / 120.0 / M_PI / H / H / H : 7.0 / 478.0 / M_PI / H / H;
	MassFluid = rho_zero * (DIM == 3 ? dp * dp * dp : dp * dp);
	MassBound = rho_zero * (DIM == 3 ? dp * dp * dp : dp * dp);
	cbar = coeff_sound * umax;
	B = rho_zero * cbar * cbar / gamma_;
	Pbackground = Bfactor * B;
	eta = nu * rho_zero;
	Re = umax * 2.0 * 5.2 / nu;

	gravity = getVectorNorm(gravity_vector);

	for (int dim = 0; dim < DIM; dim++)
	{
		if (bc[dim] == NON_PERIODIC) // non periodic, fluid covered by boundary
		{
			length[dim] = dp * (Nfluid_big[dim]);
			length_small[dim] = dp * (Nfluid_small[dim]);
			length_big[dim] = dp * (Nfluid_big[dim]);
			sz[dim] = Nfluid_big[dim] + 2 * (Nboundary_big[dim] + 1);
			offset_domain[dim] = (0.5 + Nboundary_big[dim]) * dp;

			if (Nboundary_big[dim] != 0)
				sz_aux[dim] = 2 * Nfluid_big[dim] - 1 + 2 * (2 * Nboundary_big[dim] + 1 + 1);
			else // for a direction with no boundary particles we dont need to add anything
				sz_aux[dim] = sz[dim];

			if (BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
				offset_recipient[dim] = 0.25 * Nboundary_big[dim] * dp;
			else if (BC_TYPE == NO_SLIP)
				offset_recipient[dim] = Nboundary_big[dim] * dp;
		}
		else // periodic, open ended
		{
			Nfluid_big[dim] -= 1;
			length[dim] = dp * (Nfluid_big[dim] + Nfluid_small[dim]);
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
							{length[0] + offset_domain[0],
							 length[1] + offset_domain[1]});

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
							   {length[0] + offset_recipient[0] + offset_periodic_recipient[0],
								length[1] + offset_recipient[1] + offset_periodic_recipient[1]});

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
	Ghost<DIM, double> g(r_threshold);

	// create particle object
	particles vd_loc(0, domain, bc, g, DEC_GRAN(512));
	vd = vd_loc;

	// correct the number of particles in case of periodicity, we substracted 1 before to accomodate the periodic boundary
	for (int dim = 0; dim < DIM; dim++)
	{
		if (bc[dim] == PERIODIC)
		{
			Nfluid_big[dim] += 1;
			length[dim] += dp;
		}
	}

	// Write constants on file
	std::string customString = "";
	WriteParameters(Nfluid_big, v_cl, customString);

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
	if (PROBES_ENABLED)
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
			probe_filenames.push_back("probes_" + std::to_string(k) + "_" + filename);
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
		vd.template getLastProp<rho>() = rho_zero;
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
		vd.template getLastProp<rho>() = rho_zero;
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

	if (BC_TYPE == NEW_NO_SLIP)
	{
		holes.add(recipient_hole_small);
		holes.add(recipient_hole_big);
		holes.add(CornerHole_New);
		sz[0] = sz_aux[0];
		sz[1] = sz_aux[1];
	}
	else if (BC_TYPE == NO_SLIP)
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

			if (BC_TYPE == NEW_NO_SLIP)
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
			vd.template getLastProp<rho>() = rho_zero;
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
					vd.template getLastProp<velocity>()[xyz] = vw_bottom.get(xyz);
				}
				else if (position.get(1) > length[1] - dp / 4.0) // top wall
				{
					vd.template getLastProp<velocity>()[xyz] = vw_top.get(xyz);
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
										"arc_length"});
	vd.setPropNames(names);
}
int main(int argc, char *argv[])
{

	// initialize the library
	openfpm_init(&argc, &argv);

	// create a Vcluster object ( like MPI communicator )
	Vcluster<> &v_cl = create_vcluster();

	// Create a particle vector
	particles vd;
	std::vector<std::pair<probe_particles, int>> vp;

	Obstacle *obstacle_ptr = nullptr;
	if (SCENARIO == STEP)
	{
		CreateParticleGeometryStep(vd, vp, v_cl);
	}
	else if (SCENARIO == TAYLOR_COUETTE)
	{
		CreateParticleGeometryTaylorCouette(vd, vp, v_cl, obstacle_ptr);
	}
	else
	{
		CreateParticleGeometry(vd, vp, v_cl, obstacle_ptr);
	}
	vd.map();
	vd.write_frame(filename, 0, WRITER);

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
	vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();

	auto NN = vd.getCellList(r_threshold);
	// vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();
	// vd.updateCellList(NN);

	if (BC_TYPE == NO_SLIP) // set up boundary particle velocity
	{
		calc_boundary(vd, NN);
		vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();
	}
	else if (BC_TYPE == NEW_NO_SLIP) // Set up fluid vector and normal vector of the boundary particles
	{
		calcFluidVec(vd, NN);
		vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();
		calcNormalVec(vd, NN);
		vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();
		calcCurvature(vd, NN, v_cl);
		vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();
	}

	// Evolve
	size_t write = 0;
	size_t it = 0;
	size_t it_reb = 0;
	double t = 0.0;
	std::ofstream avgvelstream("avgvel.csv");
	bool calc_drag = false;
	double cylinder_force = 0.0;

	while (t <= t_end)
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
		if (write < (t + dt) * write_const)
		{
			calc_drag = true;
		}

		// Integrate one time step
		kick_drift_int(vd, NN, dt, v_cl, cylinder_force, calc_drag);

		// increment time
		t += dt;
		if (write < t * write_const)
		{
			// // sensor calculation require ghost and update cell-list
			if (PROBES_ENABLED)
			{
				vd.map();
				vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();
				vd.updateCellList(NN);
				for (int k = 0; k < vp.size(); k++)
				{
					probe_particles &probe = vp[k].first;
					// probe.map();
					int &component = vp[k].second;
					sensor_velocity_comp(vd, probe, v_cl, NN, component, obstacle_ptr);
					probe.write_frame(probe_filenames[k], write, WRITER);
				}
			}

			vd.deleteGhost();
			vd.write_frame(filename, write, WRITER);
			vd.ghost_get<type, rho, pressure, velocity, v_transport, normal_vector, curvature_boundary, arc_length>();
			write++;
			computeAverageVelocity(vd, v_cl, t, avgvelstream, cylinder_force);
			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "TIME: " << t << "  write " << cnt << std::endl;
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

	openfpm_finalize();
}
