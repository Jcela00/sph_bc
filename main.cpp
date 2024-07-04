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
#define CYLINDER_ARRAY_LONG 5
#define SQUARE 6

// TYPE OF KERNERL
#define CUBIC 0
#define QUINTIC 1

// TYPE OF DENSITY COMPUTATION
#define DENSITY_SUMMATION 0
#define DENSITY_DIFFERENTIAL 1

//////// CODE PARAMETERS /////////////////////////////////////////////////////////////
// Dimensionality of the problem changes the normalizaiton of the kernel
const int dimensions = 2;

// BC and viscosity we are using in the simulation
const int BC_TYPE = NEW_NO_SLIP;
const int SCENARIO = SQUARE;
const int KERNEL = QUINTIC;
const int DENSITY_TYPE = DENSITY_SUMMATION;
const int WRITER = VTK_WRITER; // VTK_WRITER or CSV_WRITER

// DECLARATION OF GLOBAL VARIABLES
// Output file name
std::string filename;

// Initial particle spacing
double dp;
// Factor relating H (smoothing length) and dp (particle spacing)
const double Hconst = 1.0;
// Smoothing length
double H;
// Size of the kernel support
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

// Gravity vector and magnitude
Point<3, double> gravity_vector;
double gravity;
// Wall velocity
Point<3, double> vw_top;
Point<3, double> vw_bottom;

// Mass of the fluid and boundary particles M=rho*V, V=dp^3 or V = dp^2
double MassFluid;
double MassBound;

// Kinematic viscosity
double nu;
// Dynamic viscosity
double eta;

// Factor relating Pbackground and B, Pbackground = Bfactor*B
const double Bfactor = 1.0;
// Background pressure in the transport force
double Pbackground;

// End simulation time
const double t_end = 100.0;
// Constant used to define time integration
const double CFLnumber = 0.2;
// Minimum T
const double DtMin = 0.00001;
// Controls otput file frequency, low means less frequent
const int write_const = 10;

//////// ALIAS FOR THE PARTICLE PROPERTIES //////////////////////////////////////////
// Properties
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

typedef vector_dist<3, double, aggregate<size_t, double, double, double, double[3], double[3], double[3], double[3], double[3], double, double>> particles;
//                                       |         |     |        |        |        |           |		     |             |		|       |
//                                       |         |     |        |        |        |           |		     |			   |		|       |
//                                    type        rho   pressure delta   force     velocity   pb force   v_transport    normal	curvature   arc_length
//                                                              density
// TEMPLATE TO COPY PASTE AND	 SYNC ALL
// vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();
struct ModelCustom
{
	template <typename Decomposition, typename vector>
	inline void addComputation(Decomposition &dec,
							   vector &vd,
							   size_t v,
							   size_t p)
	{
		if (vd.template getProp<type>(p) == FLUID)
			dec.addComputationCost(v, 4);
		else
			dec.addComputationCost(v, 3);
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
inline void EqState(particles &vd)
{
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto a = it.get();

		double rho_a = vd.template getProp<rho>(a);
		double rho_frac = rho_a / rho_zero;

		vd.template getProp<pressure>(a) = B * (std::pow(rho_frac, gamma_) - 1.0) + xi;

		++it;
	}
}

// Inverted equation of state, compute density given pressure, particle wise
inline double InvEqState_particle(const double pressure)
{
	return rho_zero * std::pow(((pressure - xi) / B + 1.0), 1.0 / gamma_);
}
inline double EqState_particle(const double rho)
{
	return B * (std::pow(rho / rho_zero, gamma_) - 1.0) + xi;
}

// Vector utilities
inline double dotProduct(const Point<3, double> &v, const Point<3, double> &w)
{
	return v.get(0) * w.get(0) + v.get(1) * w.get(1) + v.get(2) * w.get(2);
}
inline double getVectorNorm(const Point<3, double> &v)
{
	return sqrt(v.get(0) * v.get(0) + v.get(1) * v.get(1) + v.get(2) * v.get(2));
}
inline void normalizeVector(Point<3, double> &v)
{
	const double norm = getVectorNorm(v);
	v = v / norm;
}
Point<3, double> getPerpendicularUnit(const Point<3, double> &v)
{
	Point<3, double> perp;
	perp.get(0) = -v.get(1);
	perp.get(1) = v.get(0);
	perp.get(2) = v.get(2);
	normalizeVector(perp);
	return perp;
}
// Kernel functions
inline double Cubic_W(double r)
{
	r /= H;
	if (r >= 0.0 && r < 1.0)
		return Kcubic * (1.0 - 1.5 * r * r + 0.75 * r * r * r);
	else if (r >= 1.0 && r < 2.0)
		return Kcubic * (0.25 * (2.0 - r) * (2.0 - r) * (2.0 - r));
	else
		return 0.0;
}
inline double Quintic_W(double r)
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
inline double Wab(double r)
{
	if (KERNEL == CUBIC)
		return Cubic_W(r);
	else if (KERNEL == QUINTIC)
		return Quintic_W(r);
}

inline Point<3, double> Grad_Cubic_W(const Point<3, double> &dx, const double r)
{
	Point<3, double> DW;
	const double q = r / H;

	const double c1 = Kcubic * (-3.0 / H);
	const double d1 = Kcubic * (9.0 / 4.0 / H);
	const double c2 = Kcubic * (-3.0 / 4.0 / H);

	double factor;
	if (q > 0.0 && q < 1.0)
		factor = c1 * q + d1 * q * q / r;
	else if (q >= 1.0 && q < 2.0)
		factor = c2 * (2.0 - q) * (2.0 - q) / r;
	else
		factor = 0.0;

	DW.get(0) = factor * dx.get(0);
	DW.get(1) = factor * dx.get(1);
	DW.get(2) = factor * dx.get(2);

	return DW;
}

inline Point<3, double> Grad_Quintic_W(const Point<3, double> &dx, const double r)
{
	Point<3, double> DW;
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

	DW.get(0) = factor * dx.get(0);
	DW.get(1) = factor * dx.get(1);
	DW.get(2) = factor * dx.get(2);

	return DW;
}
inline Point<3, double> DWab(const Point<3, double> &dx, const double r)
{
	if (KERNEL == CUBIC)
		return Grad_Cubic_W(dx, r);
	else if (KERNEL == QUINTIC)
		return Grad_Quintic_W(dx, r);
}

inline Point<3, double> Pi_physical(const Point<3, double> &dr, const double &r, const Point<3, double> &dv, const Point<3, double> &dW)
{
	return eta * (dv * dotProduct(dr, dW)) / r / r; // / 1.44;
}

inline double PressureForce(const double &rhoa, const double &rhob, const double &prsa, const double &prsb)
{
	return -1.0 * (rhob * prsa + rhoa * prsb) / (rhoa + rhob); // / 1.44;
}

std::array<Point<3, double>, 3> GetBoundaryPositions(const Point<3, double> &r, const Point<3, double> &normal)
{
	// aplies offset to a vector and returns the vectors pointing at the dummy wall particles

	Point<3, double> offset_1 = -0.5 * normal * dp; // offset from wall to first particle
	Point<3, double> offset_2 = -1.0 * normal * dp; // offset from first to second particle, and from second to third
	Point<3, double> r1 = r + offset_1;				// first dummy particle
	Point<3, double> r2 = r1 + offset_2;			// second dummy particle
	Point<3, double> r3 = r2 + offset_2;			// third dummy particle

	std::array<Point<3, double>, 3> r_dummy = {r1, r2, r3};

	return r_dummy;
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
			Point<3, double> xa = vd.getPos(a);

			// initialize sum
			Point<3, double> r_fluid_sum = {0.0, 0.0, 0.0};
			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

			// iterate the neighborhood particles
			while (Np.isNext() == true)
			{
				// Key of b particle
				const unsigned long b = Np.get();

				if (vd.getProp<type>(b) == FLUID)
				{
					// Get the position xb of the particle b
					const Point<3, double> xb = vd.getPos(b);

					// Get the vector pointing at fluid from boundary
					Point<3, double> dr = xb - xa;

					// take the norm squared of this vector
					const double r2 = norm2(dr);

					// If the particles interact ...
					if (r2 < r_threshold * r_threshold)
					{
						normalizeVector(dr);
						double W = Wab(sqrt(r2));
						r_fluid_sum += dr * W;
					}
				}

				++Np;
			}
			normalizeVector(r_fluid_sum);
			// store in force because we are not using it for boundary particles
			vd.template getProp<force>(a)[0] = r_fluid_sum.get(0);
			vd.template getProp<force>(a)[1] = r_fluid_sum.get(1);
			vd.template getProp<force>(a)[2] = r_fluid_sum.get(2);
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
			Point<3, double> xa = vd.getPos(a);

			// get vector that points at fluid
			Point<3, double> Rfluid = vd.getProp<force>(a);

			// initialize sum
			Point<3, double> n_sum = {0.0, 0.0, 0.0};
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
					const Point<3, double> xb = vd.getPos(b);

					// Get the vector pointing at b from a
					const Point<3, double> dr = xa - xb;

					// take the norm squared of this vector
					const double r2 = norm2(dr);

					// If the particles interact ...
					if (r2 < r_threshold * r_threshold)
					{
						// get perpendicular vector to dr
						Point<3, double> perp = getPerpendicularUnit(dr); // this is normalized

						// get scalar product of perp and Rfluid
						double perp_dot_fluid = dotProduct(perp, Rfluid);

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
			normalizeVector(n_sum);
			// store in normal vector
			vd.template getProp<normal_vector>(a)[0] = n_sum.get(0);
			vd.template getProp<normal_vector>(a)[1] = n_sum.get(1);
			vd.template getProp<normal_vector>(a)[2] = n_sum.get(2);
		}

		++part;
	}
}

template <typename CellList>
void calcCurvature(particles &vd, CellList &NN)
{
	// This function computes the curvature of the boundary particles from the divergence of the normal vector

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
			Point<3, double> xa = vd.getPos(a);

			// get normal of a
			Point<3, double> normal_a = vd.getProp<normal_vector>(a);

			// initialize sums
			double K_sum = 0.0;
			double w_sum = 0.0;

			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

			// iterate the neighborhood particles
			while (Np.isNext() == true)
			{
				// Key of b particle
				const unsigned long b = Np.get();

				// if (a == b) skip this particle
				// if (a.getKey() == b)
				// {
				// 	++Np;
				// 	continue;
				// }

				if (vd.getProp<type>(b) == BOUNDARY)
				{
					// Get the position xb of the particle b
					const Point<3, double> xb = vd.getPos(b);

					// Get the vector pointing at a from b
					const Point<3, double> dr = xa - xb;

					// take the norm squared of this vector
					const double r2 = norm2(dr);

					// If the particles interact ...
					if (r2 < r_threshold * r_threshold)
					{
						// get normal vector of b
						Point<3, double> normal_b = vd.getProp<normal_vector>(b);

						// evaluate kernel
						double W = Wab(sqrt(r2));
						// evaluate kernel gradient
						Point<3, double> dW = DWab(dr, sqrt(r2));

						K_sum += dotProduct(normal_b - normal_a, dW);
						w_sum += W;
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
			Point<3, double> xa = vd.getPos(a);

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
				const Point<3, double> xb = vd.getPos(b);

				// Get the vector pointing at xa from xb
				const Point<3, double> dr = xa - xb;

				// take the norm squared of this vector
				const double r2 = norm2(dr);

				// If the particles interact ...
				if (r2 < r_threshold * r_threshold)
				{
					if (vd.getProp<type>(b) == FLUID)
					{
						// calculate distance
						const double r = sqrt(r2);

						// evaluate kernel
						const double w = Wab(r);

						W_sum += w;
					}
					else
					{
						if (BC_TYPE == NO_SLIP)
						{
							// calculate distance
							const double r = sqrt(r2);

							// evaluate kernel
							const double w = Wab(r);

							W_sum += w;
						}
						else if (BC_TYPE == NEW_NO_SLIP) // need to evaluate kernel at dummy particles
						{
							// xw.get(0) this is the x coordinate of the wall
							Point<3, double> normal = vd.getProp<normal_vector>(b);

							// Apply offsets to dr to get 3 vectrors pointing to dummy particles
							std::array<Point<3, double>, 3> R_dummy = GetBoundaryPositions(-1.0 * dr, normal);

							const Point<3, double> r1 = R_dummy[0];
							const Point<3, double> r2 = R_dummy[1];
							const Point<3, double> r3 = R_dummy[2];

							const double W1 = Wab(sqrt(norm2(r1)));
							const double W2 = Wab(sqrt(norm2(r2)));
							const double W3 = Wab(sqrt(norm2(r3)));

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
			Point<3, double> xa = vd.getPos(a);

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
				const Point<3, double> xb = vd.getPos(b);

				// Get the vector pointing at xa from xb
				const Point<3, double> dr = xa - xb;

				// take the norm squared of this vector
				const double r2 = norm2(dr);

				// If the particles interact ...
				if (r2 < r_threshold * r_threshold)
				{

					if (vd.getProp<type>(b) == FLUID)
					{
						// calculate distance
						const double r = sqrt(r2);

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
							// xw.get(0) this is the x coordinate of the wall
							const Point<3, double> normal = vd.getProp<normal_vector>(b);

							// Apply offsets to dr to get 3 vectrors pointing to dummy particles
							const std::array<Point<3, double>, 3> R_dummy = GetBoundaryPositions(-1.0 * dr, normal);
							const double kappa = vd.getProp<curvature_boundary>(b);
							const double dxwall = vd.getProp<arc_length>(b);

							const double Vol1 = 0.5 * (2.0 * dp + dp * dp * kappa - 2.0 * 1.0 * dp * dp * kappa) * dxwall;
							const double Vol2 = 0.5 * (2.0 * dp + dp * dp * kappa - 2.0 * 2.0 * dp * dp * kappa) * dxwall;
							const double Vol3 = 0.5 * (2.0 * dp + dp * dp * kappa - 2.0 * 3.0 * dp * dp * kappa) * dxwall;

							const double mass1 = Vol1 * rho_zero;
							const double mass2 = Vol2 * rho_zero;
							const double mass3 = Vol3 * rho_zero;

							const Point<3, double> r1 = R_dummy[0];
							const Point<3, double> r2 = R_dummy[1];
							const Point<3, double> r3 = R_dummy[2];

							const double W1 = Wab(getVectorNorm(r1));
							const double W2 = Wab(getVectorNorm(r2));
							const double W3 = Wab(getVectorNorm(r3));

							rho_sum += (W1 * mass1 + W2 * mass2 + W3 * mass3);
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
			Point<3, double> xb = vd.getPos(b);

			// Get the Velocity of the boundary particle
			Point<3, double> vb = vd.getProp<velocity>(b);

			// initialize sums
			Point<3, double> sum_vW;
			sum_vW.get(0) = 0.0;
			sum_vW.get(1) = 0.0;
			sum_vW.get(2) = 0.0;
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
				Point<3, double> xf = vd.getPos(f);

				// Get the velocity of the fluid particle
				Point<3, double> vf = vd.getProp<velocity>(f);

				// Get the density of the fluid particle
				double rhof = vd.getProp<rho>(f);

				// Get the pressure of the fluid particle
				double Pf = vd.getProp<pressure>(f);

				// Get the vector pointing at xb from xf rwf
				Point<3, double> dr = xb - xf;

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
					sum_vW.get(0) += w * vf.get(0);
					sum_vW.get(1) += w * vf.get(1);
					sum_vW.get(2) += w * vf.get(2);

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
				vd.template getProp<v_transport>(b)[0] = 2.0 * vd.template getProp<velocity>(b)[0] - sum_vW.get(0) / sum_W;
				vd.template getProp<v_transport>(b)[1] = 2.0 * vd.template getProp<velocity>(b)[1] - sum_vW.get(1) / sum_W;
				vd.template getProp<v_transport>(b)[2] = 2.0 * vd.template getProp<velocity>(b)[2] - sum_vW.get(2) / sum_W;
				// Set the pressure of the boundary particle b
				vd.template getProp<pressure>(b) = sum_pW / sum_W;
				// Compute density from inverted Eq of state
				vd.template getProp<rho>(b) = InvEqState_particle(vd.template getProp<pressure>(b));
			}
			else
			{
				vd.template getProp<v_transport>(b)[0] = 2.0 * vd.template getProp<velocity>(b)[0];
				vd.template getProp<v_transport>(b)[1] = 2.0 * vd.template getProp<velocity>(b)[1];
				vd.template getProp<v_transport>(b)[2] = 2.0 * vd.template getProp<velocity>(b)[2];
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
								 const Point<3, double> &xf,
								 const Point<3, double> &vf,
								 const Point<3, double> &xw,
								 const Point<3, double> &r_wall_to_fluid,
								 const double &dist2,
								 unsigned long &boundary_key)
{

	// Points from fluid to wall
	Point<3, double> r_fluid_to_wall = -1.0 * r_wall_to_fluid;

	const double massw = MassBound;
	const double rhow = vd.getProp<rho>(boundary_key);
	const double Pw = vd.getProp<pressure>(boundary_key);
	const Point<3, double> vw = vd.getProp<velocity>(boundary_key);
	const Point<3, double> vtf = vd.getProp<v_transport>(fluid_key);
	const Point<3, double> vdiff_f = vtf - vf;

	// Get normal vector
	Point<3, double> normal = vd.getProp<normal_vector>(boundary_key);
	Point<3, double> tangential = {normal.get(1), normal.get(0), 0.0};
	if (tangential.get(1) < 0.0)
	{
		tangential.get(1) = -tangential.get(1);
	}

	// get curvature, and arc length
	double kappa = vd.getProp<curvature_boundary>(boundary_key);
	double dxwall = vd.getProp<arc_length>(boundary_key);

	// Project va on tangential and normal directions
	double vt = dotProduct(vf, tangential);
	double vn = dotProduct(vf, normal);

	// vertical distance from fluid particle to wall
	double lf = dotProduct(r_fluid_to_wall, normal);
	lf = (lf < 0.0 ? -1.0 * lf : lf); // absolute value

	// Get array of vectors from fluid to boundary particles, get norm, and get normal distance to wall
	std::array<Point<3, double>, 3> r_boundary = GetBoundaryPositions(r_fluid_to_wall, normal);
	std::array<double, 3> r_boundary_norm = {getVectorNorm(r_boundary[0]), getVectorNorm(r_boundary[1]), getVectorNorm(r_boundary[2])};
	std::array<double, 3> lwall = {0.5 * dp, 1.5 * dp, 2.5 * dp};

	// Initialize arrays for boundary particles
	std::array<double, 3> p_boundary = {0.0, 0.0, 0.0};
	std::array<double, 3> rho_boundary = {0.0, 0.0, 0.0};
	std::array<Point<3, double>, 3> v_boundary = {Point<3, double>{0.0, 0.0, 0.0}, Point<3, double>{0.0, 0.0, 0.0}, Point<3, double>{0.0, 0.0, 0.0}};
	std::array<double, 3> Volume_boundary = {0.0, 0.0, 0.0};
	std::array<double, 3> Mass_boundary = {0.0, 0.0, 0.0};
	double g_normal = dotProduct(gravity_vector, normal);

	size_t interact_count = 0;
	const double GradATerm = 0.5 * (rhof * dotProduct(vf, vdiff_f));

	for (int i = 0; i < 3; i++)
	{
		if (r_boundary_norm[i] < r_threshold)
		{

			// compute volume of boundary particle, this gives density and pressure
			Volume_boundary[i] = 0.5 * (2.0 * dp + dp * dp * kappa - 2.0 * (i + 1.0) * dp * dp * kappa) * dxwall;
			Mass_boundary[i] = Volume_boundary[i] * rho_zero;

			// rho_boundary[i] = MassBound / Volume_boundary[i];
			// p_boundary[i] = EqState_particle(rho_boundary[i]);

			// compute velocity of boundary particle
			v_boundary[i] = 2.0 * vw - vt * (lwall[i] / lf) * tangential - vn * (lwall[i] / lf) * normal;
			p_boundary[i] = Pf + rhof * g_normal * (lf + lwall[i]);
			rho_boundary[i] = InvEqState_particle(p_boundary[i]);

			// flip sing of r_boundary to get vector pointing from boundary to fluid (Force routines use the vector pointing from b to a)
			r_boundary[i] = -1.0 * r_boundary[i];

			// Evaluate kernel gradient
			const Point<3, double> DW = DWab(r_boundary[i], r_boundary_norm[i]);

			// Compute forces
			const Point<3, double> v_rel = vf - v_boundary[i];
			const double Va2 = (massf / rhof) * (massf / rhof);
			const double Vb2 = (Mass_boundary[i] / rho_boundary[i]) * (Mass_boundary[i] / rho_boundary[i]);

			// const double Vb2 = Volume_boundary[i] * Volume_boundary[i];

			const Point<3, double> ViscosityTerm = Pi_physical(r_boundary[i], getVectorNorm(r_boundary[i]), v_rel, DW);
			const double PressureTerm = PressureForce(rhof, rho_boundary[i], Pf, p_boundary[i]);
			vd.getProp<force>(fluid_key)[0] += (Va2 + Vb2) * ((GradATerm + PressureTerm) * DW.get(0) + ViscosityTerm.get(0)) / massf;
			vd.getProp<force>(fluid_key)[1] += (Va2 + Vb2) * ((GradATerm + PressureTerm) * DW.get(1) + ViscosityTerm.get(1)) / massf;
			vd.getProp<force>(fluid_key)[2] += (Va2 + Vb2) * ((GradATerm + PressureTerm) * DW.get(2) + ViscosityTerm.get(2)) / massf;

			vd.getProp<force_transport>(fluid_key)[0] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(0) / massf;
			vd.getProp<force_transport>(fluid_key)[1] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(1) / massf;
			vd.getProp<force_transport>(fluid_key)[2] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(2) / massf;

			if (DENSITY_TYPE == DENSITY_DIFFERENTIAL)
			{
				// vd.getProp<drho>(fluid_key) += MassBound * dotProduct(v_rel, DW);
				vd.getProp<drho>(fluid_key) += Mass_boundary[i] * dotProduct(v_rel, DW);
			}
		}
	}
}

void interact_fluid_fluid(particles &vd,
						  const vect_dist_key_dx &a_key,
						  const double &massa,
						  const double &rhoa,
						  const double &Pa,
						  const Point<3, double> &xa,
						  const Point<3, double> &va,
						  const Point<3, double> &xb,
						  const Point<3, double> &r_ab,
						  const double &r2,
						  const unsigned long &b_key)
{

	const double massb = MassFluid;
	const double rhob = vd.getProp<rho>(b_key);
	const double Pb = vd.getProp<pressure>(b_key);
	const Point<3, double> vb = vd.getProp<velocity>(b_key);

	const double r = sqrt(r2);

	const Point<3, double> v_rel = va - vb;

	const Point<3, double> DW = DWab(r_ab, r);

	const double Va2 = (massa / rhoa) * (massa / rhoa);
	const double Vb2 = (massb / rhob) * (massb / rhob);

	const Point<3, double> ViscosityTerm = Pi_physical(r_ab, r, v_rel, DW);
	const double PressureTerm = PressureForce(rhoa, rhob, Pa, Pb);

	const Point<3, double> vta = vd.getProp<v_transport>(a_key);
	const Point<3, double> vtb = vd.getProp<v_transport>(b_key);
	const Point<3, double> vdiff_a = vta - va;
	const Point<3, double> vdiff_b = vtb - vb;

	double GradATerm = 0.5 * (rhoa * dotProduct(va, vdiff_a) + rhob * dotProduct(vb, vdiff_b));

	vd.getProp<force>(a_key)[0] += (Va2 + Vb2) * ((GradATerm + PressureTerm) * DW.get(0) + ViscosityTerm.get(0)) / massa;
	vd.getProp<force>(a_key)[1] += (Va2 + Vb2) * ((GradATerm + PressureTerm) * DW.get(1) + ViscosityTerm.get(1)) / massa;
	vd.getProp<force>(a_key)[2] += (Va2 + Vb2) * ((GradATerm + PressureTerm) * DW.get(2) + ViscosityTerm.get(2)) / massa;

	vd.getProp<force_transport>(a_key)[0] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(0) / massa;
	vd.getProp<force_transport>(a_key)[1] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(1) / massa;
	vd.getProp<force_transport>(a_key)[2] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(2) / massa;

	if (DENSITY_TYPE == DENSITY_DIFFERENTIAL)
	{
		vd.getProp<drho>(a_key) += massb * dotProduct(v_rel, DW);
	}
}

void interact_fluid_boundary_old(particles &vd,
								 vect_dist_key_dx &fluid_key,
								 const double &massf,
								 const double &rhof,
								 const double &Pf,
								 const Point<3, double> &xf,
								 const Point<3, double> &vf,
								 const Point<3, double> &xb,
								 const Point<3, double> &r_ab,
								 const double &r2,
								 unsigned long &boundary_key)
{
	const double massb = MassBound;
	const double rhob = vd.getProp<rho>(boundary_key);
	const double Pb = vd.getProp<pressure>(boundary_key);
	const Point<3, double> vb = vd.getProp<velocity>(boundary_key);
	const Point<3, double> vb_noslip = vd.getProp<v_transport>(boundary_key);

	const Point<3, double> r_fb = xf - xb;
	const double r = sqrt(r2);

	const Point<3, double> v_rel = vf - vb;
	const Point<3, double> v_rel_aux = vf - vb_noslip;

	const Point<3, double> DW = DWab(r_ab, r);

	const double Va2 = (massf / rhof) * (massf / rhof);
	const double Vb2 = (massb / rhob) * (massb / rhob);

	const Point<3, double> ViscosityTerm = Pi_physical(r_ab, r, v_rel_aux, DW);
	const double PressureTerm = PressureForce(rhof, rhob, Pf, Pb);

	const Point<3, double> vtf = vd.getProp<v_transport>(fluid_key);
	const Point<3, double> vdiff_f = vtf - vf;
	// Boundary particles have no transport velocity difference

	const double GradATerm = 0.5 * (rhof * dotProduct(vf, vdiff_f));

	vd.getProp<force>(fluid_key)[0] += (Va2 + Vb2) * ((GradATerm + PressureTerm) * DW.get(0) + ViscosityTerm.get(0)) / massf;
	vd.getProp<force>(fluid_key)[1] += (Va2 + Vb2) * ((GradATerm + PressureTerm) * DW.get(1) + ViscosityTerm.get(1)) / massf;
	vd.getProp<force>(fluid_key)[2] += (Va2 + Vb2) * ((GradATerm + PressureTerm) * DW.get(2) + ViscosityTerm.get(2)) / massf;

	vd.getProp<force_transport>(fluid_key)[0] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(0) / massf;
	vd.getProp<force_transport>(fluid_key)[1] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(1) / massf;
	vd.getProp<force_transport>(fluid_key)[2] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(2) / massf;

	if (DENSITY_TYPE == DENSITY_DIFFERENTIAL)
	{
		vd.getProp<drho>(fluid_key) += massb * dotProduct(v_rel, DW);
	}
}

template <typename CellList>
inline void calc_forces(particles &vd, CellList &NN)
{
	vector_dist_iterator part = vd.getDomainIterator();

	// Update the cell-list
	vd.updateCellList(NN);

	// For each particle ...
	while (part.isNext())
	{
		// get key of particle a

		vect_dist_key_dx a = part.get();

		// We threat FLUID particle differently from BOUNDARY PARTICLES
		// Boundary particles dont need any force computation

		if (vd.getProp<type>(a) == FLUID) // INTERACTION OF A FLUID PARTICLE WITH ITS NEIGHBORHOOD
		{

			// Get the position xp of the particle
			Point<3, double> xa = vd.getPos(a);

			// Take the mass of the particle dependently if it is FLUID or BOUNDARY
			double massa = MassFluid;

			// Get the density and pressure of the of the particle a
			double rhoa = vd.getProp<rho>(a);
			double Pa = vd.getProp<pressure>(a);

			// Get the Velocity of the particle a
			Point<3, double> va = vd.getProp<velocity>(a);

			// Reset the force counter (0 + gravity)
			vd.template getProp<force>(a)[0] = gravity_vector.get(0);
			vd.template getProp<force>(a)[1] = gravity_vector.get(1);
			vd.template getProp<force>(a)[2] = gravity_vector.get(2);

			vd.template getProp<force_transport>(a)[0] = 0.0;
			vd.template getProp<force_transport>(a)[1] = 0.0;
			vd.template getProp<force_transport>(a)[2] = 0.0;

			vd.template getProp<drho>(a) = 0.0;

			// Get an iterator over the neighborhood particles of p
			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

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
				Point<3, double> xb = vd.getPos(b);

				// Get the distance between a and b
				// in fluid - boundary its xf-xb
				Point<3, double> dr = xa - xb;

				// take the norm (squared) of this vector
				double r2 = norm2(dr);

				// if they interact
				if (r2 < r_threshold * r_threshold)
				{
					if (vd.getProp<type>(b) == BOUNDARY)
					{
						if (BC_TYPE == NO_SLIP)
							interact_fluid_boundary_old(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b);
						else if (BC_TYPE == NEW_NO_SLIP)
							interact_fluid_boundary_new(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b);
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

void max_acceleration_and_velocity(particles &vd, double &max_acc, double &max_vel)
{
	// Calculate the maximum acceleration
	auto part = vd.getDomainIterator();

	while (part.isNext())
	{
		auto a = part.get();

		Point<3, double> acc(vd.getProp<force>(a));
		double acc2 = norm2(acc);

		Point<3, double> vel(vd.getProp<velocity>(a));
		double vel2 = norm2(vel);

		if (vel2 >= max_vel)
			max_vel = vel2;

		if (acc2 >= max_acc)
			max_acc = acc2;

		++part;
	}
	max_acc = sqrt(max_acc);
	max_vel = sqrt(max_vel);

	Vcluster<> &v_cl = create_vcluster();
	v_cl.max(max_acc);
	v_cl.max(max_vel);
	v_cl.execute();
}

double calc_deltaT(particles &vd)
{
	double Maxacc = 0.0;
	double Maxvel = 0.0;
	max_acceleration_and_velocity(vd, Maxacc, Maxvel);

	double dt_u = 0.25 * H / (cbar + abs(Maxvel));
	double dt_visc = 0.25 * H * H / (nu);
	double dt_g = 0.25 * sqrt(H / getVectorNorm(gravity_vector));
	double dt = CFLnumber * std::min({dt_u, dt_visc, dt_g});

	if (dt < DtMin)
		dt = DtMin;

	return dt;
}

size_t cnt = 0;

template <typename CellList>
void kick_drift_int(particles &vd, CellList &NN, const double dt, Vcluster<> &v_cl)
{
	// particle iterator
	auto part = vd.getDomainIterator();

	const double dt_2 = dt * 0.5;

	// For each particle ...
	while (part.isNext())
	{
		// get particle a key
		vect_dist_key_dx a = part.get();

		// if the particle is boundary skip
		if (vd.template getProp<type>(a) == BOUNDARY)
		{
			++part;
			continue;
		}

		vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + dt_2 * vd.template getProp<force>(a)[0];
		vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + dt_2 * vd.template getProp<force>(a)[1];
		vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + dt_2 * vd.template getProp<force>(a)[2];

		vd.template getProp<v_transport>(a)[0] = vd.template getProp<velocity>(a)[0] + dt_2 * vd.template getProp<force_transport>(a)[0];
		vd.template getProp<v_transport>(a)[1] = vd.template getProp<velocity>(a)[1] + dt_2 * vd.template getProp<force_transport>(a)[1];
		vd.template getProp<v_transport>(a)[2] = vd.template getProp<velocity>(a)[2] + dt_2 * vd.template getProp<force_transport>(a)[2];

		vd.getPos(a)[0] += dt * vd.template getProp<v_transport>(a)[0];
		vd.getPos(a)[1] += dt * vd.template getProp<v_transport>(a)[1];
		vd.getPos(a)[2] += dt * vd.template getProp<v_transport>(a)[2];

		if (DENSITY_TYPE == DENSITY_DIFFERENTIAL)
			vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);

		++part;
	}
	vd.map();
	v_cl.barrier();
	vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();

	if (DENSITY_TYPE == DENSITY_SUMMATION)
		calc_density(vd, NN);

	vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();
	v_cl.barrier();

	// Calculate pressure from the density
	EqState(vd);
	vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();
	v_cl.barrier();

	if (BC_TYPE == NO_SLIP)
	{
		calc_boundary(vd, NN);
		vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();
	}
	v_cl.barrier();
	calc_forces(vd, NN);
	v_cl.barrier();
	vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();

	// particle iterator
	auto part2 = vd.getDomainIterator();

	// For each particle ...
	while (part2.isNext())
	{
		// ... a
		auto a = part2.get();

		// if the particle is boundary skip
		if (vd.template getProp<type>(a) == BOUNDARY)
		{
			++part2;
			continue;
		}

		// if (DENSITY_TYPE == DENSITY_DIFFERENTIAL)
		// 	vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);

		vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + dt_2 * vd.template getProp<force>(a)[0];
		vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + dt_2 * vd.template getProp<force>(a)[1];
		vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + dt_2 * vd.template getProp<force>(a)[2];

		vd.template getProp<v_transport>(a)[0] = vd.template getProp<velocity>(a)[0] + dt_2 * vd.template getProp<force_transport>(a)[0];
		vd.template getProp<v_transport>(a)[1] = vd.template getProp<velocity>(a)[1] + dt_2 * vd.template getProp<force_transport>(a)[1];
		vd.template getProp<v_transport>(a)[2] = vd.template getProp<velocity>(a)[2] + dt_2 * vd.template getProp<force_transport>(a)[2];

		++part2;
	}

	vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();
	v_cl.barrier();

	// increment the iteration counter
	cnt++;
}

void SetFilename(std::string &filename, std::string custom_string = "")
{
	filename = "";
	if (SCENARIO == POISEUILLE)
	{
		filename = "Poiseuille";
	}
	else if (SCENARIO == COUETTE)
	{
		filename = "Couette";
	}
	else if (SCENARIO == HYDROSTATIC)
	{
		filename = "Hydrostatic";
	}
	else if (SCENARIO == CYLINDER_ARRAY)
	{
		filename = "CylinderArray";
	}
	else if (SCENARIO == CYLINDER_LATTICE)
	{
		filename = "CylinderLattice";
	}
	else if (SCENARIO == CYLINDER_ARRAY_LONG)
	{
		filename = "CylinderArrayLong";
	}
	else if (SCENARIO == SQUARE)
	{
		filename = "Square";
	}
	if (BC_TYPE == NO_SLIP)
	{
		filename += "_OLD_BC";
	}
	else if (BC_TYPE == NEW_NO_SLIP)
	{
		filename += "_NEW_BC";
	}
	if (KERNEL == CUBIC)
	{
		filename += "_Cubic";
	}
	else if (KERNEL == QUINTIC)
	{
		filename += "_Quintic";
	}
	if (DENSITY_TYPE == DENSITY_SUMMATION)
	{
		filename += "_Summation";
	}
	else if (DENSITY_TYPE == DENSITY_DIFFERENTIAL)
	{
		filename += "_Differential";
	}

	filename += ("_" + custom_string);
}

void WriteParameters(size_t Nfluid[3], double length[3], std::string customString)
{
	const double Lx = length[0] + 0.5 * dp; // channel height
	const double Ly = length[1];			// channel width
	const double Lz = length[2];			// channel length ( +0.5dp due to periodicity)

	std::string reynolds_size_name = std::to_string(Nfluid[0]) + "_" + std::to_string(Nfluid[1]);
	SetFilename(filename, reynolds_size_name + "_" + customString);

	std::string constants_filename = "Constants_" + filename + ".txt";
	std::ofstream file(constants_filename);

	file << "gravity = {" << gravity_vector.get(0) << ", " << gravity_vector.get(1) << ", " << gravity_vector.get(2) << "}" << std::endl;
	file << "v_topwall = {" << vw_top.get(0) << ", " << vw_top.get(1) << ", " << vw_top.get(2) << "}" << std::endl;
	file << "v_bottomwall = {" << vw_bottom.get(0) << ", " << vw_bottom.get(1) << ", " << vw_bottom.get(2) << "}" << std::endl;
	file << "Lx, Ly, Lz = " << Lx << ", " << Ly << ", " << Lz << std::endl;
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
	file << "SCENARIO = " << SCENARIO << std::endl;
	file << "BC_TYPE = " << BC_TYPE << std::endl;
	file << "KERNEL = " << KERNEL << std::endl;
	file << "DENSITY_TYPE = " << DENSITY_TYPE << std::endl;
}

void AddCylinderNewBC(particles &vd, Point<3, double> CylinderCentre, double CylinderRadius)
{
	const double perimeter = 2.0 * M_PI * CylinderRadius;
	const int Np_cylinder = ceil(perimeter / dp);
	const double dtheta = 2.0 * M_PI / Np_cylinder;
	const double dxwall = dtheta * CylinderRadius;
	double theta = 0.0;

	while (theta < 2.0 * M_PI)
	{
		Point<3, double> cylinder_particle{CylinderCentre.get(0) + CylinderRadius * cos(theta), CylinderCentre.get(1) + CylinderRadius * sin(theta), CylinderCentre.get(2)};
		Point<3, double> normal = {cos(theta), sin(theta), 0.0};
		vd.add();
		vd.getLastPos()[0] = cylinder_particle.get(0);
		vd.getLastPos()[1] = cylinder_particle.get(1);
		vd.getLastPos()[2] = cylinder_particle.get(2);

		vd.template getLastProp<type>() = BOUNDARY;
		vd.template getLastProp<pressure>() = 0.0;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<drho>() = 0.0;

		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;

		vd.template getLastProp<force_transport>()[0] = 0.0;
		vd.template getLastProp<force_transport>()[1] = 0.0;
		vd.template getLastProp<force_transport>()[2] = 0.0;

		vd.template getLastProp<v_transport>()[0] = 0.0;
		vd.template getLastProp<v_transport>()[1] = 0.0;
		vd.template getLastProp<v_transport>()[2] = 0.0;

		// vd.template getLastProp<normal_vector>()[0] = normal.get(0);
		// vd.template getLastProp<normal_vector>()[1] = normal.get(1);
		// vd.template getLastProp<normal_vector>()[2] = normal.get(2);

		vd.template getLastProp<curvature_boundary>() = 0.0; // 1.0 / CylinderRadius;
		vd.template getLastProp<arc_length>() = dxwall;
		theta += dtheta;
	}
}

void AddSquareNewBC(particles &vd, Point<3, double> SquareCentre, int integerSideLength)
{
	const double sideLength = (integerSideLength - 1) * dp;
	const Point<3, double> LowerLeftCorner = {SquareCentre.get(0) - sideLength / 2.0, SquareCentre.get(1) - sideLength / 2.0, SquareCentre.get(2)};
	const Point<3, double> UpperLeftCorner = {SquareCentre.get(0) - sideLength / 2.0, SquareCentre.get(1) + sideLength / 2.0, SquareCentre.get(2)};
	const Point<3, double> LowerRightCorner = {SquareCentre.get(0) + sideLength / 2.0, SquareCentre.get(1) - sideLength / 2.0, SquareCentre.get(2)};

	// upper and lower walls
	for (int k = 0; k < integerSideLength; k++)
	{
		Point<3, double> Xoffset = {k * dp, 0.0, 0.0};
		Point<3, double> LowerWall = LowerLeftCorner + Xoffset;
		Point<3, double> UpperWall = UpperLeftCorner + Xoffset;

		vd.add();
		vd.getLastPos()[0] = LowerWall.get(0); //+ ((double)rand() / RAND_MAX - 0.5) * dp;
		vd.getLastPos()[1] = LowerWall.get(1); //+ ((double)rand() / RAND_MAX - 0.5) * dp;
		vd.getLastPos()[2] = LowerWall.get(2);

		vd.template getLastProp<type>() = BOUNDARY;
		vd.template getLastProp<pressure>() = 0.0;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<drho>() = 0.0;

		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;

		vd.template getLastProp<force_transport>()[0] = 0.0;
		vd.template getLastProp<force_transport>()[1] = 0.0;
		vd.template getLastProp<force_transport>()[2] = 0.0;

		vd.template getLastProp<v_transport>()[0] = 0.0;
		vd.template getLastProp<v_transport>()[1] = 0.0;
		vd.template getLastProp<v_transport>()[2] = 0.0;

		vd.template getLastProp<normal_vector>()[0] = 0.0;
		vd.template getLastProp<normal_vector>()[1] = 0.0;
		vd.template getLastProp<normal_vector>()[2] = 0.0;

		vd.template getLastProp<curvature_boundary>() = 0.0;
		vd.template getLastProp<arc_length>() = dp;

		vd.add();
		vd.getLastPos()[0] = UpperWall.get(0);
		vd.getLastPos()[1] = UpperWall.get(1);
		vd.getLastPos()[2] = UpperWall.get(2);

		vd.template getLastProp<type>() = BOUNDARY;
		vd.template getLastProp<pressure>() = 0.0;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<drho>() = 0.0;

		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;

		vd.template getLastProp<force_transport>()[0] = 0.0;
		vd.template getLastProp<force_transport>()[1] = 0.0;
		vd.template getLastProp<force_transport>()[2] = 0.0;

		vd.template getLastProp<v_transport>()[0] = 0.0;
		vd.template getLastProp<v_transport>()[1] = 0.0;
		vd.template getLastProp<v_transport>()[2] = 0.0;

		vd.template getLastProp<normal_vector>()[0] = 0.0;
		vd.template getLastProp<normal_vector>()[1] = 0.0;
		vd.template getLastProp<normal_vector>()[2] = 0.0;

		vd.template getLastProp<curvature_boundary>() = 0.0;
		vd.template getLastProp<arc_length>() = dp;
	}
	// left and right walls

	for (int k = 1; k < integerSideLength - 1; k++)
	{
		Point<3, double> Yoffset = {0.0, k * dp, 0.0};
		Point<3, double> LeftWall = LowerLeftCorner + Yoffset;
		Point<3, double> RightWall = LowerRightCorner + Yoffset;

		vd.add();
		vd.getLastPos()[0] = LeftWall.get(0);
		vd.getLastPos()[1] = LeftWall.get(1);
		vd.getLastPos()[2] = LeftWall.get(2);

		vd.template getLastProp<type>() = BOUNDARY;
		vd.template getLastProp<pressure>() = 0.0;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<drho>() = 0.0;

		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;

		vd.template getLastProp<force_transport>()[0] = 0.0;
		vd.template getLastProp<force_transport>()[1] = 0.0;
		vd.template getLastProp<force_transport>()[2] = 0.0;

		vd.template getLastProp<v_transport>()[0] = 0.0;
		vd.template getLastProp<v_transport>()[1] = 0.0;
		vd.template getLastProp<v_transport>()[2] = 0.0;

		vd.template getLastProp<normal_vector>()[0] = 0.0;
		vd.template getLastProp<normal_vector>()[1] = 0.0;
		vd.template getLastProp<normal_vector>()[2] = 0.0;

		vd.template getLastProp<curvature_boundary>() = 0.0;
		vd.template getLastProp<arc_length>() = dp;

		vd.add();
		vd.getLastPos()[0] = RightWall.get(0);
		vd.getLastPos()[1] = RightWall.get(1);
		vd.getLastPos()[2] = RightWall.get(2);

		vd.template getLastProp<type>() = BOUNDARY;
		vd.template getLastProp<pressure>() = 0.0;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<drho>() = 0.0;

		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;

		vd.template getLastProp<force_transport>()[0] = 0.0;
		vd.template getLastProp<force_transport>()[1] = 0.0;
		vd.template getLastProp<force_transport>()[2] = 0.0;

		vd.template getLastProp<v_transport>()[0] = 0.0;
		vd.template getLastProp<v_transport>()[1] = 0.0;
		vd.template getLastProp<v_transport>()[2] = 0.0;

		vd.template getLastProp<normal_vector>()[0] = 0.0;
		vd.template getLastProp<normal_vector>()[1] = 0.0;
		vd.template getLastProp<normal_vector>()[2] = 0.0;

		vd.template getLastProp<curvature_boundary>() = 0.0;
		vd.template getLastProp<arc_length>() = dp;
	}
}

particles CreateParticleGeometry(Vcluster<> &v_cl)
{

	// Physical size of the fluid domain, it goes from (0,0,0) to (length[0],length[1],length[2])
	// First particle will always be placed at (dp/2,dp/2,dp/2) and the last particle will be placed at (length[0]-dp/2,length[1]-dp/2,length[2]-dp/2)
	double length[3];

	// Size of the virtual cartesian grid that defines where to place the particles
	size_t sz[3];

	// In the case of the new bc we need particles at the wall, for this we need sz_aux
	// We want to put one virtual grid point between each pair of the old ones,
	// so that the new spacing is dp/2, and we can put a fluid particle exactly at the wall
	size_t sz_aux[3];
	size_t sz_aux_sq[3];

	// Boundary conditions
	size_t bc[3];

	// Number of boundary particles in each direction
	size_t Nboundary[3];

	// Number of fluid particles in each direction
	size_t Nfluid[3];

	double CylinderRadius;

	if (SCENARIO == POISEUILLE)
	{
		Nfluid[0] = 20;
		Nfluid[1] = 40;
		Nfluid[2] = 1;
		dp = 1.0 / Nfluid[1];
		H = Hconst * dp;
		r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
		Kcubic = (dimensions == 3) ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;
		Kquintic = (dimensions == 3) ? 1.0 / 120.0 / M_PI / H / H / H : 7.0 / 478.0 / M_PI / H / H;

		rho_zero = 1.0;
		nu = 0.01;
		eta = nu * rho_zero;

		MassFluid = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);
		MassBound = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);

		gravity_vector = {0.1, 0.0, 0.0};
		gravity = sqrt(norm2(gravity_vector));
		vw_top = {0.0, 0.0, 0.0};
		vw_bottom = {0.0, 0.0, 0.0};

		umax = gravity_vector.get(0) * 1.0 * 1.0 / (8.0 * nu);
		cbar = coeff_sound * umax;
		B = gamma_ * cbar * cbar / rho_zero;
		Pbackground = Bfactor * B;
		Re = umax * 1.0 / nu;

		bc[0] = PERIODIC;
		bc[1] = NON_PERIODIC;
		bc[2] = NON_PERIODIC;
		size_t Nbound_x = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = 0;
		Nboundary[1] = Nbound_x;
		Nboundary[2] = 0;
	}
	else if (SCENARIO == COUETTE)
	{
		Nfluid[0] = 20;
		Nfluid[1] = 40;
		Nfluid[2] = 1;
		dp = 1.0 / Nfluid[1];
		H = Hconst * dp;
		r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
		Kcubic = (dimensions == 3) ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;
		Kquintic = (dimensions == 3) ? 1.0 / 120.0 / M_PI / H / H / H : 7.0 / 478.0 / M_PI / H / H;

		rho_zero = 1.0;
		nu = 0.01;
		eta = nu * rho_zero;

		MassFluid = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);
		MassBound = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);

		gravity_vector = {0.0, 0.0, 0.0};
		gravity = sqrt(norm2(gravity_vector));
		vw_top = {0.0, 1.25, 0.0};
		vw_bottom = {0.0, 0.0, 0.0};

		umax = vw_top.get(1);
		cbar = coeff_sound * umax;
		B = gamma_ * cbar * cbar / rho_zero;
		Pbackground = Bfactor * B;
		Re = umax * 1.0 / nu;

		bc[0] = PERIODIC;
		bc[1] = NON_PERIODIC;
		bc[2] = NON_PERIODIC;
		size_t Nbound_x = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = 0;
		Nboundary[1] = Nbound_x;
		Nboundary[2] = 0;
	}
	else if (SCENARIO == HYDROSTATIC)
	{
		Nfluid[0] = 40;
		Nfluid[1] = 40;
		Nfluid[2] = 1;
		dp = 1.0 / Nfluid[0];
		H = Hconst * dp;
		r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
		Kcubic = (dimensions == 3) ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;
		Kquintic = (dimensions == 3) ? 1.0 / 120.0 / M_PI / H / H / H : 7.0 / 478.0 / M_PI / H / H;

		rho_zero = 1.0;
		nu = 0.1;
		eta = nu * rho_zero;

		MassFluid = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);
		MassBound = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);

		gravity_vector = {0.0, 0.1, 0.0};
		gravity = sqrt(norm2(gravity_vector));
		vw_top = {0.0, 0.0, 0.0};
		vw_bottom = {0.0, 0.0, 0.0};

		umax = 1;
		cbar = coeff_sound * umax;
		B = gamma_ * cbar * cbar / rho_zero;
		Pbackground = Bfactor * B;
		Re = umax * 1.0 / nu;

		bc[0] = NON_PERIODIC;
		bc[1] = NON_PERIODIC;
		bc[2] = NON_PERIODIC;
		size_t Nbound_x = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = Nbound_x;
		Nboundary[1] = Nbound_x;
		Nboundary[2] = 0;
	}
	else if (SCENARIO == CYLINDER_ARRAY)
	{

		// TESTED PARAMETERS
		Nfluid[0] = 60;
		Nfluid[1] = 40;
		Nfluid[2] = 1;
		dp = 1.0 / Nfluid[1];
		H = Hconst * dp;
		r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
		Kcubic = (dimensions == 3) ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;
		Kquintic = (dimensions == 3) ? 1.0 / 120.0 / M_PI / H / H / H : 7.0 / 478.0 / M_PI / H / H;
		CylinderRadius = (Nfluid[1] / 4.0) * H;

		rho_zero = 1.0;
		nu = 0.1;
		eta = nu * rho_zero;

		MassFluid = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);
		MassBound = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);

		gravity_vector = {0.1, 0.0, 0.0};

		gravity = sqrt(norm2(gravity_vector));
		vw_top = {0.0, 0.0, 0.0};
		vw_bottom = {0.0, 0.0, 0.0};

		umax = gravity_vector.get(0) * 1.0 * 1.0 / (8.0 * nu);

		cbar = coeff_sound * umax;
		B = gamma_ * cbar * cbar / rho_zero;
		Pbackground = Bfactor * B;

		Re = umax * 0.02 / nu;

		bc[0] = PERIODIC;
		bc[1] = NON_PERIODIC;
		bc[2] = NON_PERIODIC;
		size_t Nbound_x = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = 0;
		Nboundary[1] = Nbound_x;
		Nboundary[2] = 0;

		// PAPER PARAMETERS

		// Nfluid[0] = 72;
		// Nfluid[1] = 48;
		// Nfluid[2] = 1;
		// dp = 4.0 * 0.02 / Nfluid[1];
		// H = Hconst * dp;
		// r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
		// Kcubic = (dimensions == 3) ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;
		// Kquintic = (dimensions == 3) ? 1.0 / 120.0 / M_PI / H / H / H : 7.0 / 478.0 / M_PI / H / H;
		// CylinderRadius = (Nfluid[1] / 4.0) * dp;

		// rho_zero = 1000.0;
		// eta = 0.1;
		// nu = eta / rho_zero;

		// MassFluid = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);
		// MassBound = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);

		// gravity_vector = {2.5 * 1e-4, 0.0, 0.0};
		// gravity = getVectorNorm(gravity_vector);
		// vw_top = {0.0, 0.0, 0.0};
		// vw_bottom = {0.0, 0.0, 0.0};

		// // umax = gravity_vector.get(1) * 1.0 * 1.0 / (8.0 * nu);
		// cbar = 0.1 * std::sqrt(gravity * CylinderRadius);
		// B = gamma_ * cbar * cbar / rho_zero;
		// Pbackground = Bfactor * B;

		// Re = 2.4 * 1e-4;

		// bc[0] = NON_PERIODIC;
		// bc[1] = PERIODIC;
		// bc[2] = NON_PERIODIC;
		// size_t Nbound_x = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		// Nboundary[0] = 0;
		// Nboundary[1] = Nbound_x;
		// Nboundary[2] = 0;
	}
	else if (SCENARIO == CYLINDER_LATTICE)
	{

		// TESTED PARAMETERS
		Nfluid[0] = 60;
		Nfluid[1] = 40;
		Nfluid[2] = 1;
		dp = 1.0 / Nfluid[1];
		H = Hconst * dp;
		r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
		Kcubic = (dimensions == 3) ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;
		Kquintic = (dimensions == 3) ? 1.0 / 120.0 / M_PI / H / H / H : 7.0 / 478.0 / M_PI / H / H;
		CylinderRadius = (Nfluid[1] / 4.0) * H;

		rho_zero = 1.0;
		nu = 0.01;
		eta = nu * rho_zero;

		MassFluid = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);
		MassBound = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);

		gravity_vector = {0.1, 0.0, 0.0};

		gravity = sqrt(norm2(gravity_vector));
		vw_top = {0.0, 0.0, 0.0};
		vw_bottom = {0.0, 0.0, 0.0};

		umax = gravity_vector.get(0) * 1.0 * 1.0 / (8.0 * nu);
		cbar = coeff_sound * umax;
		B = gamma_ * cbar * cbar / rho_zero;
		Pbackground = Bfactor * B;

		Re = umax * 0.02 / nu;

		bc[0] = PERIODIC;
		bc[1] = PERIODIC;
		bc[2] = NON_PERIODIC;
		size_t Nbound_x = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = 0;
		Nboundary[1] = Nbound_x;
		Nboundary[2] = 0;
	}
	else if (SCENARIO == CYLINDER_ARRAY_LONG)
	{

		// TESTED PARAMETERS
		Nfluid[0] = 120;
		Nfluid[1] = 50;
		Nfluid[2] = 1;
		dp = 1.0 / Nfluid[1];
		H = Hconst * dp;
		r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
		Kcubic = (dimensions == 3) ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;
		Kquintic = (dimensions == 3) ? 1.0 / 120.0 / M_PI / H / H / H : 7.0 / 478.0 / M_PI / H / H;
		CylinderRadius = (Nfluid[1] / 4.0) * H;

		rho_zero = 1.0;
		nu = 0.01;
		eta = nu * rho_zero;

		MassFluid = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);
		MassBound = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);

		gravity_vector = {0.1, 0.0, 0.0};

		gravity = sqrt(norm2(gravity_vector));
		vw_top = {0.0, 0.0, 0.0};
		vw_bottom = {0.0, 0.0, 0.0};

		umax = gravity_vector.get(0) * 1.0 * 1.0 / (8.0 * nu);
		cbar = coeff_sound * umax;
		B = gamma_ * cbar * cbar / rho_zero;
		Pbackground = Bfactor * B;

		Re = umax * 0.02 / nu;

		bc[0] = PERIODIC;
		bc[1] = NON_PERIODIC;
		bc[2] = NON_PERIODIC;
		size_t Nbound_x = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = 0;
		Nboundary[1] = Nbound_x;
		Nboundary[2] = 0;
	}
	else if (SCENARIO == SQUARE)
	{
		Nfluid[0] = 60;
		Nfluid[1] = 40;
		Nfluid[2] = 1;
		dp = 1.0 / Nfluid[1];
		H = Hconst * dp;
		r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
		Kcubic = (dimensions == 3) ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;
		Kquintic = (dimensions == 3) ? 1.0 / 120.0 / M_PI / H / H / H : 7.0 / 478.0 / M_PI / H / H;

		rho_zero = 1.0;
		nu = 0.01;
		eta = nu * rho_zero;

		MassFluid = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);
		MassBound = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);

		gravity_vector = {0.1, 0.0, 0.0};
		gravity = sqrt(norm2(gravity_vector));
		vw_top = {0.0, 0.0, 0.0};
		vw_bottom = {0.0, 0.0, 0.0};

		umax = gravity_vector.get(0) * 1.0 * 1.0 / (8.0 * nu);
		cbar = coeff_sound * umax;
		B = gamma_ * cbar * cbar / rho_zero;
		Pbackground = Bfactor * B;
		Re = umax * 1.0 / nu;

		bc[0] = PERIODIC;
		bc[1] = NON_PERIODIC;
		bc[2] = NON_PERIODIC;
		size_t Nbound_x = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		Nboundary[0] = 0;
		Nboundary[1] = Nbound_x;
		Nboundary[2] = 0;
	}

	// Now define the boxes and the grid of the domain

	// We define the boxes in terms of offstes with respect to the fluid box that goes from 0 to length
	// offset to add to the domain box to create the correct particle positions
	double offset_domain[3] = {0.0, 0.0, 0.0};
	// offset to add to the recipient box, if non periodic it has to go from - Nboundary*dp  to length + Nboundary*dp
	double offset_recipient[3] = {0.0, 0.0, 0.0};
	// In case of periodic boundary conditions, we need to add an asymetric offset to the right, top or front of the domain
	double offset_periodic_fluid[3] = {0.0, 0.0, 0.0};
	// For the new BC we also need an offset
	double offset_periodic_recipient[3] = {0.0, 0.0, 0.0};

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

	for (int dim = 0; dim < 3; dim++)
	{

		if (bc[dim] == NON_PERIODIC) // non periodic, fluid covered by boundary
		{
			length[dim] = dp * Nfluid[dim];
			sz[dim] = Nfluid[dim] + 2 * (Nboundary[dim] + 1);
			offset_domain[dim] = (0.5 + Nboundary[dim]) * dp;

			if (Nboundary[dim] != 0)
				sz_aux[dim] = 2 * Nfluid[dim] - 1 + 2 * (2 * Nboundary[dim] + 1 + 1);
			else // for a direction with no boundary particles we dont need to add anything
				sz_aux[dim] = sz[dim];

			sz_aux_sq[dim] = sz_aux[dim];

			if (BC_TYPE == NEW_NO_SLIP) // Nboundary should only be 0 or 1 if we are using the new bc
				offset_recipient[dim] = 0.25 * Nboundary[dim] * dp;
			else if (BC_TYPE == NO_SLIP)
				offset_recipient[dim] = Nboundary[dim] * dp;
		}
		else // periodic, open ended
		{
			Nfluid[dim] -= 1;
			length[dim] = dp * Nfluid[dim];

			sz[dim] = Nfluid[dim] + 2;
			sz_aux[dim] = sz[dim];
			sz_aux_sq[dim] = 2 * Nfluid[dim] - 1 + 2 * (1 + 1);

			offset_domain[dim] = 0.5 * dp;
			offset_periodic_fluid[dim] = 0.75 * dp;
			offset_periodic_recipient[dim] = 0.85 * dp;
		}
	}

	// Define the boxes
	Box<3, double> domain({-offset_domain[0],
						   -offset_domain[1],
						   -offset_domain[2]},
						  {length[0] + offset_domain[0],
						   length[1] + offset_domain[1],
						   length[2] + offset_domain[2]});

	Box<3, double> fluid_box({0.0,
							  0.0,
							  0.0},
							 {length[0] + offset_periodic_fluid[0],
							  length[1] + offset_periodic_fluid[1],
							  length[2] + offset_periodic_fluid[2]});

	Box<3, double> recipient({-offset_recipient[0],
							  -offset_recipient[1],
							  -offset_recipient[2]},
							 {length[0] + offset_recipient[0] + offset_periodic_recipient[0],
							  length[1] + offset_recipient[1] + offset_periodic_recipient[1],
							  length[2] + offset_recipient[2] + offset_periodic_recipient[2]});

	// Will only be used in the new bc
	Box<3, double> recipient_hole({offset_recipient[0],
								   offset_recipient[1],
								   offset_recipient[2]},
								  {length[0] - offset_recipient[0] + offset_periodic_fluid[0],
								   length[1] - offset_recipient[1] + offset_periodic_fluid[1],
								   length[2] - offset_recipient[2] + offset_periodic_fluid[2]});

	// extended boundary around the domain, and the processor domain
	Ghost<3, double> g(r_threshold + H);

	// create particle object
	particles vd(0, domain, bc, g, DEC_GRAN(512));

	// correct the number of particles in case of periodicity, we substracted 1 before to accomodate the periodic boundary
	for (int dim = 0; dim < 3; dim++)
	{
		if (bc[dim] == PERIODIC)
		{
			Nfluid[dim] += 1;
			length[dim] += dp;
		}
	}

	// Write constants on file
	std::string customString = "4proc";
	WriteParameters(Nfluid, length, customString);

	// Set cylindrical object parameters
	Point<3, double> CylinderCentre = {length[0] / 2.0, length[1] / 2.0, length[2] / 2.0};

	// Set square obstacle parameters

	int integerSideLength = 9;
	Point<3, double> SquareCentre = CylinderCentre;
	Box<3, double> SquareOutside{{SquareCentre.get(0) - (((double)integerSideLength - 1.0) * dp) / 2.0,
								  SquareCentre.get(1) - (((double)integerSideLength - 1.0) * dp) / 2.0,
								  SquareCentre.get(2)},
								 {SquareCentre.get(0) + (((double)integerSideLength - 1.0) * dp) / 2.0,
								  SquareCentre.get(1) + (((double)integerSideLength - 1.0) * dp) / 2.0,
								  SquareCentre.get(2)}};

	// Box<3, double> SquareInside{{SquareCentre.get(0) - Ln2 * dp,
	// 							 SquareCentre.get(1) - Ln2 * dp,
	// 							 SquareCentre.get(2)},
	// 							{SquareCentre.get(0) + Ln2 * dp,
	// 							 SquareCentre.get(1) + Ln2 * dp,
	// 							 SquareCentre.get(2)}};

	if (SCENARIO == CYLINDER_ARRAY_LONG)
		CylinderCentre = {length[0] / 2.0, length[1] / 2.0, length[2] / 6.0};

	if (BC_TYPE == NEW_NO_SLIP && (SCENARIO == CYLINDER_ARRAY || SCENARIO == CYLINDER_LATTICE || SCENARIO == CYLINDER_ARRAY_LONG))
	{
		if (v_cl.getProcessUnitID() == 0) // Add the cylinder as boundary particles only on processor 0
		{
			AddCylinderNewBC(vd, CylinderCentre, CylinderRadius);
		}
	}

	if (BC_TYPE == NEW_NO_SLIP && (SCENARIO == SQUARE))
	{
		if (v_cl.getProcessUnitID() == 0) // Add the square as boundary particles only on processor 0
		{
			AddSquareNewBC(vd, SquareCentre, integerSideLength);
		}
	}

	Sphere<3, double> Cylinder(CylinderCentre, CylinderRadius);
	Sphere<3, double> Cylinder_aux(CylinderCentre, CylinderRadius + 0.5 * dp);
	Sphere<3, double> *Cylinder_ptr = nullptr;

	if (BC_TYPE == NEW_NO_SLIP)
		Cylinder_ptr = &Cylinder_aux;
	else if (BC_TYPE == NO_SLIP)
		Cylinder_ptr = &Cylinder;

	const double Lx = length[0] + 0.5 * dp; // channel height
	const double Ly = length[1] + 0.5 * dp; // channel width
	const double Lz = length[2] + 0.5 * dp; // channel length ( +0.5dp due to periodicity)

	// return an iterator to the fluid particles to add to vd
	auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

	// for each particle inside the fluid box ...
	while (fluid_it.isNext())
	{

		Point<3, double> iterator_posion = fluid_it.get();

		if (SCENARIO == CYLINDER_ARRAY || SCENARIO == CYLINDER_LATTICE || SCENARIO == CYLINDER_ARRAY_LONG)
		{
			if ((*Cylinder_ptr).isInside(iterator_posion)) // if inside the cylinder region
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
			else // if not inside just add fluid particles as usual
			{
				// ... add a particle ...
				vd.add();
				vd.template getLastProp<type>() = FLUID;
			}
		}
		else if (SCENARIO == SQUARE)
		{
			if (SquareOutside.isInside(iterator_posion)) // if inside the cylinder region
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
			else // if not inside just add fluid particles as usual
			{
				// ... add a particle ...
				vd.add();
				vd.template getLastProp<type>() = FLUID;
			}
		}
		else // if no cylinder at all just add fluid particles
		{
			// ... add a particle ...
			vd.add();
			vd.template getLastProp<type>() = FLUID;
		}

		// Set position
		// If periodic bc, we substract a small value to the last particle so that its not exactly at the boundary
		// otherwise we can get holes in the fluid region
		if (bc[0] == PERIODIC && iterator_posion.get(0) == Lx)
			vd.getLastPos()[0] = iterator_posion.get(0) - 1e-10 * dp;
		else
			vd.getLastPos()[0] = iterator_posion.get(0);

		if (bc[1] == PERIODIC && iterator_posion.get(1) == Ly)
			vd.getLastPos()[1] = iterator_posion.get(1) - 1e-10 * dp;
		else
			vd.getLastPos()[1] = iterator_posion.get(1);

		if (bc[2] == PERIODIC && iterator_posion.get(2) == Lz)
			vd.getLastPos()[2] = iterator_posion.get(2) - 1e-10 * dp;
		else
			vd.getLastPos()[2] = iterator_posion.get(2);

		// Set properties
		vd.template getLastProp<pressure>() = 0.0;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<drho>() = 0.0;

		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;

		vd.template getLastProp<force_transport>()[0] = 0.0;
		vd.template getLastProp<force_transport>()[1] = 0.0;
		vd.template getLastProp<force_transport>()[2] = 0.0;

		vd.template getLastProp<v_transport>()[0] = 0.0;
		vd.template getLastProp<v_transport>()[1] = 0.0;
		vd.template getLastProp<v_transport>()[2] = 0.0;

		// Not needed but filled anyway
		vd.template getLastProp<normal_vector>()[0] = 0.0;
		vd.template getLastProp<normal_vector>()[1] = 0.0;
		vd.template getLastProp<normal_vector>()[2] = 0.0;

		vd.template getLastProp<curvature_boundary>() = 0.0;
		vd.template getLastProp<arc_length>() = 0.0;

		// next fluid particle
		++fluid_it;
	}

	// Now place solid particles
	openfpm::vector<Box<3, double>> holes;
	// openfpm::vector<Box<3, double>> holes_square;
	// holes_square.add(SquareInside);

	if (BC_TYPE == NEW_NO_SLIP)
	{
		holes.add(recipient_hole);
		sz[0] = sz_aux[0];
		sz[1] = sz_aux[1];
		sz[2] = sz_aux[2];
	}
	else if (BC_TYPE == NO_SLIP)
		holes.add(fluid_box);

	Box<3, double> hole_get = holes.get(0);
	PointIteratorSkin<3, double, CartDecomposition<3, double>> bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);
	// PointIteratorSkin<3, double, CartDecomposition<3, double>> square_bound_box = DrawParticles::DrawSkin(vd, sz_aux_sq, domain, holes_square, SquareOutside);

	// if (SCENARIO == SQUARE && BC_TYPE == NEW_NO_SLIP)
	// {
	// 	while (square_bound_box.isNext())
	// 	{
	// 		Point<3, double> position = square_bound_box.get();

	// 		// if (BC_TYPE == NEW_NO_SLIP)
	// 		// {
	// 		// 	// Check if x and z coordinates are multiples of dp, keep multiples, discard the rest
	// 		// 	double remx = fmod(position.get(0), dp);
	// 		// 	double remz = fmod(position.get(1), dp);
	// 		// 	double tol = 0.5 * dp * 10e-2;

	// 		// 	if (remx > tol && remx < dp - tol)
	// 		// 	{
	// 		// 		++square_bound_box;
	// 		// 		continue;
	// 		// 	}
	// 		// 	if (remz > tol && remz < dp - tol)
	// 		// 	{
	// 		// 		++square_bound_box;
	// 		// 		continue;
	// 		// 	}
	// 		// }
	// 		vd.add();

	// 		vd.getLastPos()[0] = square_bound_box.get().get(0);
	// 		vd.getLastPos()[1] = square_bound_box.get().get(1);
	// 		vd.getLastPos()[2] = square_bound_box.get().get(2);

	// 		vd.template getLastProp<type>() = BOUNDARY;
	// 		vd.template getLastProp<rho>() = rho_zero;
	// 		vd.template getLastProp<pressure>() = 0.0;
	// 		vd.template getLastProp<drho>() = 0.0;

	// 		vd.template getLastProp<force>()[0] = 0.0;
	// 		vd.template getLastProp<force>()[1] = 0.0;
	// 		vd.template getLastProp<force>()[2] = 0.0;

	// 		vd.template getLastProp<force_transport>()[0] = 0.0;
	// 		vd.template getLastProp<force_transport>()[1] = 0.0;
	// 		vd.template getLastProp<force_transport>()[2] = 0.0;

	// 		vd.template getLastProp<v_transport>()[0] = 0.0;
	// 		vd.template getLastProp<v_transport>()[1] = 0.0;
	// 		vd.template getLastProp<v_transport>()[2] = 0.0;

	// 		vd.template getLastProp<normal_vector>()[0] = 0.0;
	// 		vd.template getLastProp<normal_vector>()[1] = 0.0;
	// 		vd.template getLastProp<normal_vector>()[2] = 0.0;

	// 		vd.template getLastProp<velocity>()[0] = vw_bottom.get(0);
	// 		vd.template getLastProp<velocity>()[1] = vw_bottom.get(1);
	// 		vd.template getLastProp<velocity>()[2] = vw_bottom.get(2);

	// 		vd.template getLastProp<velocity>()[0] = vw_top.get(0);
	// 		vd.template getLastProp<velocity>()[1] = vw_top.get(1);
	// 		vd.template getLastProp<velocity>()[2] = vw_top.get(2);

	// 		vd.template getLastProp<curvature_boundary>() = 0.0;

	// 		vd.template getLastProp<arc_length>() = dp;

	// 		++square_bound_box;
	// 	}
	// }
	if (bc[0] != PERIODIC || bc[1] != PERIODIC) // no walls in all periodic scenario
	{
		while (bound_box.isNext())
		{
			Point<3, double> position = bound_box.get();

			// periodic bc, with no boundary particles in y direction has a bug, it puts 3 extra particles outside in the y direction
			// When running on multiple cores, with this we check if particle is outside the recipient box
			// Another bug places boundary particles in the correct plane, but inside the fluid box;
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
			// if (bc[0] == PERIODIC && position.get(0) > dp / 2.0 && position.get(0) < length[0] - dp / 2.0)
			// {
			// 	++bound_box;
			// 	continue;
			// }

			if (BC_TYPE == NEW_NO_SLIP && (bc[0] == NON_PERIODIC && bc[1] == NON_PERIODIC && bc[2] == NON_PERIODIC))
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

			vd.getLastPos()[0] = bound_box.get().get(0);
			vd.getLastPos()[1] = bound_box.get().get(1);
			vd.getLastPos()[2] = bound_box.get().get(2);

			vd.template getLastProp<type>() = BOUNDARY;
			vd.template getLastProp<rho>() = rho_zero;
			vd.template getLastProp<pressure>() = 0.0;
			vd.template getLastProp<drho>() = 0.0;

			vd.template getLastProp<force>()[0] = 0.0;
			vd.template getLastProp<force>()[1] = 0.0;
			vd.template getLastProp<force>()[2] = 0.0;

			vd.template getLastProp<force_transport>()[0] = 0.0;
			vd.template getLastProp<force_transport>()[1] = 0.0;
			vd.template getLastProp<force_transport>()[2] = 0.0;

			vd.template getLastProp<v_transport>()[0] = 0.0;
			vd.template getLastProp<v_transport>()[1] = 0.0;
			vd.template getLastProp<v_transport>()[2] = 0.0;

			vd.template getLastProp<normal_vector>()[0] = 0.0;
			vd.template getLastProp<normal_vector>()[1] = 0.0;
			vd.template getLastProp<normal_vector>()[2] = 0.0;

			if (position.get(0) < dp / 4.0) // bottom wall
			{
				vd.template getLastProp<velocity>()[0] = vw_bottom.get(0);
				vd.template getLastProp<velocity>()[1] = vw_bottom.get(1);
				vd.template getLastProp<velocity>()[2] = vw_bottom.get(2);

				// vd.template getLastProp<normal_vector>()[0] += 1.0;
				// vd.template getLastProp<normal_vector>()[1] += 0.0;
				// vd.template getLastProp<normal_vector>()[2] += 0.0;

				vd.template getLastProp<curvature_boundary>() = 0.0;

				vd.template getLastProp<arc_length>() = dp;
			}
			if (position.get(0) > length[0] - dp / 4.0) // top wall
			{
				vd.template getLastProp<velocity>()[0] = vw_top.get(0);
				vd.template getLastProp<velocity>()[1] = vw_top.get(1);
				vd.template getLastProp<velocity>()[2] = vw_top.get(2);

				// vd.template getLastProp<normal_vector>()[0] += -1.0;
				// vd.template getLastProp<normal_vector>()[1] += 0.0;
				// vd.template getLastProp<normal_vector>()[2] += 0.0;

				vd.template getLastProp<curvature_boundary>() = 0.0;

				vd.template getLastProp<arc_length>() = dp;
			}
			if (position.get(1) < dp / 4.0) // left wall
			{
				vd.template getLastProp<velocity>()[0] = 0.0;
				vd.template getLastProp<velocity>()[1] = 0.0;
				vd.template getLastProp<velocity>()[2] = 0.0;

				// vd.template getLastProp<normal_vector>()[0] += 0.0;
				// vd.template getLastProp<normal_vector>()[1] += 1.0;
				// vd.template getLastProp<normal_vector>()[2] += 0.0;

				vd.template getLastProp<curvature_boundary>() = 0.0;

				vd.template getLastProp<arc_length>() = dp;
			}
			if (position.get(1) > length[1] - dp / 4.0) // right wall
			{
				vd.template getLastProp<velocity>()[0] = 0.0;
				vd.template getLastProp<velocity>()[1] = 0.0;
				vd.template getLastProp<velocity>()[2] = 0.0;

				// vd.template getLastProp<normal_vector>()[0] += 0.0;
				// vd.template getLastProp<normal_vector>()[1] += -1.0;
				// vd.template getLastProp<normal_vector>()[2] += 0.0;

				vd.template getLastProp<curvature_boundary>() = 0.0;

				vd.template getLastProp<arc_length>() = dp;
			}

			// lower left corner
			if (position.get(0) < dp / 4.0 && position.get(1) < dp / 4.0)
			{

				vd.template getLastProp<curvature_boundary>() = 0.0 * dp;
			}
			// lower right corner
			if (position.get(0) > length[0] - dp / 4.0 && position.get(1) < dp / 4.0)
			{
				vd.template getLastProp<curvature_boundary>() = 0.0 * dp;
			}
			// upper left corner
			if (position.get(0) < dp / 4.0 && position.get(1) > length[1] - dp / 4.0)
			{
				vd.template getLastProp<curvature_boundary>() = 0.0 * dp;
			}
			// upper right corner
			if (position.get(0) > length[0] - dp / 4.0 && position.get(1) > length[1] - dp / 4.0)
			{
				vd.template getLastProp<curvature_boundary>() = 0.0 * dp;
			}

			// normalize for corner particles that have normal vector with norm 2
			// Point<3, double> normalvec = vd.template getLastProp<normal_vector>();
			// double norm = getVectorNorm(normalvec);
			// vd.template getLastProp<normal_vector>()[0] /= norm;
			// vd.template getLastProp<normal_vector>()[1] /= norm;
			// vd.template getLastProp<normal_vector>()[2] /= norm;

			++bound_box;
		}
	}
	return vd;
}

int main(int argc, char *argv[])
{

	// initialize the library
	openfpm_init(&argc, &argv);
	Vcluster<> &v_cl = create_vcluster();
	particles vd = CreateParticleGeometry(v_cl);
	vd.write("FIRST");

	vd.map();
	vd.write("SECOND");

	openfpm::vector<std::string> names({"type", "rho", "pressure", "drho", "force", "velocity", "force_transport", "v_transport", "normal", "curvature", "arc_length"});
	vd.setPropNames(names);
	// Now that we fill the vector with particles
	ModelCustom md;

	// vd.addComputationCosts(md);
	vd.getDecomposition().decompose();
	vd.map();

	vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();

	vd.write("THIRD");

	auto NN = vd.getCellList(r_threshold + H);

	// if (DENSITY_TYPE == DENSITY_SUMMATION)
	// {
	// 	// Adjust mass to have the correct density
	// 	fix_mass(vd, NN);
	// }

	// Evolve
	size_t write = 0;
	size_t it = 0;
	size_t it_reb = 0;
	double t = 0.0;

	vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();
	v_cl.barrier();
	if (BC_TYPE == NO_SLIP) // set up boundary particle velocity
	{
		calc_boundary(vd, NN);
	}
	else if (BC_TYPE == NEW_NO_SLIP) // Set up fluid vector and normal vector of the boundary particles
	{
		calcFluidVec(vd, NN);
		vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();
		calcNormalVec(vd, NN);
		vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();
		calcCurvature(vd, NN);
	}
	vd.write("CURVATURE");

	// Compute forces for the first iteration
	vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();
	calc_forces(vd, NN);
	vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();

	while (t <= t_end)
	{
		timer it_time;

		//// Do rebalancing every 200 timesteps
		it_reb++;
		if (it_reb == 1000000)
		{
			vd.map();

			it_reb = 0;
			ModelCustom md;
			vd.addComputationCosts(md);
			vd.getDecomposition().decompose();
			// just get decomposition
			vd.map();

			if (v_cl.getProcessUnitID() == 0)
				std::cout << "REBALANCED " << std::endl;
		}
		vd.map();

		// Calculate dt for time stepping
		double dt = calc_deltaT(vd);

		// Integrate one time step
		kick_drift_int(vd, NN, dt, v_cl);

		// increment time
		t += dt;

		if (write < t * write_const)
		{
			// vd.map();

			vd.deleteGhost();
			vd.write_frame(filename, write, WRITER);
			vd.ghost_get<type, rho, pressure, drho, force, velocity, force_transport, v_transport, normal_vector, curvature_boundary, arc_length>();

			write++;
			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "TIME: " << t << "  write " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "   " << cnt << std::endl;
			}
		}
		else
		{
			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "TIME: " << t << "  " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "   " << cnt << std::endl;
			}
		}
	}

	openfpm_finalize();
}
