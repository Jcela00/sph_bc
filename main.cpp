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

// TYPE OF KERNERL
#define CUBIC 0
#define QUINTIC 1

//////// CODE PARAMETERS /////////////////////////////////////////////////////////////
// Dimensionality of the problem changes the normalizaiton of the kernel
const int dimensions = 2;

// BC and viscosity we are using in the simulation
const int BC_TYPE = NEW_NO_SLIP;
const int SCENARIO = POISEUILLE;
const int KERNEL = QUINTIC;
const int WRITER = VTK_WRITER; // VTK_WRITER or CSV_WRITER

// Output file name, filled later
std::string filename;

////////////////////////////////////////////////////////////////////////////////////

//////// SPATIAL CONSTANTS /////////////////////////////////////////////////////////
// initial spacing between particles dp in the formulas
const double dp = 0.025;

// Factor relating H (smoothing length) and dp ( particle spacing)
const double Hconst = 1.0;

// H, smoothing length
const double H = Hconst * dp;

const double r_threshold = (KERNEL == CUBIC ? 2.0 * H : 3.0 * H);
////////////////////////////////////////////////////////////////////////////////////

//////// EQ OF STATE CONSTANTS //////////////////////////////////////////////////////
// Eq of state constant, filled later
double B = 0.0;

// Constant used for the sound speed, number of times the max velocity
const double coeff_sound = 10.0;

// Sound speed (calculated later)
double cbar = 0.0;

// gamma in the formulas
const double gamma_ = 1.0;
const double xi = 0.0;

////////////////////////////////////////////////////////////////////////////////////

//////// PHYSICAL CONSTANTS ////////////////////////////////////////////////////////
// Reference density
const double rho_zero = 1.0;

// Minimum, Maximum Rho allowed
const double RhoMin = 0.7 * rho_zero;
const double RhoMax = 1.3 * rho_zero;

// Mass of the fluid and boundary particles M=rho*V, V=dp^3 or V = dp^2
const double MassOld = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);
double MassFluid = MassOld;
double MassBound = MassOld;

// Gravity vector and magnitude
Point<3, double> gravity_vector;
double gravity;

// Wall Normal
Point<3, double> normal_bottom_wall = {1.0, 0.0, 0.0};

// Wall velocity
Point<3, double> vw_top;
Point<3, double> vw_bottom;

// Tensile correction constant, filled later
double W_dap = 0.0;

// Viscosity we want to use, independent of the type of viscosity term
const double nu = 0.1;

// Viscosity coefficient of the artificial viscosity (alpha) filled later if using artificial viscosity
double visco = 500;

// Dynamic viscosity, used when using physical viscosity term
double eta = nu * rho_zero;

// Eta in the formulas
const double Eta2 = 0.01 * H * H;

const double initial_perturbation = 0.0;

double Pbackground = 10; // filled later
const double Bfactor = 1.0;

////////////////////////////////////////////////////////////////////////////////////

//////// TIME CONSTANTS //////////////////////////////////////////////////////////////
// End simulation time
const double t_end = 10;
// Constant used to define time integration
const double CFLnumber = 0.1;
// Minimum T
const double DtMin = 0.00001;
// Controls otput file frequency, low means less frequent
const int write_const = 10;
////////////////////////////////////////////////////////////////////////////////////

//////// ALIAS FOR THE PARTICLE PROPERTIES //////////////////////////////////////////
// Properties
// FLUID or BOUNDARY
const size_t type = 0;
// Density
const int rho = 1;
// Density at step n-1
const int rho_prev = 2;
// Pressure
const int Pressure = 3;
// Delta rho calculated in the force calculation
const int drho = 4;
// calculated force
const int force = 5;
// velocity
const int velocity = 6;
// velocity at previous step
const int velocity_prev = 7;
// Background pressure force
const int Background_P = 8;
//
const int v_transport = 9;

typedef vector_dist<3, double, aggregate<size_t, double, double, double, double, double[3], double[3], double[3], double[3], double[3]>> particles;
//                                       |     |        |        |        |        |         |            |		     |          |
//                                       |     |        |        |        |        |         |            |		     |			|
//                                     type  density   density  Pressure delta   force     velocity    velocity	  Pb force   v_transport
//                                                      at n-1           density                       at n - 1

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
inline void EqState(particles &vd)
{
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto a = it.get();

		double rho_a = vd.template getProp<rho>(a);
		double rho_frac = rho_a / rho_zero;

		vd.template getProp<Pressure>(a) = B * (std::pow(rho_frac, gamma_) - 1.0) + xi;

		++it;
	}
}

// Inverted equation of state, compute density given pressure, particle wise
inline double InvEqState(const double pressure)
{
	return rho_zero * std::pow(((pressure - xi) / B + 1.0), 1.0 / gamma_);
}

// Vector utilities
inline double dotProduct(const Point<3, double> &v, const Point<3, double> &w)
{
	return v.get(0) * w.get(0) + v.get(1) * w.get(1) + v.get(2) * w.get(2);
}
inline double norm_own(const Point<3, double> &v)
{
	return sqrt(v.get(0) * v.get(0) + v.get(1) * v.get(1) + v.get(2) * v.get(2));
}

// Kernel functions
// Normalization constant for the kernel
const double Kcubic = (dimensions == 3) ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;
const double Kquintic = (dimensions == 3) ? 1.0 / 120.0 / M_PI / H / H / H : 7.0 / 478.0 / M_PI / H / H;

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
	if (q >= 0.0 && q < 1.0)
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

	if (q >= 0.0 && q < 1.0)
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

// OLD TENSILE CONSTANTS 0.01:-0.2

inline double Tensile(const double r, const double rhoa, const double rhob, const double prs1, const double prs2)
{
	// Evaluate kernel
	const double wab = Wab(r);

	// Tensile correction.
	double fab = wab * W_dap;
	fab *= fab;
	fab *= fab;					 // fab=fab^4
	const double epsilon = 0.06; // (Monaghan 1999) for n=4 H = 1.5dp
	const double tensilp1 = (prs1 / (rhoa * rhoa)) * (prs1 > 0.0 ? epsilon : 0.00);
	const double tensilp2 = (prs2 / (rhob * rhob)) * (prs2 > 0.0 ? epsilon : 0.00);

	return (fab * (tensilp1 + tensilp2));
}

inline double Pi_artificial(const Point<3, double> &dr, const double rr2, const Point<3, double> &dv, const double rhoa, const double rhob, const double massb)
{
	const double dot = dotProduct(dr, dv);
	const double dot_rr2 = dot / (rr2 + Eta2);

	const float amubar = H * dot_rr2;
	const float robar = (rhoa + rhob) * 0.5;
	const float pi_visc = (-visco * cbar * amubar / robar);

	return pi_visc;
}

inline Point<3, double> Pi_physical(const Point<3, double> &dr, const double &r, const Point<3, double> &dv, const Point<3, double> &dW)
{
	return eta * (dv * dotProduct(dr, dW)) / r / r / 1.44;
}

inline double PressureForce(const double &rhoa, const double &rhob, const double &prsa, const double &prsb)
{
	return -1.0 * (rhob * prsa + rhoa * prsb) / (rhoa + rhob) / 1.44;
}

std::array<Point<3, double>, 3> GetDummyPositions(const Point<3, double> &r, const Point<3, double> &normal)
{
	Point<3, double> offset_1 = -0.5 * normal * dp; // offset from wall to first particle
	Point<3, double> offset_2 = -1.0 * normal * dp; // offset from first to second particle, and from second to third
	Point<3, double> r1 = r + offset_1;				// first dummy particle
	Point<3, double> r2 = r1 + offset_2;			// second dummy particle
	Point<3, double> r3 = r2 + offset_2;			// third dummy particle

	std::array<Point<3, double>, 3> r_dummy = {r1, r2, r3};

	return r_dummy;
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

				// if (a == b) skip this particle
				if (a.getKey() == b)
				{
					++Np;
					continue;
				};

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
							bool is_bottom_wall = (xb.get(0) < dp ? true : false);
							Point<3, double> normal = (is_bottom_wall ? normal_bottom_wall : -1.0 * normal_bottom_wall);

							// Apply offsets to dr to get 3 vectrors pointing to dummy particles
							std::array<Point<3, double>, 3> R_dummy = GetDummyPositions(-1.0 * dr, normal);

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
	// This function fills the value of velocity_prev for the boundary particles, with the velocity for no slip BC

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

				// if (a == b) skip this particle
				if (a.getKey() == b)
				{
					++Np;
					continue;
				};

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

						rho_sum += w;
					}
					else
					{
						if (BC_TYPE == NO_SLIP)
						{
							// calculate distance
							const double r = sqrt(r2);

							// evaluate kernel
							const double w = Wab(r);

							rho_sum += w;
						}
						else if (BC_TYPE == NEW_NO_SLIP) // need to evaluate kernel at dummy particles
						{
							// xw.get(0) this is the x coordinate of the wall
							bool is_bottom_wall = (xb.get(0) < dp ? true : false);
							Point<3, double> normal = (is_bottom_wall ? normal_bottom_wall : -1.0 * normal_bottom_wall);

							// Apply offsets to dr to get 3 vectrors pointing to dummy particles
							std::array<Point<3, double>, 3> R_dummy = GetDummyPositions(-1.0 * dr, normal);

							const Point<3, double> r1 = R_dummy[0];
							const Point<3, double> r2 = R_dummy[1];
							const Point<3, double> r3 = R_dummy[2];

							const double W1 = Wab(sqrt(norm2(r1)));
							const double W2 = Wab(sqrt(norm2(r2)));
							const double W3 = Wab(sqrt(norm2(r3)));

							rho_sum += (W1 + W2 + W3);
						}
					}
				}

				++Np;
			}
			if (rho_sum != 0.0)
			{
				vd.template getProp<rho>(a) = rho_sum * massa;
				// std::cout << "sum w: " << rho_sum << "  rho_sum*dp^2: " << rho_sum * H * H << " density: " << vd.template getProp<rho>(a) << std::endl;
				// std::cout << "Density accounted for " << countpart << " particles" << std::endl;
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
	// This function fills the value of velocity_prev for the boundary particles, with the velocity for no slip BC

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
				double Pf = vd.getProp<Pressure>(f);

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
				vd.template getProp<velocity_prev>(b)[0] = 2.0 * vd.template getProp<velocity>(b)[0] - sum_vW.get(0) / sum_W;
				vd.template getProp<velocity_prev>(b)[1] = 2.0 * vd.template getProp<velocity>(b)[1] - sum_vW.get(1) / sum_W;
				vd.template getProp<velocity_prev>(b)[2] = 2.0 * vd.template getProp<velocity>(b)[2] - sum_vW.get(2) / sum_W;
				// Set the pressure of the boundary particle b
				vd.template getProp<Pressure>(b) = sum_pW / sum_W;
				// Compute density from inverted Eq of state
				vd.template getProp<rho>(b) = InvEqState(vd.template getProp<Pressure>(b));
			}
			else
			{
				vd.template getProp<velocity_prev>(b)[0] = 2.0 * vd.template getProp<velocity>(b)[0];
				vd.template getProp<velocity_prev>(b)[1] = 2.0 * vd.template getProp<velocity>(b)[1];
				vd.template getProp<velocity_prev>(b)[2] = 2.0 * vd.template getProp<velocity>(b)[2];
				vd.template getProp<Pressure>(b) = 0.0;
				vd.template getProp<rho>(b) = rho_zero;
			}
		}

		++part;
	}
}

// function to know the type of variables currently at auto
template <class T>
std::string
type_name()
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

	const double massw = MassBound;
	const double rhow = vd.getProp<rho>(boundary_key);
	const double Pw = vd.getProp<Pressure>(boundary_key);
	const Point<3, double> vw = vd.getProp<velocity>(boundary_key);
	const Point<3, double> vtf = vd.getProp<v_transport>(fluid_key);
	const Point<3, double> vdiff_f = vtf - vf;

	// xw.get(0) this is the x coordinate of the wall
	bool is_bottom_wall = (xw.get(0) < dp ? true : false);
	Point<3, double> normal = (is_bottom_wall ? normal_bottom_wall : -1.0 * normal_bottom_wall);
	Point<3, double> tangential = {normal.get(2), 0.0, normal.get(0)};
	if (!is_bottom_wall)
	{
		tangential.get(2) = -tangential.get(2);
	}
	Point<3, double> r_fluid_to_wall = -1.0 * r_wall_to_fluid; // Points to wall from fluid
	// r_fluid_to_wall.get(1) = 0.0;

	// double vt = (vf.get(0) * tangential.get(0) + vf.get(1) * tangential.get(1) + vf.get(2) * tangential.get(2));
	// double vn = (vf.get(0) * normal.get(0) + vf.get(1) * normal.get(1) + vf.get(2) * normal.get(2));

	// Project va on tangential and normal directions
	double vt = dotProduct(vf, tangential);
	double vn = dotProduct(vf, normal);

	double lf = dotProduct(r_fluid_to_wall, normal); // vertical distance from fluid particle to wall
	lf = (lf < 0.0 ? -1.0 * lf : lf);				 // absolute value

	std::array<Point<3, double>, 3> r_dummy = GetDummyPositions(r_fluid_to_wall, normal);

	const Point<3, double> r1 = r_dummy[0]; // Vector from fluid particle to first dummy particle
	double r1_norm = sqrt(norm2(r1));		// Norm of r1
	double lw1 = 0.5 * dp;					// Distance from wall to first dummy particle

	const Point<3, double> r2 = r_dummy[1]; // Vector from fluid particle to second dummy particle
	double r2_norm = sqrt(norm2(r2));		// Norm of r2
	double lw2 = 1.5 * dp;					// Distance from wall to second dummy particle

	const Point<3, double> r3 = r_dummy[2]; // Vector from fluid particle to third dummy particle
	double r3_norm = sqrt(norm2(r3));		// Norm of r3
	double lw3 = 2.5 * dp;					// Distance from wall to third dummy particle

	Point<3, double> v1 = {0.0, 0.0, 0.0}; // Velocity of first dummy particle
	Point<3, double> v2 = {0.0, 0.0, 0.0}; // Velocity of second dummy particle
	Point<3, double> v3 = {0.0, 0.0, 0.0}; // Velocity of third dummy particle
	double p1, p2, p3;					   // Pressure of the dummy particles
	double rho1, rho2, rho3;			   // Density of the dummy particles

	std::array<double, 3> p_dummy = {0.0, 0.0, 0.0};
	std::array<double, 3> rho_dummy = {0.0, 0.0, 0.0};
	std::array<Point<3, double>, 3> DW_dummy = {Point<3, double>{0.0, 0.0, 0.0}, Point<3, double>{0.0, 0.0, 0.0}, Point<3, double>{0.0, 0.0, 0.0}};
	std::array<Point<3, double>, 3> v_dummy = {Point<3, double>{0.0, 0.0, 0.0}, Point<3, double>{0.0, 0.0, 0.0}, Point<3, double>{0.0, 0.0, 0.0}};

	double g_normal = dotProduct(gravity_vector, normal);
	int nonzeros = 0;
	if (r1_norm < r_threshold)
	{
		v1 = 2.0 * vw - vt * (lw1 / lf) * tangential - vn * (lw1 / lf) * normal;
		p1 = Pf + rhof * g_normal * (lf + lw1);
		rho1 = InvEqState(p1);

		r_dummy[0] = -1.0 * r_dummy[0]; // Force routines use the vector pointing from b to a
		Point<3, double> DW = DWab(r_dummy[0], r1_norm);

		p_dummy[0] = p1;
		rho_dummy[0] = rho1;
		v_dummy[0] = v1;
		DW_dummy[0] = DW;

		nonzeros++;
	}
	else if (r2_norm < r_threshold)
	{
		v2 = 2.0 * vw - vt * (lw2 / lf) * tangential - vn * (lw2 / lf) * normal;
		p2 = Pf + rhof * g_normal * (lf + lw2);
		rho2 = InvEqState(p2);

		r_dummy[1] = -1.0 * r_dummy[1];
		Point<3, double> DW = DWab(r_dummy[1], r2_norm);

		p_dummy[1] = p2;
		rho_dummy[1] = rho2;
		v_dummy[1] = v2;
		DW_dummy[1] = DW;

		nonzeros++;
	}
	else if (r3_norm < r_threshold)
	{
		v3 = 2.0 * vw - vt * (lw3 / lf) * tangential - vn * (lw3 / lf) * normal;
		p3 = Pf + rhof * g_normal * (lf + lw3);
		rho3 = InvEqState(p3);

		r_dummy[2] = -1.0 * r_dummy[2];
		Point<3, double> DW = DWab(r_dummy[2], r3_norm);

		p_dummy[2] = p3;
		rho_dummy[2] = rho3;
		v_dummy[2] = v3;
		DW_dummy[2] = DW;

		nonzeros++;
	}

	const double GradATerm = 0.5 * (rhof * dotProduct(vf, vdiff_f));

	for (size_t layer = 0; layer < nonzeros; layer++)
	{
		Point<3, double> v_rel = vf - v_dummy[layer];
		Point<3, double> DW = DW_dummy[layer];
		const double Va2 = (massf / rhof) * (massf / rhof);
		const double Vb2 = (MassBound / rho_dummy[layer]) * (MassBound / rho_dummy[layer]);

		Point<3, double> ViscosityTerm = Pi_physical(r_dummy[layer], sqrt(norm2(r_dummy[layer])), v_rel, DW_dummy[layer]);
		double PressureTerm = PressureForce(rhof, rho_dummy[layer], Pf, p_dummy[layer]);

		vd.getProp<force>(fluid_key)[0] += (Va2 + Vb2) * ((GradATerm + PressureTerm) * DW.get(0) + ViscosityTerm.get(0)) / massf;
		vd.getProp<force>(fluid_key)[1] += (Va2 + Vb2) * ((GradATerm + PressureTerm) * DW.get(1) + ViscosityTerm.get(1)) / massf;
		vd.getProp<force>(fluid_key)[2] += (Va2 + Vb2) * ((GradATerm + PressureTerm) * DW.get(2) + ViscosityTerm.get(2)) / massf;

		vd.getProp<Background_P>(fluid_key)[0] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(0) / massf;
		vd.getProp<Background_P>(fluid_key)[1] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(1) / massf;
		vd.getProp<Background_P>(fluid_key)[2] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(2) / massf;

		// MassBound * (vf.get(0) * (W_dummy[comp]).get(0) + vf.get(1) * (W_dummy[comp]).get(1) + vf.get(2) * (W_dummy[comp]).get(2));
		// vd.getProp<drho>(fluid_key) += MassBound * dotProduct(vf, DW_dummy[comp]);
		// vd.getProp<drho>(fluid_key) += MassBound * dotProduct(v_rel, DW_dummy[comp]);
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
	const double Pb = vd.getProp<Pressure>(b_key);
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

	vd.getProp<Background_P>(a_key)[0] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(0) / massa;
	vd.getProp<Background_P>(a_key)[1] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(1) / massa;
	vd.getProp<Background_P>(a_key)[2] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(2) / massa;

	// vd.getProp<drho>(a_key) += (massb) * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
	// vd.getProp<drho>(a_key) += massb * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
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
	const double Pb = vd.getProp<Pressure>(boundary_key);
	const Point<3, double> vb = vd.getProp<velocity>(boundary_key);
	const Point<3, double> vb_noslip = vd.getProp<velocity_prev>(boundary_key);

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

	vd.getProp<Background_P>(fluid_key)[0] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(0) / massf;
	vd.getProp<Background_P>(fluid_key)[1] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(1) / massf;
	vd.getProp<Background_P>(fluid_key)[2] += -1.0 * (Va2 + Vb2) * (Pbackground)*DW.get(2) / massf;

	// vd.getProp<drho>(a_key) += (massb) * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
	// vd.getProp<drho>(a_key) += massb * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
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
			double Pa = vd.getProp<Pressure>(a);

			// Get the Velocity of the particle a
			Point<3, double> va = vd.getProp<velocity>(a);

			// Reset the force counter (0 + gravity)
			vd.template getProp<force>(a)[0] = gravity_vector.get(0);
			vd.template getProp<force>(a)[1] = gravity_vector.get(1);
			vd.template getProp<force>(a)[2] = gravity_vector.get(2);

			vd.template getProp<Background_P>(a)[0] = 0.0;
			vd.template getProp<Background_P>(a)[1] = 0.0;
			vd.template getProp<Background_P>(a)[2] = 0.0;

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
	double dt_g = 0.25 * sqrt(H / gravity_vector.norm());
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

		vd.template getProp<v_transport>(a)[0] = vd.template getProp<velocity>(a)[0] + dt_2 * vd.template getProp<Background_P>(a)[0];
		vd.template getProp<v_transport>(a)[1] = vd.template getProp<velocity>(a)[1] + dt_2 * vd.template getProp<Background_P>(a)[1];
		vd.template getProp<v_transport>(a)[2] = vd.template getProp<velocity>(a)[2] + dt_2 * vd.template getProp<Background_P>(a)[2];

		vd.getPos(a)[0] += dt * vd.template getProp<v_transport>(a)[0];
		vd.getPos(a)[1] += dt * vd.template getProp<v_transport>(a)[1];
		vd.getPos(a)[2] += dt * vd.template getProp<v_transport>(a)[2];

		++part;
	}
	// v_cl.barrier();
	vd.map();
	// v_cl.barrier();
	vd.ghost_get<type, rho, Pressure, velocity, velocity_prev, v_transport>();
	// v_cl.barrier();

	calc_density(vd, NN);
	// Calculate pressure from the density
	EqState(vd);
	// vd.ghost_get<type, rho, rho_prev, Pressure, drho, force, velocity, velocity_prev, Background_P, v_transport>();
	// v_cl.barrier();
	vd.ghost_get<type, rho, Pressure, velocity, velocity_prev, v_transport>();
	// v_cl.barrier();

	if (BC_TYPE == NO_SLIP)
	{
		calc_boundary(vd, NN);
		vd.ghost_get<type, rho, Pressure, velocity, velocity_prev, v_transport>();
	}
	calc_forces(vd, NN);
	// v_cl.barrier();

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

		// vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);

		vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + dt_2 * vd.template getProp<force>(a)[0];
		vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + dt_2 * vd.template getProp<force>(a)[1];
		vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + dt_2 * vd.template getProp<force>(a)[2];

		vd.template getProp<v_transport>(a)[0] = vd.template getProp<velocity>(a)[0] + dt_2 * vd.template getProp<Background_P>(a)[0];
		vd.template getProp<v_transport>(a)[1] = vd.template getProp<velocity>(a)[1] + dt_2 * vd.template getProp<Background_P>(a)[1];
		vd.template getProp<v_transport>(a)[2] = vd.template getProp<velocity>(a)[2] + dt_2 * vd.template getProp<Background_P>(a)[2];

		++part2;
	}

	vd.ghost_get<type, rho, Pressure, velocity, velocity_prev, v_transport>();

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

	filename += ("_" + custom_string);
}

int main(int argc, char *argv[])
{

	// initialize the library
	openfpm_init(&argc, &argv);
	SetFilename(filename, "test");
	Box<3, double> domain;
	Box<3, double> fluid_box;
	Box<3, double> recipient;
	Box<3, double> recipient_hole;

	// Physical size of the fluid domain, it goes from (0,0,0) to (length[0],length[1],length[2])
	// First particle will always be placed at (dp/2,dp/2,dp/2) and the last particle will be placed at (length[0]-dp/2,length[1]-dp/2,length[2]-dp/2)
	double length[3];

	// Size of the virtual grid that defines where to place the particles
	size_t sz[3];

	// Here we define the boundary conditions of our problem
	size_t bc[3];

	// In the case of the new bc we need particles at the wall, for this we need sz_aux
	// We want to put one virtual grid point between each pair of the old ones,
	// so that the new spacing is dp/2, and we can put a fluid particle exactly at the wall
	size_t sz_aux[3];
	size_t Np_boundary[3];
	if (SCENARIO == POISEUILLE || SCENARIO == COUETTE || SCENARIO == HYDROSTATIC)
	{
		if (SCENARIO == POISEUILLE)
		{
			gravity_vector = {0.0, 0.0, 0.1};
			gravity = sqrt(norm2(gravity_vector));
			vw_top = {0.0, 0.0, 0.0};
			vw_bottom = {0.0, 0.0, 0.0};
			bc[0] = NON_PERIODIC;
			bc[1] = NON_PERIODIC;
			bc[2] = PERIODIC;
			size_t Nbound_x = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
			Np_boundary[0] = Nbound_x;
			Np_boundary[1] = 0;
			Np_boundary[2] = 0;
		}
		else if (SCENARIO == COUETTE)
		{
			gravity_vector = {0.0, 0.0, 0.0};
			gravity = sqrt(norm2(gravity_vector));
			vw_top = {0.0, 0.0, 1.25};
			vw_bottom = {0.0, 0.0, 0.0};
			bc[0] = NON_PERIODIC;
			bc[1] = NON_PERIODIC;
			bc[2] = PERIODIC;
			size_t Nbound_x = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
			Np_boundary[0] = Nbound_x;
			Np_boundary[1] = 0;
			Np_boundary[2] = 0;
		}
		else if (SCENARIO == HYDROSTATIC)
		{
			gravity_vector = {0.0, 0.0, -0.1};
			gravity = sqrt(norm2(gravity_vector));
			vw_top = {0.0, 0.0, 0.0};
			vw_bottom = {0.0, 0.0, 0.0};
			bc[0] = NON_PERIODIC;
			bc[1] = NON_PERIODIC;
			bc[2] = NON_PERIODIC;
			size_t Nbound_x = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
			Np_boundary[0] = Nbound_x;
			Np_boundary[1] = 0;
			Np_boundary[2] = Nbound_x;
		}

		// Number of fluid and boundary particles in the x, y and z direction
		size_t Np_fluid[3] = {40, 1, 20};

		// bc[0] = NON_PERIODIC;
		// bc[1] = NON_PERIODIC;
		// bc[2] = PERIODIC;

		// offset to add to the domain box to create the correct particle positions, if non periodic it has to go from -dp/2 - Np_boundary*dp to Lx+dp/2 + Np_boundary*dp
		double offset_domain[3];
		// offset to add to the recipient box, if non periodic it has to go from - Np_boundary*dp  to length + Np_boundary*dp
		double offset_recipient[3];
		// In case of periodic boundary conditions, we need to add an asymetric offset to the right,top or front of the domain
		double offset_periodic_fluid[3];
		double offset_periodic_recipient[3];

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
			length[dim] = dp * Np_fluid[dim];
			// non periodic, fluid covered by boundary
			if (bc[dim] == NON_PERIODIC)
			{
				sz[dim] = Np_fluid[dim] + 2 * (Np_boundary[dim] + 1);
				offset_domain[dim] = (0.5 + Np_boundary[dim]) * dp;
				offset_periodic_fluid[dim] = 0.0;
				offset_periodic_recipient[dim] = 0.0;

				if (Np_boundary[dim] != 0)
					sz_aux[dim] = 2 * sz[dim] - 1;
				else // for a direction with no boundary particles we dont need to add anything
					sz_aux[dim] = sz[dim];

				if (BC_TYPE == NEW_NO_SLIP) // Np_boundary should only be 0 or 1 if we are using the new bc
					offset_recipient[dim] = 0.25 * Np_boundary[dim] * dp;
				else
					offset_recipient[dim] = Np_boundary[dim] * dp;
			}
			// periodic, open ended
			else
			{
				sz[dim] = Np_fluid[dim] + 2;
				sz_aux[dim] = sz[dim];
				offset_domain[dim] = 0.5 * dp;
				offset_recipient[dim] = 0.0;
				offset_periodic_fluid[dim] = 0.75 * dp;
				offset_periodic_recipient[dim] = 0.85 * dp;
			}
		}

		// Define the boxes
		Box<3, double> box_domain({-offset_domain[0],
								   -offset_domain[1],
								   -offset_domain[2]},
								  {length[0] + offset_domain[0],
								   length[1] + offset_domain[1],
								   length[2] + offset_domain[2]});

		Box<3, double> box_fluid({0.0,
								  0.0,
								  0.0},
								 {length[0] + offset_periodic_fluid[0],
								  length[1] + offset_periodic_fluid[1],
								  length[2] + offset_periodic_fluid[2]});

		// Box<3, double> box_recipient({-offset_recipient[0],
		// 							  -offset_recipient[1],
		// 							  -offset_recipient[2]},
		// 							 {length[0] + offset_recipient[0] + offset_periodic_recipient[0],
		// 							  length[1] + offset_recipient[1] + offset_periodic_recipient[1],
		// 							  length[2] + offset_recipient[2] + offset_periodic_recipient[2]});
		double offset_recipient_top;
		if (SCENARIO == HYDROSTATIC)
		{
			offset_recipient_top = offset_recipient[2]; // 0.0;
		}
		else
		{
			offset_recipient_top = offset_recipient[2];
		}

		Box<3, double> box_recipient({-offset_recipient[0],
									  -offset_recipient[1],
									  -offset_recipient[2]},
									 {length[0] + offset_recipient[0] + offset_periodic_recipient[0],
									  length[1] + offset_recipient[1] + offset_periodic_recipient[1],
									  length[2] + offset_recipient_top + offset_periodic_recipient[2]});

		// Will only be used in the new bc
		Box<3, double> box_recipient_hole({offset_recipient[0],
										   offset_recipient[1],
										   offset_recipient[2]},
										  {length[0] - offset_recipient[0] + offset_periodic_fluid[0],
										   length[1] - offset_recipient[1] + offset_periodic_fluid[1],
										   length[2] - offset_recipient[2] + offset_periodic_fluid[2]});

		domain = box_domain;
		fluid_box = box_fluid;
		recipient = box_recipient;
		recipient_hole = box_recipient_hole;
	}

	// last position of fluid particles in z coordinate is at
	// length[2] + 0.5*dp
	double z_end = length[2] + 0.5 * dp;

	// Fill W_dap
	W_dap = 1.0 / Wab(dp);

	// extended boundary around the domain, and the processor domain
	Ghost<3, double> g(r_threshold + H);

	particles vd(0, domain, bc, g, DEC_GRAN(512));

	// return an iterator to the fluid particles to add to vd
	auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

	// Fill needed constants
	double L = length[0];
	double umax;

	if (SCENARIO == POISEUILLE)
	{
		umax = gravity_vector.get(2) * L * L / (8.0 * nu);
	}
	else if (SCENARIO == COUETTE)
	{
		umax = vw_top.get(2);
	}
	else if (SCENARIO == HYDROSTATIC)
	{
		umax = 1;
	}

	cbar = coeff_sound * umax;
	B = cbar * cbar * rho_zero / gamma_;
	Pbackground = B * Bfactor;
	double Re = L * umax / nu;

	std::string constants_filename = filename + "_Constants_" + ".txt";
	std::ofstream file(constants_filename);
	file << "gravity_x:  " << gravity_vector.get(0) << std::endl;
	file << "gravity_y:  " << gravity_vector.get(1) << std::endl;
	file << "gravity_z:  " << gravity_vector.get(2) << std::endl;
	file << "Re: " << Re << std::endl;
	file << "L: " << L << std::endl;
	file << "nu: " << nu << std::endl;
	file << "umax: " << umax << std::endl;
	file << "alpha: " << visco << std::endl;
	file << "cbar: " << cbar << std::endl;
	file << "B: " << B << std::endl;
	file << "Pbackground: " << Pbackground << std::endl;
	file << "rho_zero: " << rho_zero << std::endl;
	file << "gamma: " << gamma_ << std::endl;
	file << "dp: " << dp << std::endl;
	file << "xi: " << xi << std::endl;
	file << "H: " << H << std::endl;
	file << "H/dp: " << H / dp << std::endl;
	file << "CFLnumber: " << CFLnumber << std::endl;
	file << "BC_TYPE: " << BC_TYPE << std::endl;
	file << "SCENARIO: " << SCENARIO << std::endl;
	file << "MassOld: " << MassOld << std::endl;
	file << "MassFluid_OLD: " << MassFluid << std::endl;
	file << "MassBoundary_OLD: " << MassBound << std::endl;

	// for each particle inside the fluid box ...
	const double profile_parameter = gravity_vector.get(2) / (2.0 * nu);
	while (fluid_it.isNext())
	{

		// ... add a particle ...
		vd.add();

		// ... and set it position ...
		vd.getLastPos()[0] = fluid_it.get().get(0);
		vd.getLastPos()[1] = fluid_it.get().get(1);

		if (fluid_it.get().get(2) == z_end)
		{
			vd.getLastPos()[2] = fluid_it.get().get(2) - 1e-10 * dp;
		}
		else
		{
			vd.getLastPos()[2] = fluid_it.get().get(2);
		}

		// and its type.
		vd.template getLastProp<type>() = FLUID;

		// We also initialize the density of the particle and the hydro-static pressure given by
		//
		// rho_zero*g*h = P
		//
		// rho_p = (P/B + 1)^(1/Gamma) * rho_zero
		//

		vd.template getLastProp<Pressure>() = 0.0; // rho_zero * gravity_vector.get(2) * fluid_it.get().get(2) - rho_zero * gravity_vector.get(2) * length[2] / 2.0;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<rho_prev>() = rho_zero;
		vd.template getLastProp<velocity>()[0] = initial_perturbation * (-1.0 + (double)2.0 * rand() / RAND_MAX) * (L - std::abs(fluid_it.get().get(0) - L / 2.0));
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = initial_perturbation * (-1.0 + (double)2.0 * rand() / RAND_MAX) * (L - std::abs(fluid_it.get().get(0) - L / 2.0));
		// profile_parameter * fluid_it.get().get(0) * (L - fluid_it.get().get(0));

		// initial_perturbation * (-1.0 + (double)2.0 * rand() / RAND_MAX) * (L - std::abs(fluid_it.get().get(0) - L / 2.0));

		// profile_parameter * fluid_it.get().get(0) * (L - fluid_it.get().get(0));

		vd.template getLastProp<velocity_prev>()[0] = 0.0;
		vd.template getLastProp<velocity_prev>()[1] = 0.0;
		vd.template getLastProp<velocity_prev>()[2] = 0.0;

		vd.template getLastProp<force>()[0] = 0.0;
		vd.template getLastProp<force>()[1] = 0.0;
		vd.template getLastProp<force>()[2] = 0.0;
		// profile_parameter *fluid_it.get().get(0) * (L - fluid_it.get().get(0));

		// next fluid particle
		++fluid_it;
	}

	// Recipient

	openfpm::vector<Box<3, double>> holes;

	if (BC_TYPE == NEW_NO_SLIP)
	{
		holes.add(recipient_hole);
		sz[0] = sz_aux[0];
		sz[1] = sz_aux[1];
		sz[2] = sz_aux[2];
	}
	else
	{
		holes.add(fluid_box);
	}
	auto bound_box = DrawParticles::DrawSkin(vd, sz, domain, holes, recipient);
	double epsilon = 0.001;

	while (bound_box.isNext())
	{
		Point<3, double> position = bound_box.get();

		// periodic bc, with no boundary particles in y direction has a bug, it puts 3 extra particles outside in the y direction
		// When running on multiple cores, with this we check if particle is outside the recipient box
		// Another bug places boundary particles in the correct plane, but inside the fluid box;
		if ((!recipient.isInside((position))))
		{
			++bound_box;
			continue;
		}

		// if ((fluid_box.isInside(position) && !(position.get(2) < z_end + epsilon && position.get(2) > z_end - epsilon)))
		// {
		// 	++bound_box;
		// 	continue;
		// }

		vd.add();

		vd.getLastPos()[0] = bound_box.get().get(0);
		vd.getLastPos()[1] = bound_box.get().get(1);
		vd.getLastPos()[2] = bound_box.get().get(2);

		vd.template getLastProp<type>() = BOUNDARY;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<rho_prev>() = rho_zero;
		vd.template getLastProp<Pressure>() = 0.0; // rho_zero * gravity_vector.get(2) * bound_box.get().get(2) - rho_zero * gravity_vector.get(2) * length[2] / 2.0;

		if (bound_box.get().get(0) < 0) // bottom wall
		{
			vd.template getLastProp<velocity>()[0] = vw_bottom.get(0);
			vd.template getLastProp<velocity>()[1] = vw_bottom.get(1);
			vd.template getLastProp<velocity>()[2] = vw_bottom.get(2);

			vd.template getLastProp<force>()[0] = 0.0;
			vd.template getLastProp<force>()[1] = 0.0;
			vd.template getLastProp<force>()[2] = 0.0;
		}
		else if (bound_box.get().get(0) > length[0]) // top wall
		{
			vd.template getLastProp<velocity>()[0] = vw_top.get(0);
			vd.template getLastProp<velocity>()[1] = vw_top.get(1);
			vd.template getLastProp<velocity>()[2] = vw_top.get(2);
		}

		vd.template getLastProp<velocity_prev>()[0] = 0.0;
		vd.template getLastProp<velocity_prev>()[1] = 0.0;
		vd.template getLastProp<velocity_prev>()[2] = 0.0;

		++bound_box;
		// std::cout << "Boundary particle " << count << " at position x=" << vd.getLastPos()[0] << " y=" << vd.getLastPos()[1] << " z=" << vd.getLastPos()[2] << std::endl;
	}
	openfpm::vector<std::string> names({"type", "rho", "rho_prev", "pressure", "drho", "force", "velocity", "velocity_prev", "Background_P", "v_transport"});
	vd.setPropNames(names);

	vd.map();

	// vd.write_frame("FIRST", 0, WRITER);

	// Now that we fill the vector with particles
	ModelCustom md;

	vd.addComputationCosts(md);
	vd.getDecomposition().decompose();
	vd.map();

	vd.ghost_get<type, rho, Pressure, velocity, velocity_prev, v_transport>();

	auto NN = vd.getCellList(r_threshold + H);

	// Adjust mass to have the correct density
	fix_mass(vd, NN);
	file << "MassFluid: " << MassFluid << std::endl;
	file << "MassBoundary: " << MassBound << std::endl;

	// Evolve
	size_t write = 0;
	size_t it = 0;
	size_t it_reb = 0;
	double t = 0.0;

	// Compute forces for the first iteration

	vd.ghost_get<type, rho, Pressure, velocity, velocity_prev, v_transport>();

	if (BC_TYPE == NO_SLIP)
	{
		calc_boundary(vd, NN);
	}
	vd.ghost_get<type, rho, Pressure, velocity, velocity_prev, v_transport>();
	calc_forces(vd, NN);
	// sync all to ghost
	vd.ghost_get<type, rho, drho, Pressure, force, velocity, velocity_prev, Background_P, v_transport>();

	while (t <= t_end)
	{
		Vcluster<> &v_cl = create_vcluster();
		timer it_time;

		//// Do rebalancing every 200 timesteps
		it_reb++;
		if (it_reb == 200)
		{
			vd.map();

			it_reb = 0;
			ModelCustom md;
			vd.addComputationCosts(md);
			vd.getDecomposition().decompose();

			if (v_cl.getProcessUnitID() == 0)
				std::cout << "REBALANCED " << std::endl;
		}
		// vd.map();

		// Calculate dt for time stepping
		double dt = calc_deltaT(vd);

		// Integrate one time step
		kick_drift_int(vd, NN, dt, v_cl);

		// increment time
		t += dt;

		if (write < t * write_const)
		{
			// vd.map();

			// vd.deleteGhost();
			vd.write_frame(filename, write, WRITER);
			// vd.ghost_get<type, rho, Pressure, velocity, velocity_prev, v_transport>();

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
