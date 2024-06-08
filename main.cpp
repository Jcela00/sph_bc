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
#define FREE_SLIP 1
#define NEW_NO_SLIP 2

// Type of viscosity
#define ARTIFICIAL_VISCOSITY 0
#define PHYSICAL_VISCOSITY 1

// Type of pressure force
#define OLD_PRESSURE 0
#define NEW_PRESSURE 1
#define NO_PRESSURE 2

// TYPE OF SCENARIO
#define POISEUILLE 0
#define COUETTE 1

// TENSILE ENABLED
#define TENSILE_ENABLED 1
#define TENSILE_DISABLED 0

//////// CODE PARAMETERS /////////////////////////////////////////////////////////////
// Dimensionality of the problem changes the normalizaiton of the kernel
const int dimensions = 2;

// BC and viscosity we are using in the simulation
const int BC_TYPE = NO_SLIP;
const int VISC_TYPE = PHYSICAL_VISCOSITY;
const int PRESSURE_TYPE = NEW_PRESSURE;
const int SCENARIO = POISEUILLE;
const int WRITER = VTK_WRITER; // VTK_WRITER or CSV_WRITER
const int TENSILE = TENSILE_ENABLED;

// Output file name

std::string filename = "test";

// std::string filename = "Poiseuille_NoP_ArtVisc"; // ALREADY USED
// std::string filename = "Poiseuille_OldP_ArtVisc"; // AlREADY USED
// std::string filename = "Poiseuille_NewP_ArtVisc"; // ALREADY USED
// std::string filename = "Poiseuille_NoP_PhysVisc"; // ALREADY USED
// std::string filename = "Poiseuille_OldP_PhysVisc"; // ALREADY USED
// std::string filename = "Poiseuille_NewP_PhysVisc"; // ALREADY USED

// std::string filename = "Couette_NoP_ArtVisc"; // ALREADY USED
// std::string filename = "Couette_NewP_ArtVisc"; // ALREADY USED
// std::string filename = "Couette_OldP_ArtVisc"; // AlREADY USED
// Controls otput file frequency, low means less frequent
const int write_const = 10;
////////////////////////////////////////////////////////////////////////////////////

//////// SPATIAL CONSTANTS /////////////////////////////////////////////////////////
// initial spacing between particles dp in the formulas
const double dp = 0.025;

// Factor relating H (smoothing length) and dp ( particle spacing)
const double Hconst = sqrt(3.0);

// H, smoothing length
const double H = Hconst * dp;
////////////////////////////////////////////////////////////////////////////////////

//////// EQ OF STATE CONSTANTS //////////////////////////////////////////////////////
// Eq of state constant, filled later
double B = 0.0;

// Constant used for the sound speed, number of times the max velocity
const double coeff_sound = 10.0;

// Sound speed (calculated later)
double cbar = 0.0;

// gamma in the formulas
const double gamma_ = 7.0;
const double xi = 0.0;

////////////////////////////////////////////////////////////////////////////////////

//////// PHYSICAL CONSTANTS ////////////////////////////////////////////////////////
// Reference density
const double rho_zero = 1.0;

// Minimum, Maximum Rho allowed
const double RhoMin = 0.7 * rho_zero;
const double RhoMax = 1.3 * rho_zero;

// Mass of the fluid and boundary particles M=rho*V, V=dp^3 or V = dp^2
const double MassFluid = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);
const double MassBound = rho_zero * (dimensions == 3 ? dp * dp * dp : dp * dp);

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
const double nu = 0.01;

// Viscosity coefficient of the artificial viscosity (alpha) filled later if using artificial viscosity
double visco = 500;

// Dynamic viscosity, used when using physical viscosity term
double eta = nu * rho_zero;

// Eta in the formulas
const double Eta2 = 0.01 * H * H;

const double initial_perturbation = 0.0;

double Pbackground = 10; // filled later

////////////////////////////////////////////////////////////////////////////////////

//////// TIME CONSTANTS //////////////////////////////////////////////////////////////
// End simulation time
const double t_end = 100;
// Constant used to define time integration
const double CFLnumber = 0.1;
// Minimum T
const double DtMin = 0.00001;
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

typedef vector_dist<3, double, aggregate<size_t, double, double, double, double, double[3], double[3], double[3], double[3]>> particles;
//                                       |      |        |          |            |            |         |            |		   |
//                                       |      |        |          |            |            |         |            |		   |
//                                     type   density   density    Pressure    delta       force     velocity    velocity	background pressure force
//                                                      at n-1                 density                           at n - 1

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
// Particle wise, inverted equation of state, ie compute density given pressure
inline double InvEqState(const double pressure)
{
	return rho_zero * std::pow(((pressure - xi) / B + 1.0), 1.0 / gamma_);
}

// Normalization constant for the kernel
// In 3D the normalization is 1/(pi*H^3) in 2d is 10/(7*pi*H^2)
const double a2 = (dimensions == 3) ? 1.0 / (M_PI * H * H * H) : 10.0 / (7.0 * M_PI * H * H);

// Kernel function
inline double Wab(double r)
{
	r /= H;
	if (r < 1.0)
		return (1.0 - (3.0 / 2.0) * r * r + (3.0 / 4.0) * r * r * r) * a2;
	else if (r < 2.0)
	{
		return ((1.0 / 4.0) * (2.0 - r) * (2.0 - r) * (2.0 - r)) * a2;
	}
	else
		return 0.0;
}

inline void DWab(const Point<3, double> &dx, Point<3, double> &DW, const double r)
{
	const double c1 = (-3.0 / H) * a2;
	const double d1 = (9.0 / 4.0 / H) * a2;
	const double c2 = (-3.0 / 4.0 / H) * a2;
	const double qq = r / H;

	double qq2 = qq * qq;
	double fac1 = (c1 * qq + d1 * qq2) / r;
	double b1 = (qq < 1.0) ? 1.0 : 0.0;

	double wqq = (2.0 - qq);
	double fac2 = c2 * wqq * wqq / r;
	double b2 = (qq >= 1.0 && qq < 2.0) ? 1.0 : 0.0;

	double factor = (b1 * fac1 + b2 * fac2);

	DW.get(0) = factor * dx.get(0);
	DW.get(1) = factor * dx.get(1);
	DW.get(2) = factor * dx.get(2);
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

	if (TENSILE == TENSILE_ENABLED)
		return (fab * (tensilp1 + tensilp2));
	else
		return 0.0;
}

inline double Pi_artificial(const Point<3, double> &dr, const double rr2, const Point<3, double> &dv, const double rhoa, const double rhob, const double massb, double &visc)
{
	const double dot = dr.get(0) * dv.get(0) + dr.get(1) * dv.get(1) + dr.get(2) * dv.get(2);
	const double dot_rr2 = dot / (rr2 + Eta2);
	visc = std::max(dot_rr2, visc);

	const float amubar = H * dot_rr2;
	const float robar = (rhoa + rhob) * 0.5;
	const float pi_visc = (-visco * cbar * amubar / robar);

	return pi_visc;
}

inline Point<3, double> Pi_physical(const Point<3, double> &dr, const double r, const Point<3, double> &dv, const Point<3, double> &dW, const double rhoa, const double rhob, const double massa, const double massb, double &visc)
{
	const Point<3, double> e_ab = dr / r; // unit vector from a to b
	const double dot = e_ab.get(0) * dW.get(0) + e_ab.get(1) * dW.get(1) + e_ab.get(2) * dW.get(2);

	double Va2 = (massa / rhoa) * (massa / rhoa);
	double Vb2 = (massb / rhob) * (massb / rhob);
	Point<3, double> pi_visc = eta * ((Va2 + Vb2) * dot * dv) / r;
	const double normvisc = sqrt(norm2(pi_visc));
	visc = std::max(normvisc, visc);
	return pi_visc;
}

inline double PressureForce(const double rhoa, const double rhob, const double prsa, const double prsb, const double massa, const double massb)
{
	return -massb * (prsa + prsb) / (rhoa * rhob);
}
inline double PressureForceNew(const double rhoa, const double rhob, const double prsa, const double prsb, const double massa, const double massb)
{
	const double Va2 = (massa / rhoa) * (massa / rhoa);
	const double Vb2 = (massb / rhob) * (massb / rhob);
	const double pbar = (rhob * prsa + rhoa * prsb) / (rhoa + rhob);

	return (-1.0 / massa) * (Va2 + Vb2) * pbar;
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
				if (r2 < 4.0 * H * H)
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
					const double dot = gravity_vector.get(0) * dr.get(0) + gravity_vector.get(1) * dr.get(1) + gravity_vector.get(2) * dr.get(2);
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

// template <typename CellList>
void interact_fluid_boundary_new(particles &vd,
								 vect_dist_key_dx fluid_key,
								 const double massf,
								 const double rhof,
								 const double Pf,
								 const Point<3, double> xf,
								 const Point<3, double> vf,
								 const Point<3, double> &xw,
								 const Point<3, double> &r_fw,
								 const double &dist2,
								 unsigned long &boundary_key,
								 double &max_visc)
{

	const double massw = MassBound;
	const double rhow = vd.getProp<rho>(boundary_key);
	const double Pw = vd.getProp<Pressure>(boundary_key);
	const Point<3, double> vw = vd.getProp<velocity>(boundary_key);

	// xw.get(0) this is the x coordinate of the wall
	bool is_bottom_wall = (xw.get(0) < dp ? true : false);
	Point<3, double> normal = (is_bottom_wall ? normal_bottom_wall : -1.0 * normal_bottom_wall);
	Point<3, double> tangential = {normal.get(2), 0.0, normal.get(0)};
	if (!is_bottom_wall)
	{
		tangential.get(2) = -tangential.get(2);
	}
	Point<3, double> r_wf = -1.0 * r_fw; // Points to wall from fluid

	// Project va on tangential and normal directions
	double vt = (vf.get(0) * tangential.get(0) + vf.get(1) * tangential.get(1) + vf.get(2) * tangential.get(2));
	double vn = (vf.get(0) * normal.get(0) + vf.get(1) * normal.get(1) + vf.get(2) * normal.get(2));

	double lf = r_wf.get(0) * normal.get(0) + r_wf.get(1) * normal.get(1) + r_wf.get(2) * normal.get(2); // vertical distance from fluid particle to wall
	lf = (lf < 0.0 ? -1.0 * lf : lf);																	 // absolute value

	Point<3, double> offset_1 = -0.5 * normal * dp; // offset from wall to first particle
	Point<3, double> offset_2 = -1.0 * normal * dp; // offset from first to second particle, and from second to third

	Point<3, double> r1 = r_wf + offset_1; // Vector from fluid particle to first dummy particle
	double r1_norm = sqrt(norm2(r1));	   // Norm of r1
	double lw1 = dp / 2.0;				   // Distance from wall to first dummy particle

	Point<3, double> r2 = r1 + offset_2; // Vector from fluid particle to second dummy particle
	double r2_norm = sqrt(norm2(r2));	 // Norm of r2
	double lw2 = 1.5 * dp;				 // Distance from wall to second dummy particle

	Point<3, double> r3 = r2 + offset_2; // Vector from fluid particle to third dummy particle
	double r3_norm = sqrt(norm2(r3));	 // Norm of r3
	double lw3 = 2.5 * dp;				 // Distance from wall to third dummy particle

	Point<3, double> v1 = {0.0, 0.0, 0.0}; // Velocity of first dummy particle
	Point<3, double> v2 = {0.0, 0.0, 0.0}; // Velocity of second dummy particle
	Point<3, double> v3 = {0.0, 0.0, 0.0}; // Velocity of third dummy particle
	double p1, p2, p3;					   // Pressure of the dummy particles
	double rho1, rho2, rho3;			   // Density of the dummy particles

	std::vector<Point<3, double>> r_dummy;
	std::vector<Point<3, double>> v_dummy;
	std::vector<double> P_dummy;
	std::vector<double> rho_dummy;
	std::vector<Point<3, double>> W_dummy;
	if (r1_norm < 2.0 * H)
	{
		v1 = 2.0 * vw - vt * (lw1 / lf) * tangential + vn * (lw1 / lf) * normal;
		double dot = gravity_vector.get(0) * normal.get(0) + gravity_vector.get(1) * normal.get(1) + gravity_vector.get(2) * normal.get(2);
		p1 = Pf + rhof * dot * (lf + lw1);
		rho1 = InvEqState(p1);
		r1 = -1.0 * r1; // force routines use r pointing from wall to fluid
		Point<3, double> DW;
		DWab(r1, DW, r1_norm);

		r_dummy.push_back(r1);
		v_dummy.push_back(v1);
		P_dummy.push_back(p1);
		rho_dummy.push_back(rho1);
		W_dummy.push_back(Wab(r1_norm));
	}
	else if (r2_norm < 2.0 * H)
	{
		v2 = 2.0 * vw - vt * (lw2 / lf) * tangential + vn * (lw2 / lf) * normal;
		double dot = gravity_vector.get(0) * normal.get(0) + gravity_vector.get(1) * normal.get(1) + gravity_vector.get(2) * normal.get(2);
		p2 = Pf + rhof * dot * (lf + lw2);
		rho2 = InvEqState(p2);
		r2 = -1.0 * r2;
		Point<3, double> DW;
		DWab(r2, DW, r2_norm);

		r_dummy.push_back(r2);
		v_dummy.push_back(v2);
		P_dummy.push_back(p2);
		rho_dummy.push_back(rho2);
		W_dummy.push_back(Wab(r2_norm));
	}
	else if (r3_norm < 2.0 * H)
	{
		v3 = 2.0 * vw - vt * (lw3 / lf) * tangential + vn * (lw3 / lf) * normal;
		double dot = gravity_vector.get(0) * normal.get(0) + gravity_vector.get(1) * normal.get(1) + gravity_vector.get(2) * normal.get(2);
		p3 = Pf + rhof * dot * (lf + lw3);
		rho3 = InvEqState(p3);
		r3 = -1.0 * r3;

		Point<3, double> DW;
		DWab(r3, DW, r3_norm);

		r_dummy.push_back(r3);
		v_dummy.push_back(v3);
		P_dummy.push_back(p3);
		rho_dummy.push_back(rho3);
		W_dummy.push_back(Wab(r3_norm));
	}

	for (size_t comp = 0; comp < v_dummy.size(); comp++)
	{

		Point<3, double> PhysicalViscosityTerm;
		double ArtificalViscosityTerm;
		double PressureTerm;

		double TensileTerm = -MassBound * Tensile(sqrt(norm2(r_dummy[comp])), rhof, rho_dummy[comp], Pf, P_dummy[comp]);

		// Compute pressure term depending on the type of pressure
		if (PRESSURE_TYPE == OLD_PRESSURE)
			PressureTerm = PressureForce(rhof, rho_dummy[comp], Pf, P_dummy[comp], massf, MassBound);
		else if (PRESSURE_TYPE == NEW_PRESSURE)
			PressureTerm = PressureForceNew(rhof, rho_dummy[comp], Pf, P_dummy[comp], massf, MassBound);
		else if (PRESSURE_TYPE == NO_PRESSURE)
			PressureTerm = 0.0;

		double factor = PressureTerm + TensileTerm;

		if (VISC_TYPE == PHYSICAL_VISCOSITY)
		{
			PhysicalViscosityTerm = Pi_physical(r_dummy[comp], sqrt(norm2(r_dummy[comp])), v_dummy[comp], W_dummy[comp], rhof, rho_dummy[comp], massf, MassBound, max_visc) / massf;
			vd.getProp<force>(fluid_key)[0] += PhysicalViscosityTerm.get(0);
			vd.getProp<force>(fluid_key)[1] += PhysicalViscosityTerm.get(1);
			vd.getProp<force>(fluid_key)[2] += PhysicalViscosityTerm.get(2);
		}
		else if (VISC_TYPE == ARTIFICIAL_VISCOSITY)
		{
			ArtificalViscosityTerm = -MassBound * Pi_artificial(r_dummy[comp], sqrt(norm2(r_dummy[comp])), v_dummy[comp], rhof, rho_dummy[comp], massf, max_visc);
			factor += ArtificalViscosityTerm;
		}
		vd.getProp<force>(fluid_key)[0] += factor * (W_dummy[comp]).get(0);
		vd.getProp<force>(fluid_key)[1] += factor * (W_dummy[comp]).get(1);
		vd.getProp<force>(fluid_key)[2] += factor * (W_dummy[comp]).get(2);

		vd.getProp<drho>(fluid_key) += MassBound * (vf.get(0) * (W_dummy[comp]).get(0) + vf.get(1) * (W_dummy[comp]).get(1) + vf.get(2) * (W_dummy[comp]).get(2));
	}

	// if (is_bottom_wall)
	// {
	// 	std::cout << "Bottom wall" << std::endl;
	// 	std::cout << "normal: " << normal.get(0) << " " << normal.get(1) << " " << normal.get(2) << std::endl;
	// 	std::cout << "tangential: " << tangential.get(0) << " " << tangential.get(1) << " " << tangential.get(2) << std::endl;
	// 	std::cout << "xw: " << xw.get(0) << " " << xw.get(1) << " " << xw.get(2) << std::endl;
	// 	std::cout << "r_fw" << r_fw.get(0) << " " << r_fw.get(1) << " " << r_fw.get(2) << std::endl;
	// 	std::cout << "offset_1: " << offset_1.get(0) << " " << offset_1.get(1) << " " << offset_1.get(2) << std::endl;
	// 	std::cout << "offset_2: " << offset_2.get(0) << " " << offset_2.get(1) << " " << offset_2.get(2) << std::endl;
	// 	std::cout << "r1: " << r1.get(0) / dp << " " << r1.get(1) / dp << " " << r1.get(2) / dp << " norm: " << r1_norm << std::endl;
	// 	std::cout << "r2: " << r2.get(0) / dp << " " << r2.get(1) / dp << " " << r2.get(2) / dp << " norm: " << r2_norm << std::endl;
	// 	std::cout << "r3: " << r3.get(0) / dp << " " << r3.get(1) / dp << " " << r3.get(2) / dp << " norm: " << r3_norm << std::endl;
	// 	std::cout << "vt: " << vt << std::endl;
	// 	std::cout << "vn: " << vn << std::endl;
	// 	std::cout << "v1: " << v1.get(0) << " " << v1.get(1) << " " << v1.get(2) << std::endl;
	// 	std::cout << "v2: " << v2.get(0) << " " << v2.get(1) << " " << v2.get(2) << std::endl;
	// 	std::cout << "v3: " << v3.get(0) << " " << v3.get(1) << " " << v3.get(2) << std::endl;

	// 	std::cout << std::endl;
	// }
	// else
	// {
	// 	std::cout << "Top wall" << std::endl;
	// 	std::cout << "normal: " << normal.get(0) << " " << normal.get(1) << " " << normal.get(2) << std::endl;
	// 	std::cout << "tangential: " << tangential.get(0) << " " << tangential.get(1) << " " << tangential.get(2) << std::endl;
	// 	std::cout << "xw: " << xw.get(0) << " " << xw.get(1) << " " << xw.get(2) << std::endl;
	// 	std::cout << "r_fw" << r_fw.get(0) << " " << r_fw.get(1) << " " << r_fw.get(2) << std::endl;
	// 	std::cout << "offset_1: " << offset_1.get(0) << " " << offset_1.get(1) << " " << offset_1.get(2) << std::endl;
	// 	std::cout << "offset_2: " << offset_2.get(0) << " " << offset_2.get(1) << " " << offset_2.get(2) << std::endl;
	// 	std::cout << "r1: " << r1.get(0) / dp << " " << r1.get(1) / dp << " " << r1.get(2) / dp << " norm: " << r1_norm << std::endl;
	// 	std::cout << "r2: " << r2.get(0) / dp << " " << r2.get(1) / dp << " " << r2.get(2) / dp << " norm: " << r2_norm << std::endl;
	// 	std::cout << "r3: " << r3.get(0) / dp << " " << r3.get(1) / dp << " " << r3.get(2) / dp << " norm: " << r3_norm << std::endl;
	// 	std::cout << "vt: " << vt << std::endl;
	// 	std::cout << "vn: " << vn << std::endl;
	// 	std::cout << "v1: " << v1.get(0) << " " << v1.get(1) << " " << v1.get(2) << std::endl;
	// 	std::cout << "v2: " << v2.get(0) << " " << v2.get(1) << " " << v2.get(2) << std::endl;
	// 	std::cout << "v3: " << v3.get(0) << " " << v3.get(1) << " " << v3.get(2) << std::endl;
	// 	std::cout << std::endl;
	// }
}

// template <typename CellList>
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
						  const unsigned long &b_key,
						  double &max_visc)
{

	const double massb = MassFluid;
	const double rhob = vd.getProp<rho>(b_key);
	const double Pb = vd.getProp<Pressure>(b_key);
	const Point<3, double> vb = vd.getProp<velocity>(b_key);

	const double r = sqrt(r2);

	const Point<3, double> v_rel = va - vb;

	Point<3, double> DW;
	DWab(r_ab, DW, r);

	Point<3, double> PhysicalViscosityTerm;
	double ArtificalViscosityTerm;
	double PressureTerm;
	double TensileTerm = -massb * Tensile(r, rhoa, rhob, Pa, Pb);

	if (PRESSURE_TYPE == OLD_PRESSURE)
		PressureTerm = PressureForce(rhoa, rhob, Pa, Pb, massa, massb);
	else if (PRESSURE_TYPE == NEW_PRESSURE)
		PressureTerm = PressureForceNew(rhoa, rhob, Pa, Pb, massa, massb);
	else if (PRESSURE_TYPE == NO_PRESSURE)
		PressureTerm = 0.0;

	double factor = PressureTerm + TensileTerm;

	if (VISC_TYPE == PHYSICAL_VISCOSITY)
	{
		PhysicalViscosityTerm = Pi_physical(r_ab, r, v_rel, DW, rhoa, rhob, massa, massb, max_visc) / massa;
		vd.getProp<force>(a_key)[0] += PhysicalViscosityTerm.get(0);
		vd.getProp<force>(a_key)[1] += PhysicalViscosityTerm.get(1);
		vd.getProp<force>(a_key)[2] += PhysicalViscosityTerm.get(2);
	}
	else if (VISC_TYPE == ARTIFICIAL_VISCOSITY)
	{
		ArtificalViscosityTerm = -massb * Pi_artificial(r_ab, r2, v_rel, rhoa, rhob, massa, max_visc);
		factor += ArtificalViscosityTerm;
	}

	vd.getProp<force>(a_key)[0] += factor * DW.get(0);
	vd.getProp<force>(a_key)[1] += factor * DW.get(1);
	vd.getProp<force>(a_key)[2] += factor * DW.get(2);

	double Va2 = (massa / rhoa) * (massa / rhoa);
	double Vb2 = (massb / rhob) * (massb / rhob);
	vd.getProp<Background_P>(a_key)[0] += (Pbackground / massa) * (Va2 + Vb2) * DW.get(0);
	vd.getProp<Background_P>(a_key)[1] += (Pbackground / massa) * (Va2 + Vb2) * DW.get(1);
	vd.getProp<Background_P>(a_key)[2] += (Pbackground / massa) * (Va2 + Vb2) * DW.get(2);

	vd.getProp<drho>(a_key) += rhoa * (massb / rhob) * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
}

// template <typename CellList>
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
								 unsigned long &boundary_key,
								 double &max_visc)
{
	const double massb = MassBound;
	const double rhob = vd.getProp<rho>(boundary_key);
	const double Pb = vd.getProp<Pressure>(boundary_key);
	const Point<3, double> vb = vd.getProp<velocity>(boundary_key);
	const Point<3, double> vb_noslip = vd.getProp<velocity_prev>(boundary_key);

	const Point<3, double> r_fb = xf - xb;
	const double r = sqrt(r2);

	const Point<3, double> v_rel = vf - vb;
	Point<3, double> v_rel_aux = vf - vb_noslip;

	Point<3, double> DW;
	DWab(r_fb, DW, r);

	Point<3, double> PhysicalViscosityTerm;
	double ArtificalViscosityTerm;
	double PressureTerm;
	double TensileTerm = -massb * Tensile(r, rhof, rhob, Pf, Pb);

	// Compute pressure term depending on the type of pressure
	if (PRESSURE_TYPE == OLD_PRESSURE)
		PressureTerm = PressureForce(rhof, rhob, Pf, Pb, massf, massb);
	else if (PRESSURE_TYPE == NEW_PRESSURE)
		PressureTerm = PressureForceNew(rhof, rhob, Pf, Pb, massf, massb);
	else if (PRESSURE_TYPE == NO_PRESSURE)
		PressureTerm = 0.0;

	double factor = PressureTerm + TensileTerm;

	if (VISC_TYPE == PHYSICAL_VISCOSITY)
	{
		PhysicalViscosityTerm = Pi_physical(r_fb, r, v_rel_aux, DW, rhof, rhob, massf, massb, max_visc) / massf;
		vd.getProp<force>(fluid_key)[0] += PhysicalViscosityTerm.get(0);
		vd.getProp<force>(fluid_key)[1] += PhysicalViscosityTerm.get(1);
		vd.getProp<force>(fluid_key)[2] += PhysicalViscosityTerm.get(2);
	}
	else if (VISC_TYPE == ARTIFICIAL_VISCOSITY)
	{
		ArtificalViscosityTerm = -massb * Pi_artificial(r_fb, r2, v_rel_aux, rhof, rhob, massf, max_visc);
		factor += ArtificalViscosityTerm;
	}

	vd.getProp<force>(fluid_key)[0] += factor * DW.get(0);
	vd.getProp<force>(fluid_key)[1] += factor * DW.get(1);
	vd.getProp<force>(fluid_key)[2] += factor * DW.get(2);

	double Va2 = (massf / rhof) * (massf / rhof);
	double Vb2 = (massb / rhob) * (massb / rhob);
	vd.getProp<Background_P>(fluid_key)[0] += (Pbackground / massf) * (Va2 + Vb2) * DW.get(0);
	vd.getProp<Background_P>(fluid_key)[1] += (Pbackground / massf) * (Va2 + Vb2) * DW.get(1);
	vd.getProp<Background_P>(fluid_key)[2] += (Pbackground / massf) * (Va2 + Vb2) * DW.get(2);

	vd.getProp<drho>(fluid_key) += rhof * (massb / rhob) * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
}

template <typename CellList>
inline void calc_forces(particles &vd, CellList &NN, double &max_visc)
{
	vector_dist_iterator part = vd.getDomainIterator();

	// Update the cell-list
	vd.updateCellList(NN);

	// For each particle ...
	while (part.isNext())
	{
		// get particle a
		vect_dist_key_dx a = part.get();

		// Get the position xp of the particle
		Point<3, double> xa = vd.getPos(a);

		// Take the mass of the particle dependently if it is FLUID or BOUNDARY
		double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;

		// Get the density and pressure of the of the particle a
		double rhoa = vd.getProp<rho>(a);
		double Pa = vd.getProp<Pressure>(a);

		// Get the Velocity of the particle a
		Point<3, double> va = vd.getProp<velocity>(a);

		// Reset the force counter (- gravity on zeta direction)
		vd.template getProp<force>(a)[0] = gravity_vector.get(0);
		vd.template getProp<force>(a)[1] = gravity_vector.get(1);
		vd.template getProp<force>(a)[2] = gravity_vector.get(2);

		vd.template getProp<Background_P>(a)[0] = 0.0;
		vd.template getProp<Background_P>(a)[1] = 0.0;
		vd.template getProp<Background_P>(a)[2] = 0.0;

		vd.template getProp<drho>(a) = 0.0;

		// We threat FLUID particle differently from BOUNDARY PARTICLES ...
		// INTERACTION OF A BOUNDARY PARTICLE WITH ITS NEIGHBORHOOD
		if (vd.getProp<type>(a) != FLUID)
		{
			// If it is a boundary particle calculate the delta rho based on equation 2
			// This require to run across the neighborhoods particles of a
			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

			// For each neighborhood particle
			while (Np.isNext() == true)
			{
				// get particle b key
				unsigned long b = Np.get();

				// if (a == b) skip this particle
				if (a.getKey() == b)
				{
					++Np;
					continue;
				};
				// if interacting with another boundary particle skip
				if (vd.getProp<type>(b) == BOUNDARY)
				{
					++Np;
					continue;
				}

				// Get the position xb of the particle
				Point<3, double> xb = vd.getPos(b);

				// Get the distance between a and b
				Point<3, double> dr = xa - xb;
				// take the norm squared of this vector
				double r2 = norm2(dr);

				// If the particles interact ...
				if (r2 < 4.0 * H * H)
				{

					// get the mass of the particle
					double massb = MassFluid;

					// Get the pressure and density of particle b
					double Pb = vd.getProp<Pressure>(b);
					double rhob = vd.getProp<rho>(b);

					// Get the velocity of the particle b
					Point<3, double> vb = vd.getProp<velocity>(b);

					// ... calculate distance
					double r = sqrt(r2);

					Point<3, double> dv = va - vb;

					Point<3, double> DW;
					DWab(dr, DW, r);

					if (VISC_TYPE == PHYSICAL_VISCOSITY)
					{
						// call the viscosity funciton just to update max_visc
						Point<3, double> dummy = Pi_physical(dr, r, dv, DW, rhoa, rhob, massa, massb, max_visc);
					}
					else if (VISC_TYPE == ARTIFICIAL_VISCOSITY)
					{
						const double dot = dr.get(0) * dv.get(0) + dr.get(1) * dv.get(1) + dr.get(2) * dv.get(2);
						const double dot_rr2 = dot / (r2 + Eta2);
						max_visc = std::max(dot_rr2, max_visc);
					}

					if (BC_TYPE == NO_SLIP || BC_TYPE == NEW_NO_SLIP)
					{
						// We dont evolve boundary particle density
						vd.getProp<drho>(a) = 0.0;
					}
					else if (BC_TYPE == FREE_SLIP)
					{
						vd.getProp<drho>(a) += massb * (dv.get(0) * DW.get(0) + dv.get(1) * DW.get(1) + dv.get(2) * DW.get(2));
					}
				}

				++Np;
			}
		}
		// INTERACTION OF A FLUID PARTICLE WITH ITS NEIGHBORHOOD
		else
		{
			// If it is a fluid particle calculate based on equation 1 and 2

			// Get an iterator over the neighborhood particles of p
			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

			// std::cout << "Type Np: " << type_name<decltype(Np)>() << '\n';

			// For each neighborhood particle
			while (Np.isNext() == true)
			{
				// get particle b
				unsigned long b = Np.get();

				// std::cout << "Type b: " << type_name<decltype(b)>() << '\n';

				// if (p == q) skip this particle
				if (a.getKey() == b)
				{
					++Np;
					continue;
				};
				// Get the position xp of the particle
				Point<3, double> xb = vd.getPos(b);

				// Get the distance between a and b
				// in fluid - boundary its xf-xb
				Point<3, double> dr = xa - xb;

				// take the norm (squared) of this vector
				double r2 = norm2(dr);

				// if they interact
				if (r2 < 4.0 * H * H)
				{
					if (vd.getProp<type>(b) == BOUNDARY)
					{
						if (BC_TYPE == NO_SLIP)
							interact_fluid_boundary_old(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b, max_visc);
						else
							interact_fluid_boundary_new(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b, max_visc);
					}
					else
					{
						interact_fluid_fluid(vd, a, massa, rhoa, Pa, xa, va, xb, dr, r2, b, max_visc);
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

double calc_deltaT(particles &vd, double ViscDtMax)
{
	double Maxacc = 0.0;
	double Maxvel = 0.0;
	max_acceleration_and_velocity(vd, Maxacc, Maxvel);

	double dt_u = 0.25 * H / (cbar + abs(Maxvel));
	double dt_visc = 0.25 * H * H / (nu);
	double dt_g = 0.25 * sqrt(H / gravity_vector.norm());
	double dt = std::min({dt_u, dt_visc, dt_g});

	if (dt < DtMin)
		dt = DtMin;

	return dt;
	// 	//-dt1 depends on force per unit mass.
	// 	const double dt_f = (Maxacc) ? sqrt(H / Maxacc) : std::numeric_limits<int>::max();

	// 	//-dt2 combines the Courant and the viscous time-step controls.
	// 	// const double dt_cv = H / (std::max(cbar, Maxvel * 10.) + H * ViscDtMax);

	// 	const double dt_cv = H / (std::max(cbar, Maxvel * 10.) + H * ViscDtMax);
	// 	//-dt new value of time step.
	// 	double dt = double(CFLnumber) * std::min(dt_f, dt_cv);
	// 	if (dt < double(DtMin))
	// 		dt = double(DtMin);

	// 	return dt;
}

size_t cnt = 0;

void verlet_int(particles &vd, double dt)
{

	// particle iterator
	auto part = vd.getDomainIterator();

	double dt205 = dt * dt * 0.5;
	double dt2 = dt * 2.0;

	// For each particle ...
	while (part.isNext())
	{
		// ... a
		auto a = part.get();

		// if the particle is boundary
		if (vd.template getProp<type>(a) == BOUNDARY)
		{

			if (BC_TYPE == FREE_SLIP)
			{
				// Update rho
				double rhop = vd.template getProp<rho>(a);
				// Update only the density
				double rhonew = vd.template getProp<rho_prev>(a) + dt2 * vd.template getProp<drho>(a);
				vd.template getProp<rho>(a) = (rhonew < rho_zero) ? rho_zero : rhonew;
				vd.template getProp<rho_prev>(a) = rhop;
			}

			++part;
			continue;
		}

		//-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
		double dx = vd.template getProp<velocity>(a)[0] * dt + (vd.template getProp<force>(a)[0] - vd.template getProp<Background_P>(a)[0]) * dt205;
		double dy = vd.template getProp<velocity>(a)[1] * dt + (vd.template getProp<force>(a)[1] - vd.template getProp<Background_P>(a)[1]) * dt205;
		double dz = vd.template getProp<velocity>(a)[2] * dt + (vd.template getProp<force>(a)[2] - vd.template getProp<Background_P>(a)[2]) * dt205;

		vd.getPos(a)[0] += dx;
		vd.getPos(a)[1] += dy;
		vd.getPos(a)[2] += dz;

		double velX = vd.template getProp<velocity>(a)[0];
		double velY = vd.template getProp<velocity>(a)[1];
		double velZ = vd.template getProp<velocity>(a)[2];
		double rhop = vd.template getProp<rho>(a);

		vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity_prev>(a)[0] + vd.template getProp<force>(a)[0] * dt2;
		vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity_prev>(a)[1] + vd.template getProp<force>(a)[1] * dt2;
		vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity_prev>(a)[2] + vd.template getProp<force>(a)[2] * dt2;
		vd.template getProp<rho>(a) = vd.template getProp<rho_prev>(a) + dt2 * vd.template getProp<drho>(a);

		vd.template getProp<velocity_prev>(a)[0] = velX;
		vd.template getProp<velocity_prev>(a)[1] = velY;
		vd.template getProp<velocity_prev>(a)[2] = velZ;
		vd.template getProp<rho_prev>(a) = rhop;

		++part;
	}

	// increment the iteration counter
	cnt++;
}

void euler_int(particles &vd, double dt)
{
	// particle iterator
	auto part = vd.getDomainIterator();

	double dt205 = dt * dt * 0.5;
	double dt2 = dt * 2.0;

	// For each particle ...
	while (part.isNext())
	{
		// ... a
		auto a = part.get();

		// if the particle is boundary
		if (vd.template getProp<type>(a) == BOUNDARY)
		{
			if (BC_TYPE == FREE_SLIP)
			{
				// Update rho
				double rhop = vd.template getProp<rho>(a);

				// Update only the density

				double rhonew = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
				vd.template getProp<rho>(a) = (rhonew < rho_zero) ? rho_zero : rhonew;

				vd.template getProp<rho_prev>(a) = rhop;
			}

			++part;
			continue;
		}

		//-Calculate displacement and update position / Calcula desplazamiento y actualiza posicion.
		double dx = vd.template getProp<velocity>(a)[0] * dt + (vd.template getProp<force>(a)[0] - vd.template getProp<Background_P>(a)[0]) * dt205;
		double dy = vd.template getProp<velocity>(a)[1] * dt + (vd.template getProp<force>(a)[1] - vd.template getProp<Background_P>(a)[1]) * dt205;
		double dz = vd.template getProp<velocity>(a)[2] * dt + (vd.template getProp<force>(a)[2] - vd.template getProp<Background_P>(a)[2]) * dt205;

		vd.getPos(a)[0] += dx;
		vd.getPos(a)[1] += dy;
		vd.getPos(a)[2] += dz;

		double velX = vd.template getProp<velocity>(a)[0];
		double velY = vd.template getProp<velocity>(a)[1];
		double velZ = vd.template getProp<velocity>(a)[2];
		double rhop = vd.template getProp<rho>(a);

		vd.template getProp<velocity>(a)[0] = vd.template getProp<velocity>(a)[0] + vd.template getProp<force>(a)[0] * dt;
		vd.template getProp<velocity>(a)[1] = vd.template getProp<velocity>(a)[1] + vd.template getProp<force>(a)[1] * dt;
		vd.template getProp<velocity>(a)[2] = vd.template getProp<velocity>(a)[2] + vd.template getProp<force>(a)[2] * dt;
		vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);

		vd.template getProp<velocity_prev>(a)[0] = velX;
		vd.template getProp<velocity_prev>(a)[1] = velY;
		vd.template getProp<velocity_prev>(a)[2] = velZ;
		vd.template getProp<rho_prev>(a) = rhop;

		++part;
	}

	// increment the iteration counter
	cnt++;
}

// void kick_drift_int(particles &vd, const double dt)
// {
// 	// particle iterator
// 	auto part = vd.getDomainIterator();

// 	const double dt_2 = dt * 0.5;

// 	// For each particle ...
// 	while (part.isNext())
// 	{
// 		// ... a
// 		auto a = part.get();

// 		// if the particle is boundary skip
// 		if (vd.template getProp<type>(a) == BOUNDARY)
// 		{
// 			// if (BC_TYPE == FREE_SLIP)
// 			// {
// 			// 	// Update rho
// 			// 	double rhop = vd.template getProp<rho>(a);

// 			// 	// Update only the density

// 			// 	double rhonew = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
// 			// 	vd.template getProp<rho>(a) = (rhonew < rho_zero) ? rho_zero : rhonew;

// 			// 	vd.template getProp<rho_prev>(a) = rhop;
// 			// }

// 			++part;
// 			continue;
// 		}

// 		//

// 		double vint_x = vd.template getProp<velocity>(a)[0] + vd.template getProp<force>(a)[0] * dt_2;
// 		double vint_y = vd.template getProp<velocity>(a)[1] + vd.template getProp<force>(a)[1] * dt_2;
// 		double vint_z = vd.template getProp<velocity>(a)[2] + vd.template getProp<force>(a)[2] * dt_2;

// 		double vtilde_x = vint_x + backgroundpressure * dt_2;
// 		double vtilde_y = vint_y + backgroundpressure * dt_2;
// 		double vtilde_z = vint_z + backgroundpressure * dt_2;

// 		vd.getPos(a)[0] += dt * vtilde_x;
// 		vd.getPos(a)[1] += dt * vtilde_y;
// 		vd.getPos(a)[2] += dt * vtilde_z;

// 		double velX = vd.template getProp<velocity>(a)[0];
// 		double velY = vd.template getProp<velocity>(a)[1];
// 		double velZ = vd.template getProp<velocity>(a)[2];
// 		double rhop = vd.template getProp<rho>(a);

// 		vd.template getProp<velocity>(a)[0] = vint_x + vd.template getProp<force>(a)[0] * dt;
// 		vd.template getProp<velocity>(a)[1] = vint_y + vd.template getProp<force>(a)[1] * dt;
// 		vd.template getProp<velocity>(a)[2] = vint_z + vd.template getProp<force>(a)[2] * dt;

// 		vd.template getProp<rho>(a) = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);

// 		vd.template getProp<velocity_prev>(a)[0] = velX;
// 		vd.template getProp<velocity_prev>(a)[1] = velY;
// 		vd.template getProp<velocity_prev>(a)[2] = velZ;
// 		vd.template getProp<rho_prev>(a) = rhop;

// 		++part;
// 	}

// 	// increment the iteration counter
// 	cnt++;
// }

template <typename Vector, typename CellList>
inline void sensor_pressure(Vector &vd,
							CellList &NN,
							openfpm::vector<openfpm::vector<double>> &press_t,
							openfpm::vector<Point<3, double>> &probes)
{
	Vcluster<> &v_cl = create_vcluster();

	press_t.add();

	for (size_t i = 0; i < probes.size(); i++)
	{
		float press_tmp = 0.0f;
		float tot_ker = 0.0;

		// if the probe is inside the processor domain
		if (vd.getDecomposition().isLocal(probes.get(i)) == true)
		{
			// Get the position of the probe i
			Point<3, double> xp = probes.get(i);

			// get the iterator over the neighbohood particles of the probes position
			auto itg = NN.getNNIterator(NN.getCell(probes.get(i)));
			while (itg.isNext())
			{
				auto q = itg.get();

				// Only the fluid particles are importants
				if (vd.template getProp<type>(q) != FLUID)
				{
					++itg;
					continue;
				}

				// Get the position of the neighborhood particle q
				Point<3, double> xq = vd.getPos(q);

				// Calculate the contribution of the particle to the pressure
				// of the probe
				double r = sqrt(norm2(xp - xq));

				double ker = Wab(r) * (MassFluid / rho_zero);

				// Also keep track of the calculation of the summed
				// kernel
				tot_ker += ker;

				// Add the total pressure contribution
				press_tmp += vd.template getProp<Pressure>(q) * ker;

				// next neighborhood particle
				++itg;
			}

			// We calculate the pressure normalizing the
			// sum over all kernels
			if (tot_ker == 0.0)
				press_tmp = 0.0;
			else
				press_tmp = 1.0 / tot_ker * press_tmp;
		}

		// This is not necessary in principle, but if you
		// want to make all processor aware of the history of the calculated
		// pressure we have to execute this
		v_cl.sum(press_tmp);
		v_cl.execute();

		// We add the calculated pressure into the history
		press_t.last().add(press_tmp);
	}
}

int main(int argc, char *argv[])
{

	// initialize the library
	openfpm_init(&argc, &argv);

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

	if (SCENARIO == POISEUILLE || SCENARIO == COUETTE)
	{
		if (SCENARIO == POISEUILLE)
		{
			gravity_vector = {0.0, 0.0, 0.1};
			gravity = sqrt(norm2(gravity_vector));
			vw_top = {0.0, 0.0, 0.0};
			vw_bottom = {0.0, 0.0, 0.0};
		}
		else if (SCENARIO == COUETTE)
		{
			gravity_vector = {0.0, 0.0, 0.0};
			gravity = sqrt(norm2(gravity_vector));
			vw_top = {0.0, 0.0, 1.25};
			vw_bottom = {0.0, 0.0, 0.0};
		}

		// Number of fluid and boundary particles in the x, y and z direction
		size_t Np_fluid[3] = {40, 1, 60};
		size_t Nbound_x = (BC_TYPE == NEW_NO_SLIP) ? 1 : 3;
		size_t Np_boundary[3] = {Nbound_x, 0, 0};

		bc[0] = NON_PERIODIC;
		bc[1] = NON_PERIODIC;
		bc[2] = PERIODIC;

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
			if (bc[dim] == 0)
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
				offset_periodic_fluid[dim] = 0.55 * dp;
				offset_periodic_recipient[dim] = 0.75 * dp;
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

		Box<3, double> box_recipient({-offset_recipient[0],
									  -offset_recipient[1],
									  -offset_recipient[2]},
									 {length[0] + offset_recipient[0] + offset_periodic_recipient[0],
									  length[1] + offset_recipient[1] + offset_periodic_recipient[1],
									  length[2] + offset_recipient[2] + offset_periodic_recipient[2]});

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
	Ghost<3, double> g(2.1 * Hconst * dp);

	particles vd(0, domain, bc, g, DEC_GRAN(512));

	// return an iterator to the fluid particles to add to vd
	auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

	// Fill needed constants
	double L = length[0];
	double computed_nu; // actual viscosity computed from alpha
	double umax;

	if (SCENARIO == POISEUILLE)
	{
		umax = gravity_vector.get(2) * L * L / (8.0 * nu);
	}
	else if (SCENARIO == COUETTE)
	{
		umax = vw_top.get(2);
	}

	cbar = coeff_sound * umax;
	B = cbar * cbar * rho_zero / gamma_;
	Pbackground = B;
	double Re = L * umax / nu;

	if (VISC_TYPE == ARTIFICIAL_VISCOSITY)
	{
		visco = nu * 2.0 * (dimensions + 2.0) / (H * cbar);
		computed_nu = (1.0 / (2.0 * (dimensions + 2.0))) * visco * H * cbar;
	}
	else if (VISC_TYPE == PHYSICAL_VISCOSITY)
	{
		computed_nu = nu;
	}

	std::string constants_filename = filename + "_Constants_" + ".txt";
	std::ofstream file(constants_filename);
	file << "gravity_x:  " << gravity_vector.get(0) << std::endl;
	file << "gravity_y:  " << gravity_vector.get(1) << std::endl;
	file << "gravity_z:  " << gravity_vector.get(2) << std::endl;
	file << "Re: " << Re << std::endl;
	file << "L: " << L << std::endl;
	file << "nu: " << nu << std::endl;
	file << "computed_nu: " << computed_nu << std::endl;
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
	file << "BC_TYPE: " << BC_TYPE << std::endl;
	file << "VISC_TYPE: " << VISC_TYPE << std::endl;
	file << "PRESSURE_TYPE: " << PRESSURE_TYPE << std::endl;
	file << "SCENARIO: " << SCENARIO << std::endl;

	// for each particle inside the fluid box ...
	const double profile_parameter = gravity_vector.get(2) / (2.0 * nu);
	while (fluid_it.isNext())
	{
		// ... add a particle ...
		vd.add();

		// ... and set it position ...
		vd.getLastPos()[0] = fluid_it.get().get(0);
		vd.getLastPos()[1] = fluid_it.get().get(1);
		vd.getLastPos()[2] = fluid_it.get().get(2);

		// and its type.
		vd.template getLastProp<type>() = FLUID;

		// We also initialize the density of the particle and the hydro-static pressure given by
		//
		// rho_zero*g*h = P
		//
		// rho_p = (P/B + 1)^(1/Gamma) * rho_zero
		//

		vd.template getLastProp<Pressure>() = 0.0;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<rho_prev>() = rho_zero;
		vd.template getLastProp<velocity>()[0] = initial_perturbation * (-1.0 + (double)2.0 * rand() / RAND_MAX) * (L - std::abs(fluid_it.get().get(0) - L / 2.0));
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = initial_perturbation * (-1.0 + (double)2.0 * rand() / RAND_MAX) * (L - std::abs(fluid_it.get().get(0) - L / 2.0));

		// profile_parameter * fluid_it.get().get(0) * (L - fluid_it.get().get(0));

		vd.template getLastProp<velocity_prev>()[0] = 0.0;
		vd.template getLastProp<velocity_prev>()[1] = 0.0;
		vd.template getLastProp<velocity_prev>()[2] = 0.0;

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
		vd.template getLastProp<Pressure>() = B;

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

	vd.map();

	// Now that we fill the vector with particles
	ModelCustom md;

	vd.addComputationCosts(md);
	vd.getDecomposition().decompose();
	vd.map();

	vd.ghost_get<type, rho, Pressure, velocity>();

	auto NN = vd.getCellList(2.0 * H);

	openfpm::vector<std::string> names({"type", "rho", "rho_prev", "pressure", "drho", "force", "velocity", "velocity_prev", "Background_P"});
	vd.setPropNames(names);

	// Evolve
	size_t write = 0;
	size_t it = 0;
	size_t it_reb = 0;
	double t = 0.0;
	while (t <= t_end)
	{
		Vcluster<> &v_cl = create_vcluster();
		timer it_time;

		////// Do rebalancing every 200 timesteps
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

		vd.map();

		// Calculate pressure from the density
		EqState(vd);

		double max_visc = 0.0;

		vd.ghost_get<type, rho, Pressure, velocity>();

		// In no slip bc, we need to compute a velocity for boundary particles
		if (BC_TYPE == NO_SLIP)
		{
			// 	Boundary condition
			calc_boundary(vd, NN);
		}

		// Calc forces
		calc_forces(vd, NN, max_visc);

		// Get the maximum viscosity term across processors
		v_cl.max(max_visc);
		v_cl.execute();

		// Calculate delta t integration
		double dt = calc_deltaT(vd, max_visc);

		// VerletStep or euler step
		it++;
		if (it < 40)
			verlet_int(vd, dt);
		else
		{
			euler_int(vd, dt);
			it = 0;
		}

		t += dt;

		if (write < t * write_const)
		{
			// sensor_pressure calculation require ghost and update cell-list
			// vd.map();
			// vd.ghost_get<type, rho, Pressure, velocity>();
			// vd.updateCellList(NN);

			// // calculate the pressure at the sensor points
			// // sensor_pressure(vd, NN, press_t, probes);

			vd.deleteGhost();
			vd.write_frame(filename, write, WRITER);
			vd.ghost_get<type, rho, Pressure, velocity>();

			write++;

			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "TIME: " << t << "  write " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "   " << cnt << "   Max visc: " << max_visc << std::endl;
			}
		}
		else
		{
			if (v_cl.getProcessUnitID() == 0)
			{
				std::cout << "TIME: " << t << "  " << it_time.getwct() << "   " << v_cl.getProcessUnitID() << "   " << cnt << "    Max visc: " << max_visc << std::endl;
			}
		}
	}

	openfpm_finalize();
}
