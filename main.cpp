#include "Vector/vector_dist.hpp"
#include "Draw/DrawParticles.hpp"
#include <math.h>
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

//////// CODE PARAMETERS /////////////////////////////////////////////////////////////
// Dimensionality of the problem changes the normalizaiton of the kernel
const int dimensions = 2;

// BC and viscosity we are using in the simulation
const int BC_TYPE = NO_SLIP;
const int VISC_TYPE = ARTIFICIAL_VISCOSITY;
const int PRESSURE_TYPE = NO_PRESSURE;
const int SCENARIO = COUETTE;
const int WRITER = VTK_WRITER; // VTK_WRITER or CSV_WRITER
// Output file name
std::string filename = "COUETTE_ArtVisc_NoPressure";

// Controls otput file frequency, low means less frequent
const int write_const = 1;
////////////////////////////////////////////////////////////////////////////////////

//////// SPATIAL CONSTANTS /////////////////////////////////////////////////////////
// initial spacing between particles dp in the formulas
const double dp = 0.025;

// Factor relating H (smoothing length) and dp ( particle spacing)
const double Hconst = std::sqrt(3.0);

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

////////////////////////////////////////////////////////////////////////////////////

//////// TIME CONSTANTS //////////////////////////////////////////////////////////////
// End simulation time
const double t_end = 1;
// Constant used to define time integration
const double CFLnumber = 0.2;
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

typedef vector_dist<3, double, aggregate<size_t, double, double, double, double, double[3], double[3], double[3]>> particles;
//                                       |      |        |          |            |            |         |            |
//                                       |      |        |          |            |            |         |            |
//                                     type   density   density    Pressure    delta       force     velocity    velocity
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

		vd.template getProp<Pressure>(a) = B * (rho_frac * rho_frac * rho_frac * rho_frac * rho_frac * rho_frac * rho_frac - 1.0);

		++it;
	}
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

inline void DWab(Point<3, double> &dx, Point<3, double> &DW, double r, bool print)
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

inline double Tensile(double r, double rhoa, double rhob, double prs1, double prs2)
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

inline double Pi_artificial(const Point<3, double> &dr, double rr2, Point<3, double> &dv, double rhoa, double rhob, double massb, double &visc)
{
	const double dot = dr.get(0) * dv.get(0) + dr.get(1) * dv.get(1) + dr.get(2) * dv.get(2);
	const double dot_rr2 = dot / (rr2 + Eta2);
	visc = std::max(dot_rr2, visc);

	const float amubar = H * dot_rr2;
	const float robar = (rhoa + rhob) * 0.5;
	const float pi_visc = (-visco * cbar * amubar / robar);

	return pi_visc;
}
inline Point<3, double> Pi_physical(const Point<3, double> &dr, double r, Point<3, double> &dv, Point<3, double> &dW, double rhoa, double rhob, double massa, double massb, double &visc)
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

// inline double interact(Point<3, double> &va, Point<3, double> &dr, bool bottomwall)
// {
// 	Point<3, double> normal_bottom_wall;
// 	normal_bottom_wall.get(0) = 1.0;
// 	normal_bottom_wall.get(1) = 0.0;
// 	normal_bottom_wall.get(2) = 0.0;

// 	Point<3, double> normal_top_wall;
// 	normal_top_wall.get(0) = -1.0;
// 	normal_top_wall.get(1) = 0.0;
// 	normal_top_wall.get(2) = 0.0;

// 	Point<3, double> normal;
// 	if (bottomwall)
// 	{
// 		normal = normal_bottom_wall;
// 	}
// 	else
// 	{
// 		normal = normal_top_wall;
// 	}

// 	Point<3, double> offset_1 = -normal * dp / 2.0; // offset from wall to first particle
// 	Point<3, double> offset_2 = -normal * dp;		// offset from first to second particle, and from second to third

// 	Point<3, double> r1 = dr + offset_1;
// 	double r1_norm = sqrt(norm2(r1));
// }

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
		auto b = part.get();

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
				auto f = Np.get();

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
				vd.template getProp<rho>(b) = pow((vd.template getProp<Pressure>(b) / B + 1.0), 1.0 / gamma_);
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

template <typename CellList>
inline void calc_forces(particles &vd, CellList &NN, double &max_visc)
{
	auto part = vd.getDomainIterator();

	// Update the cell-list
	vd.updateCellList(NN);

	// For each particle ...
	while (part.isNext())
	{
		// get particle a
		auto a = part.get();

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
				auto b = Np.get();

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
					DWab(dr, DW, r, false);

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
					{
					}

					if (BC_TYPE == NO_SLIP)
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

			// For each neighborhood particle
			while (Np.isNext() == true)
			{
				// get particle b
				auto b = Np.get();

				// if (p == q) skip this particle
				if (a.getKey() == b)
				{
					++Np;
					continue;
				};
				// Get the position xp of the particle
				Point<3, double> xb = vd.getPos(b);

				// Get the distance between a and b
				Point<3, double> dr = xa - xb;

				// take the norm (squared) of this vector
				double r2 = norm2(dr);

				// if they interact
				if (r2 < 4.0 * H * H)
				{

					// get mass of the particle b
					double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;

					// get the velocity of the particle b
					Point<3, double> vb = vd.getProp<velocity>(b);

					// If the particle is boundary, and we are using the No Slip conditon this contains
					// the velocity to be used in the momentum equation to impose BC otherwise it is useless
					Point<3, double> vb_prev = vd.getProp<velocity_prev>(b);

					// Get the pressure and density of the particle b
					double Pb = vd.getProp<Pressure>(b);
					double rhob = vd.getProp<rho>(b);

					// calculate distance
					double r = sqrt(r2);

					// In no slip case
					// interactions with fluid particles use same vrel in momentum and continuity equation
					// interactions fluid boundary use vrel_aux in momentum equation and vrel in continuity equation in the no slip case
					Point<3, double> v_rel = va - vb;
					Point<3, double> v_rel_aux;

					if (vd.getProp<type>(b) == FLUID) // fluid particle
					{								  // fluid particles use same vrel in momentum and continuity equation
						v_rel_aux = v_rel;
					}
					else // boundary particle
					{
						if (BC_TYPE == NO_SLIP) // we need to compute v_rel using the assigned velocity to the boundary particle
						{
							v_rel_aux = va - vb_prev;
						}
						else if (BC_TYPE == FREE_SLIP) // we need to set the boundary particle velocity to 0 ( it should already be)
						{
							v_rel_aux = va - 0.0;
						}
						else if (BC_TYPE == NEW_NO_SLIP)
						{

							// double v_boundary = Interact(va,dr);
							/// v_rel_aux = va - v_boundary;
						}
					}

					// Evaluate gradient of the kernel
					Point<3, double> DW;
					DWab(dr, DW, r, false);

					Point<3, double> PhysicalViscosityTerm;
					double ArtificalViscosityTerm;
					double PressureTerm;
					if (PRESSURE_TYPE == OLD_PRESSURE)
					{
						PressureTerm = PressureForce(rhoa, rhob, Pa, Pb, massa, massb);
					}
					else if (PRESSURE_TYPE == NEW_PRESSURE)
					{
						PressureTerm = PressureForceNew(rhoa, rhob, Pa, Pb, massa, massb);
					}
					else if (PRESSURE_TYPE == NO_PRESSURE)
					{
						PressureTerm = 0.0;
					}

					double TensileTerm = -massb * Tensile(r, rhoa, rhob, Pa, Pb);

					double factor = TensileTerm + PressureTerm;

					if (VISC_TYPE == PHYSICAL_VISCOSITY)
					{
						PhysicalViscosityTerm = Pi_physical(dr, r, v_rel_aux, DW, rhoa, rhob, massa, massb, max_visc) / massa;
						vd.getProp<force>(a)[0] += PhysicalViscosityTerm.get(0);
						vd.getProp<force>(a)[1] += PhysicalViscosityTerm.get(1);
						vd.getProp<force>(a)[2] += PhysicalViscosityTerm.get(2);
					}
					else if (VISC_TYPE == ARTIFICIAL_VISCOSITY)
					{
						ArtificalViscosityTerm = -massb * Pi_artificial(dr, r2, v_rel_aux, rhoa, rhob, massa, max_visc);
						factor += ArtificalViscosityTerm;
					}
					vd.getProp<force>(a)[0] += factor * DW.get(0);
					vd.getProp<force>(a)[1] += factor * DW.get(1);
					vd.getProp<force>(a)[2] += factor * DW.get(2);

					vd.getProp<drho>(a) += massb * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
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

	//-dt1 depends on force per unit mass.
	const double dt_f = (Maxacc) ? sqrt(H / Maxacc) : std::numeric_limits<int>::max();

	//-dt2 combines the Courant and the viscous time-step controls.
	// const double dt_cv = H / (std::max(cbar, Maxvel * 10.) + H * ViscDtMax);

	const double dt_cv = H / (std::max(cbar, Maxvel * 10.) + H * ViscDtMax);
	//-dt new value of time step.
	double dt = double(CFLnumber) * std::min(dt_f, dt_cv);
	if (dt < double(DtMin))
		dt = double(DtMin);

	return dt;
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
		double dx = vd.template getProp<velocity>(a)[0] * dt + vd.template getProp<force>(a)[0] * dt205;
		double dy = vd.template getProp<velocity>(a)[1] * dt + vd.template getProp<force>(a)[1] * dt205;
		double dz = vd.template getProp<velocity>(a)[2] * dt + vd.template getProp<force>(a)[2] * dt205;

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
		double dx = vd.template getProp<velocity>(a)[0] * dt + vd.template getProp<force>(a)[0] * dt205;
		double dy = vd.template getProp<velocity>(a)[1] * dt + vd.template getProp<force>(a)[1] * dt205;
		double dz = vd.template getProp<velocity>(a)[2] * dt + vd.template getProp<force>(a)[2] * dt205;

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

	bool bc_new = false;

	// Number of fluid and boundary particles in the x, y and z direction
	size_t Np_fluid[3] = {40, 1, 60};
	size_t Np_boundary[3] = {3, 0, 0};

	// Here we define the boundary conditions of our problem
	size_t bc[3] = {NON_PERIODIC, NON_PERIODIC, PERIODIC};

	// Size of the virtual grid that defines where to place the particles
	size_t sz[3];

	// In the case of the new bc we need particles at the wall, for this we need sz_aux
	// We want to put one virtual grid point between each pair of the old ones,
	// so that the new spacing is dp/2, and we can put a fluid particle exactly at the wall
	size_t sz_aux[3];

	// Physical size of the fluid domain, it goes from (0,0,0) to (length[0],length[1],length[2])
	// First particle will always be placed at (dp/2,dp/2,dp/2) and the last particle will be placed at (length[0]-dp/2,length[1]-dp/2,length[2]-dp/2)
	double length[3];
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

			if (bc_new == 1) // Np_boundary should only be 0 or 1 if we are using the new bc
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

	// last position of fluid particles in z coordinate is at
	// length[2] + 0.5*dp
	double z_end = length[2] + 0.5 * dp;

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

	// if (VISC_TYPE == ARTIFICIAL_VISCOSITY)
	// {
	// 	if (SCENARIO == POISEUILLE)
	// 	{
	// 		visco = (64.0 / 10.0) * (nu * nu / (H * gravity * L));
	// 		computed_nu = sqrt(10.0 * visco * H * gravity * L * L / 64.0);
	// 		umax = gravity_vector.get(2) * L * L / (8.0 * nu);
	// 	}
	// 	else if (SCENARIO == COUETTE)
	// 	{
	// 		umax = vw_top.get(2);
	// 		visco = (2.0 * (dimensions + 2) * nu) / (coeff_sound * H * umax);
	// 		computed_nu = (1.0 / (2 * (dimensions + 2))) * visco * H * coeff_sound * umax;
	// 	}
	// }
	// else if (VISC_TYPE == PHYSICAL_VISCOSITY)
	// {

	// }

	std::string constants_filename = "Constants_" + filename + ".txt";
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
	file << "H/dp: " << H / dp << std::endl;
	file << "computed_nu: " << computed_nu << std::endl;

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
		vd.template getLastProp<velocity>()[0] = 0.1 * (-1.0 + (double)2.0 * rand() / RAND_MAX) * (L - std::abs(fluid_it.get().get(0) - L / 2.0));
		vd.template getLastProp<velocity>()[1] = 0.0;
		// impose theoretical pouiselle flow profile for faster convergence
		vd.template getLastProp<velocity>()[2] = 0.1 * (-1.0 + (double)2.0 * rand() / RAND_MAX) * (L - std::abs(fluid_it.get().get(0) - L / 2.0));

		// profile_parameter * fluid_it.get().get(0) * (L - fluid_it.get().get(0));

		vd.template getLastProp<velocity_prev>()[0] = 0.0;
		vd.template getLastProp<velocity_prev>()[1] = 0.0;
		vd.template getLastProp<velocity_prev>()[2] = 0.0;

		// next fluid particle
		++fluid_it;
	}

	// Recipient

	openfpm::vector<Box<3, double>> holes;

	if (bc_new == true)
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
		if ((!recipient.isInside((position))) || (fluid_box.isInside(position) && !(position.get(2) < z_end + epsilon && position.get(2) > z_end - epsilon)))
		{
			++bound_box;
			continue;
		}

		vd.add();

		vd.getLastPos()[0] = bound_box.get().get(0);
		vd.getLastPos()[1] = bound_box.get().get(1);
		vd.getLastPos()[2] = bound_box.get().get(2);

		vd.template getLastProp<type>() = BOUNDARY;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<rho_prev>() = rho_zero;

		if (bound_box.get().get(0) < 0) // bottom wall
		{
			vd.template getLastProp<velocity>()[0] = vw_bottom.get(0);
			vd.template getLastProp<velocity>()[1] = vw_bottom.get(1);
			vd.template getLastProp<velocity>()[2] = vw_bottom.get(2);
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

	openfpm::vector<std::string> names({"type", "rho", "rho_prev", "pressure", "drho", "force", "velocity", "velocity_prev"});
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
