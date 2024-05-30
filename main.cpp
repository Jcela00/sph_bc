#include "Vector/vector_dist.hpp"
#include <math.h>
#include "Draw/DrawParticles.hpp"

// A constant to indicate type of particles
#define BOUNDARY 0
#define FLUID 1

// A constant to indicate the type of boundary condition
#define NO_SLIP 0
#define FREE_SLIP 1
#define NEW_NO_SLIP 2

// output file name
std::string filename = "test_prof";

// bool variable to indicate the dimensionality of the problem
// changes the normalizaiton of the kernel
// true for 3D, false for 2D
bool three_dimensional = false;
const int dimensions = three_dimensional ? 3 : 2;

// BC we are using in the simulation
const int BC_TYPE = NO_SLIP;

// initial spacing between particles dp in the formulas
// const double dp = 0.0085;
const double dp = 0.025;

// c_s in the formulas (constant used to calculate the sound speed)
// This means that the sound speed is coeff_sound times the maximum velocity
const double coeff_sound = 10.0;

// gamma in the formulas
const double gamma_ = 7.0;

// sqrt(3.0*dp*dp) support of the kernel
const double H = std::sqrt(3.0) * dp;

// Eta in the formulas
const double Eta2 = 0.01 * H * H;

// constant used to write files
const int write_const = 100;

// cbar in the formula (calculated later)
double cbar = 0.0;

// End simulation time
#ifdef TEST_RUN
const double t_end = 0.001;
#else
const double t_end = 100;
#endif

// Gravity acceleration
const double gravity = 0.1;

// double visco = 0.1;
double visco = 0.14;
// new viscosity
double eta = 0.1;

// Reference density 1000Kg/m^3
const double rho_zero = 1.0;

// Minimum Rho allowed
const double RhoMin = 0.7 * rho_zero;

// Maximum Rho allowed
const double RhoMax = 1.3 * rho_zero;

std::ofstream FileForce("force_check.txt");

// Mass of the fluid particles
// M=rho*V, V=dp^3 or V = dp^2
const double MassFluid = rho_zero * (three_dimensional ? dp * dp * dp : dp * dp);
// Mass of the boundary particles
const double MassBound = rho_zero * (three_dimensional ? dp * dp * dp : dp * dp);

// Filled later require h_swl, it is b in the formulas
double B = 0.0;

// Constant used to define time integration
const double CFLnumber = 0.2;

// Minimum T
const double DtMin = 0.00001;

// Filled in initialization
double max_fluid_height = 0.0;

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

// const double a2 = 1.0 / M_PI / H / H / H;
// in 3d the normalization is 1/(pi*h^3) in 2d is 10/(7pi*h^2)
const double a2 = three_dimensional ? 1.0 / M_PI / H / H / H : 10.0 / 7.0 / M_PI / H / H;

inline double Wab(double r)
{
	r /= H;
	if (r < 1.0)
		return (1.0 - 3.0 / 2.0 * r * r + 3.0 / 4.0 * r * r * r) * a2;
	else if (r < 2.0)
	{
		// I think there is a mistake in the original code, r^2 should be r
		// return (1.0 / 4.0 * (2.0 - r * r) * (2.0 - r * r) * (2.0 - r * r)) * a2;
		return (1.0 / 4.0 * (2.0 - r) * (2.0 - r) * (2.0 - r)) * a2;
	}
	else
		return 0.0;
}

// const double c1 = -3.0 / M_PI / H / H / H / H;
// const double d1 = 9.0 / 4.0 / M_PI / H / H / H / H;
// const double c2 = -3.0 / 4.0 / M_PI / H / H / H / H;

const double c1 = (-3.0 / H) * a2;
const double d1 = (9.0 / 4.0 / H) * a2;
const double c2 = (-3.0 / 4.0 / H) * a2;
const double a2_4 = 0.25 * a2;
// Filled later
double W_dap = 0.0;

inline void DWab(Point<3, double> &dx, Point<3, double> &DW, double r, bool print)
{
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

// Tensile correction
inline double Tensile(double r, double rhoa, double rhob, double prs1, double prs2)
{
	const double qq = r / H;
	//-Cubic Spline kernel
	double wab;
	if (r > H)
	{
		double wqq1 = 2.0 - qq;
		double wqq2 = wqq1 * wqq1;

		wab = a2_4 * (wqq2 * wqq1);
	}
	else
	{
		double wqq2 = qq * qq;
		double wqq3 = wqq2 * qq;

		wab = a2 * (1.0 - 1.5 * wqq2 + 0.75 * wqq3);
	}

	//-Tensile correction.
	double fab = wab * W_dap;
	fab *= fab;
	fab *= fab; // fab=fab^4
	const double tensilp1 = (prs1 / (rhoa * rhoa)) * (prs1 > 0 ? 0.01 : -0.2);
	const double tensilp2 = (prs2 / (rhob * rhob)) * (prs2 > 0 ? 0.01 : -0.2);

	return (fab * (tensilp1 + tensilp2));
}

inline double Pi(const Point<3, double> &dr, double rr2, Point<3, double> &dv, double rhoa, double rhob, double massb, double &visc)
{
	const double dot = dr.get(0) * dv.get(0) + dr.get(1) * dv.get(1) + dr.get(2) * dv.get(2);
	const double dot_rr2 = dot / (rr2 + Eta2);
	visc = std::max(dot_rr2, visc);

	// if (dot < 0)
	// {
	const float amubar = H * dot_rr2;
	const float robar = (rhoa + rhob) * 0.5;
	const float pi_visc = (-visco * cbar * amubar / robar);

	return pi_visc;
	// }
	// else
	// 	return 0.0;
}
inline Point<3, double> Pi_new(const Point<3, double> &dr, double r, Point<3, double> &dv, Point<3, double> &dW, double rhoa, double rhob, double massa, double massb, double &visc)
{
	const Point<3, double> e_ab = dr / r; // unit vector from a to b
	const double dot = e_ab.get(0) * dW.get(0) + e_ab.get(1) * dW.get(1) + e_ab.get(2) * dW.get(2);
	double Va2 = (massa / rhoa) * (massa / rhoa);
	double Vb2 = (massb / rhob) * (massb / rhob);
	Point<3, double> pi_visc = ((Va2 + Vb2) * dot * dv) / (r + Eta2);
	const double normvisc = sqrt(norm2(pi_visc));
	visc = std::max(normvisc, visc);
	return pi_visc;
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
		// Key of particle a
		auto a = part.get();

		// if particle boundary
		if (vd.getProp<type>(a) != FLUID)
		{

			// Get the position xa of the particle a
			Point<3, double> xa = vd.getPos(a);

			// Get the Velocity of the particle a
			Point<3, double> va = vd.getProp<velocity>(a);

			// initialize sums
			Point<3, double> sum_vW;
			sum_vW.get(0) = 0.0;
			sum_vW.get(1) = 0.0;
			sum_vW.get(2) = 0.0;
			double sum_W = 0.0;

			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));
			// iterate the neighborhood particles
			while (Np.isNext() == true)
			{
				// Key of particle b
				auto b = Np.get();

				// if (a == b) skip this particle
				if (a.getKey() == b)
				{
					++Np;
					continue;
				};
				// Skip other boundary particles
				if (vd.getProp<type>(b) != FLUID)
				{
					++Np;
					continue;
				}

				// Get the position xb of the particle b
				Point<3, double> xb = vd.getPos(b);

				// Get the velocity of the particle b
				Point<3, double> vb = vd.getProp<velocity>(b);

				// Get the distance between a and b
				Point<3, double> dr = xa - xb;

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
					sum_vW.get(0) += w * vb.get(0);
					sum_vW.get(1) += w * vb.get(1);
					sum_vW.get(2) += w * vb.get(2);
					sum_W += w;
				}

				++Np;
			}
			if (sum_W != 0.0)
			{
				// Update the velocity of the particle a
				vd.template getProp<velocity_prev>(a)[0] = -sum_vW.get(0) / sum_W;
				vd.template getProp<velocity_prev>(a)[1] = -sum_vW.get(1) / sum_W;
				vd.template getProp<velocity_prev>(a)[2] = -sum_vW.get(2) / sum_W;
			}
			else
			{
				vd.template getProp<velocity_prev>(a)[0] = 0.0;
				vd.template getProp<velocity_prev>(a)[1] = 0.0;
				vd.template getProp<velocity_prev>(a)[2] = 0.0;
			}
		}

		++part;
	}
}

template <typename CellList>
inline void calc_forces(particles &vd, CellList &NN, double &max_visc, bool flag)
{
	auto part = vd.getDomainIterator();

	// Update the cell-list
	vd.updateCellList(NN);

	// For each particle ...
	while (part.isNext())
	{
		// ... a
		auto a = part.get();

		// Get the position xp of the particle
		Point<3, double> xa = vd.getPos(a);

		// Take the mass of the particle dependently if it is FLUID or BOUNDARY
		double massa = (vd.getProp<type>(a) == FLUID) ? MassFluid : MassBound;

		// Get the density of the of the particle a
		double rhoa = vd.getProp<rho>(a);

		// Get the pressure of the particle a
		double Pa = vd.getProp<Pressure>(a);

		// Get the Velocity of the particle a
		Point<3, double> va = vd.getProp<velocity>(a);

		// Reset the force counter (- gravity on zeta direction)
		vd.template getProp<force>(a)[0] = 0.0;
		vd.template getProp<force>(a)[1] = 0.0;
		vd.template getProp<force>(a)[2] = gravity;
		vd.template getProp<drho>(a) = 0.0;

		// We threat FLUID particle differently from BOUNDARY PARTICLES ...
		if (vd.getProp<type>(a) != FLUID)
		{
			// If it is a boundary particle calculate the delta rho based on equation 2
			// This require to run across the neighborhoods particles of a
			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

			// For each neighborhood particle
			while (Np.isNext() == true)
			{
				// ... q
				auto b = Np.get();

				// Get the position xb of the particle
				Point<3, double> xb = vd.getPos(b);

				// if (a == b) skip this particle
				if (a.getKey() == b)
				{
					++Np;
					continue;
				};

				// get the mass of the particle
				double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;

				// Get the velocity of the particle b
				Point<3, double> vb = vd.getProp<velocity>(b);

				// Get the pressure and density of particle b
				double Pb = vd.getProp<Pressure>(b);
				double rhob = vd.getProp<rho>(b);

				// Get the distance between a and b
				Point<3, double> dr = xa - xb;
				// take the norm squared of this vector
				double r2 = norm2(dr);

				// If the particles interact ...
				if (r2 < 4.0 * H * H)
				{
					// ... calculate distance
					double r = sqrt(r2);

					Point<3, double> dv = va - vb;

					Point<3, double> DW;
					DWab(dr, DW, r, false);

					const double dot = dr.get(0) * dv.get(0) + dr.get(1) * dv.get(1) + dr.get(2) * dv.get(2);
					const double dot_rr2 = dot / (r2 + Eta2);
					max_visc = std::max(dot_rr2, max_visc);

					// call the viscosity funciton just to update max_visc
					// Point<3, double> dummy = Pi_new(dr, r, dv, DW, rhoa, rhob, massa, massb, max_visc);

					vd.getProp<drho>(a) += massb * (dv.get(0) * DW.get(0) + dv.get(1) * DW.get(1) + dv.get(2) * DW.get(2));
				}

				++Np;
			}
		}
		else
		{
			// If it is a fluid particle calculate based on equation 1 and 2

			// Get an iterator over the neighborhood particles of p
			auto Np = NN.getNNIterator(NN.getCell(vd.getPos(a)));

			// For each neighborhood particle
			while (Np.isNext() == true)
			{
				// ... q
				auto b = Np.get();

				// Get the position xp of the particle
				Point<3, double> xb = vd.getPos(b);

				// if (p == q) skip this particle
				if (a.getKey() == b)
				{
					++Np;
					continue;
				};

				// get mass of the particle b
				double massb = (vd.getProp<type>(b) == FLUID) ? MassFluid : MassBound;

				// get the velocity of the particle b
				Point<3, double> vb = vd.getProp<velocity>(b);

				// If the particle is boundary, this contains the velocity to be used in the momentum equation to impose BC
				// otherwise it is useless
				Point<3, double> vb_prev = vd.getProp<velocity_prev>(b);

				// Get the pressure and density of the particle b
				double Pb = vd.getProp<Pressure>(b);
				double rhob = vd.getProp<rho>(b);

				// Get the distance between a and b
				Point<3, double> dr = xa - xb;
				// take the norm of this vector
				double r2 = norm2(dr);

				// if they interact
				if (r2 < 4.0 * H * H)
				{
					double r = sqrt(r2);

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

					Point<3, double> DW;
					DWab(dr, DW, r, false);
					double fakevisc = 0.0;
					// Point<3, double> visc_term = Pi_new(dr, r, v_rel_aux, DW, rhoa, rhob, massa, massb, max_visc) * rho_zero / massb / massb;
					double PressureTerm = -massb * (vd.getProp<Pressure>(a) + vd.getProp<Pressure>(b)) / (rhoa * rhob);
					double TensileTerm = -massb * Tensile(r, rhoa, rhob, Pa, Pb);
					double ViscosityTerm = -massb * Pi(dr, r2, v_rel_aux, rhoa, rhob, massa, max_visc);
					// double factor = -massb * ((vd.getProp<Pressure>(a) + vd.getProp<Pressure>(b)) / (rhoa * rhob) + Tensile(r, rhoa, rhob, Pa, Pb) + Pi(dr, r2, v_rel_aux, rhoa, rhob, massb, max_visc));

					double factor = PressureTerm + TensileTerm + ViscosityTerm;
					vd.getProp<force>(a)[0] += factor * DW.get(0); // + visc_term.get(0);
					vd.getProp<force>(a)[1] += factor * DW.get(1); // + visc_term.get(1);
					vd.getProp<force>(a)[2] += factor * DW.get(2); // + visc_term.get(2);

					// if (flag == 1)
					// {
					// 	double Viscosity_x = ViscosityTerm * DW.get(0);
					// 	double Viscosity_y = ViscosityTerm * DW.get(1);
					// 	double Viscosity_z = ViscosityTerm * DW.get(2);
					// 	FileForce << "V_X: " << std::left << std::setw(18) << Viscosity_x << " V_Y: " << std::left << std::setw(18) << Viscosity_y << " V_Z: " << std::left << std::setw(18) << Viscosity_z << std::endl;
					// 	FileForce << "V_X: " << std::left << std::setw(18) << visc_term.get(0) << " V_Y: " << std::left << std::setw(18) << visc_term.get(1) << " V_Z: " << std::left << std::setw(18) << visc_term.get(2) << std::endl;
					// 	FileForce << "//////////////////////////////////////////////////////////////////////" << std::endl;
					// }

					// vd.getProp<drho>(a) += massb * (v_rel.get(0) * DW.get(0) + v_rel.get(1) * DW.get(1) + v_rel.get(2) * DW.get(2));
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

openfpm::vector<size_t> to_remove;

size_t cnt = 0;

void verlet_int(particles &vd, double dt)
{
	// list of the particle to remove
	to_remove.clear();

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
			// Update rho
			double rhop = vd.template getProp<rho>(a);

			// Update only the density
			vd.template getProp<velocity>(a)[0] = 0.0;
			vd.template getProp<velocity>(a)[1] = 0.0;
			vd.template getProp<velocity>(a)[2] = 0.0;
			double rhonew = vd.template getProp<rho_prev>(a) + dt2 * vd.template getProp<drho>(a);
			vd.template getProp<rho>(a) = (rhonew < rho_zero) ? rho_zero : rhonew;

			vd.template getProp<rho_prev>(a) = rhop;

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

	// // remove the particles
	// vd.remove(to_remove, 0);

	// increment the iteration counter
	cnt++;
}

void euler_int(particles &vd, double dt)
{
	// list of the particle to remove
	to_remove.clear();

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
			// Update rho
			double rhop = vd.template getProp<rho>(a);

			// Update only the density
			vd.template getProp<velocity>(a)[0] = 0.0;
			vd.template getProp<velocity>(a)[1] = 0.0;
			vd.template getProp<velocity>(a)[2] = 0.0;
			double rhonew = vd.template getProp<rho>(a) + dt * vd.template getProp<drho>(a);
			vd.template getProp<rho>(a) = (rhonew < rho_zero) ? rho_zero : rhonew;

			vd.template getProp<rho_prev>(a) = rhop;

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

	// // remove the particles
	// vd.remove(to_remove, 0);

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

	////////////////////////// WORKING CODE FOR THE OLD BC //////////////////////////////
	// Number of fluid and boundary particles in the x, y and z direction

	bool bc_new = false;
	// size_t Np_fluid[3] = {60, 1, 120};
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
	// offset to add to the fluid box, if non periodic it has to go from 0 to length
	double offset_fluid[3];
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
			offset_fluid[dim] = 0.0;
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
			offset_fluid[dim] = 0.0;
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

	Box<3, double> fluid_box({-offset_fluid[0],
							  -offset_fluid[1],
							  -offset_fluid[2]},
							 {length[0] + offset_fluid[0] + offset_periodic_fluid[0],
							  length[1] + offset_fluid[1] + offset_periodic_fluid[1],
							  length[2] + offset_fluid[2] + offset_periodic_fluid[2]});

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
	W_dap = 1.0 / Wab(H / 1.5);

	// extended boundary around the domain, and the processor domain
	Ghost<3, double> g(2.0 * H);

	particles vd(0, domain, bc, g, DEC_GRAN(512));

	// return an iterator to the fluid particles to add to vd
	auto fluid_it = DrawParticles::DrawBox(vd, sz, domain, fluid_box);

	// here we fill some of the constants needed by the simulation
	// h_swl = fluid_it.getBoxMargins().getHigh(2) - fluid_it.getBoxMargins().getLow(2);
	// B = (coeff_sound) * (coeff_sound)*gravity * h_swl * rho_zero / gamma_;
	// cbar = coeff_sound * sqrt(gravity * h_swl);

	max_fluid_height = fluid_it.getBoxMargins().getHigh(2);

	// Fill needed constants
	double L = length[0];

	double nu = sqrt(10.0 * visco * H * gravity * L * L / 64.0);
	// double nu = 0.01;
	// double visco = (64.0 / 10.0) * (nu * nu / (H * gravity * L));
	// double nu = eta;
	double umax = gravity * L * L / (8.0 * nu);
	cbar = coeff_sound * umax;
	B = cbar * cbar * rho_zero / gamma_;
	double Re = L * umax / nu;

	std::ofstream file("Constants.txt");
	file << "gravity: " << gravity << std::endl;
	file << "Re: " << Re << std::endl;
	file << "L: " << L << std::endl;
	file << "nu: " << nu << std::endl;
	file << "alpha: " << visco << std::endl;
	file << "umax: " << umax << std::endl;
	file << "cbar: " << cbar << std::endl;
	file << "B: " << B << std::endl;

	// for each particle inside the fluid box ...
	const double profile_parameter = gravity / (2.0 * nu);
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

		// No need for hydrostatic pressure in channel set up
		vd.template getLastProp<Pressure>() = 0.0; // rho_zero * gravity * (max_fluid_height - fluid_it.get().get(2));

		// vd.template getLastProp<Pressure>() = 0.0;

		// vd.template getLastProp<rho>() = pow(vd.template getLastProp<Pressure>() / B + 1, 1.0 / gamma_) * rho_zero;
		vd.template getLastProp<rho>() = rho_zero;
		vd.template getLastProp<rho_prev>() = vd.template getLastProp<rho>();
		vd.template getLastProp<velocity>()[0] = 0.0;
		vd.template getLastProp<velocity>()[1] = 0.0;
		// impose theoretical pouiselle flow profile for faster convergence
		vd.template getLastProp<velocity>()[2] = profile_parameter * fluid_it.get().get(0) * (L - fluid_it.get().get(0));

		vd.template getLastProp<velocity_prev>()[0] = 0.0;
		vd.template getLastProp<velocity_prev>()[1] = 0.0;
		vd.template getLastProp<velocity_prev>()[2] = profile_parameter * fluid_it.get().get(0) * (L - fluid_it.get().get(0));

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
		if (!recipient.isInside((position))) // with this we check if particle is outside the recipient box
		{
			++bound_box;
			continue;
		}

		if (fluid_box.isInside(position) && !(position.get(2) < z_end + epsilon && position.get(2) > z_end - epsilon))
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
		vd.template getLastProp<velocity>()[0] = 0.000;
		vd.template getLastProp<velocity>()[1] = 0.0;
		vd.template getLastProp<velocity>()[2] = 0.0;

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

	openfpm::vector<std::string> names({"type", "density", "density_prev", "pressure", "delta_density", "force", "velocity", "velocity_prev"});
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
		bool flag = 0;
		// if (t > 1.0 && t < 1.0005)
		// {
		// 	FileForce << "///////////////////////////////////////////" << std::endl;
		// 	FileForce << "NEW TIME" << std::endl;
		// 	FileForce << "///////////////////////////////////////////" << std::endl;
		// 	flag = 1;
		// }

		calc_forces(vd, NN, max_visc, flag);

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
			vd.map();
			vd.ghost_get<type, rho, Pressure, velocity>();
			vd.updateCellList(NN);

			// calculate the pressure at the sensor points
			// sensor_pressure(vd, NN, press_t, probes);

			// vd.deleteGhost();
			vd.write_frame(filename, write);
			// vd.ghost_get<type, rho, Pressure, velocity>();

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
