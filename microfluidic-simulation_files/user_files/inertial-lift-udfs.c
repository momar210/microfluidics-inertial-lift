/* UDF for computing the lateral displacement of particle suspended in fluid flow */


/**
 * References
 * Poiseuille Flow in Rectangular Channel Calculation (Umax)
 * 	https://web.archive.org/web/20220129140030/http://skinny.dk/picasso_fish/reports/Pedersen%20%282003%29%20Poiseuille%20Flow.pdf
 * Coefficient of Lift and Lift Force Calculations:
 * 	DOI: 10.1039/d1lc00225b
*/

#include "udf.h"
#include "mem.h" // cell indexing header
#include "dpm.h" // particle properties even though declared in the macro arguments
#include <math.h>

#define BETA(n,H) (2*n-1)*M_PI/H

/** TODO
 * How to pass in the powers and coefficients w/o hard coding
 * Share calculation of Pressure Gradient between UDF functions / files
 * How do we get the domain of the interior of the model?
 */

/** USER DEFINED MEMORY LOCATIONS
* 0: CLy
* 1: CLz
* 2: Max X velocity (m/s)
*/


void CALCULATE_LIFT_COEFFICIENTS(double AR, double Re, double kappa, double py, double pz, int index, double *CL)
{
  double inputs[5];
  double neuron_input_sum;
  int i, k;

  /*  ===== NEURAL NETWORK CONSTANTS ===== */
  double outputs_layer_1[20];
  static const double offsets_layer_1[20] = { -2.0664128335408196, -4.0637585987279294, 3.4178821508322792, -1.812401327910047, -0.70284484626917465,
    1.0312709391291706, -0.0512400314695292, -0.27224428449885008,
    -0.10595499530558948, 0.13798883929505978, -1.1073416447775288,
    1.2011137315220088, -0.59746206891260478, -0.48044767886307782,
    -0.36044784050974671, -1.4567495820092484, -1.1382956186553976,
    1.8542231535829585, 1.5565992288239583, 1.5054953405371005 };

  static const double weights_layer_1[100] = { 1.0113463104924822, 0.12605253090543531,
    0.015079993690685198, -0.087165422562620665, 0.070990002843449362,
    0.0022046946118252581, -0.92502021416812652, -0.11390083606583758,
    0.050632284902135095, 1.1772750924383637, -0.35952420136468893,
    0.793948851960803, 0.2873782084123912, 0.02617087943976729,
    0.098495344946317662, -1.0290738430468014, -0.18095309871089607,
    0.74162071188253287, 0.44537652446717224, 0.876249255186907,
    0.47865869190978771, 0.1030909692290797, -0.044768016855218015,
    0.53709872683483562, -0.0040048947414689732, -0.047831154406574171,
    -0.033978262551838785, 0.049655401685425489, -0.04777897750694763,
    -0.035928863838989439, -0.091112650663323858, 0.014400074092458736,
    0.085975981747414623, -0.39158462537518457, -0.22569980548290244,
    -0.24284234751418102, 0.098083648857335753, 0.023363972537760964,
    0.18008111962031295, -0.049466784226055535, 0.15013840358045971,
    0.06332281831028394, -0.025168750917782867, 0.058411428925736683,
    0.0973582404979998, -0.11727677229721047, 0.016711296228811498,
    -0.0608355346459892, 0.0069536721477733154, 0.25349194753193049,
    -0.00092930862712746172, -0.005498379628936065, -0.0021302947189261773,
    -0.21772678051382216, 0.26649322481007665, -0.0203631466496388,
    0.0972973128757761, -0.053510974162177913, -0.18517235181058825,
    -0.065799729886747815, 2.7590378669784865, -0.073329194457685723,
    -2.35424467798736, 1.1505122984614518, 0.6481210982159914, 1.23004880207753,
    -0.087479803194533648, -0.81367518041385123, 0.42870633735153735,
    -0.590107261194705, 0.4639250500694595, -0.35072224292067622,
    1.6099773842017262, 0.47096976583005062, 1.1515663993133791,
    -0.48598307383015055, -1.0362655968288494, -0.33222958483093462,
    -0.073143199502663778, 0.310904528359796, -1.3398938878203865,
    3.0308596027687575, 0.84137531907624208, 0.1855295612604155,
    0.39298668935200443, 0.3851229617574195, 1.3427019554609831,
    0.39823355471461347, -0.52180797895240927, 0.30452280476983429,
    0.20365956630006024, -0.93398693917729358, -0.53973166594524324,
    0.41153685152410419, -0.089602453460719142, -1.5163066019789919,
    -0.52624079462359619, 0.27657185166625264, -1.5501706358502383,
    -0.34644147911952977 };

  double outputs_layer_2[8];
  static const double offsets_layer_2[8] = { 1.2616439401058175, -0.88593218150199315,
    0.026411830150108496, -0.3321051287584475, -1.3888890591282435,
    -0.46068698488192128, -2.5623891336545705, 2.1171375438762543 };

  static const double weights_layer_2[160] = { 0.026873055814601709, -0.51088075200209948,
    0.2946786324235261, 0.56607976157569306, 0.0013050149630168875,
    -0.28075281088385151, -0.03744556124840559, -0.030386917926862163,
    -2.2585845223453265, -1.1537273480454431, -2.10613726494928,
    1.5980776508627674, -2.9526522004093971, 2.0216008699340455,
    -0.54232526534323144, -0.550247065646972, -4.6778207426963725,
    -0.72962693288450509, -1.9045950848145279, 0.34751148438823509,
    -0.027155184854069297, 3.1487245464402567, 0.19352024299795259,
    0.43784156537289154, -0.37281384143751273, 0.0848961463951636,
    0.011742826880873528, -0.051737091952646828, 0.65983833292636018,
    0.14609659934386149, -0.662229628724031, -0.7594965995982732,
    -1.4439823362383912, 1.2187937322356635, -0.9520950071078792,
    -2.3944518365757359, 2.0152714780384571, 1.4051089821108849,
    -1.2026066175843952, 1.4874363746946484, 0.412144976722042,
    -0.8102098663170223, 0.82513459423659252, 0.83546466302365607,
    -0.639092475827621, -0.84050300076656748, 0.28524794497328471,
    -3.9206927115325576, -0.41989179753915717, 0.548564755891428,
    -0.28983153038505893, -0.27266956963131694, -0.029807482091859502,
    0.29262680152826176, -0.10701179789880781, -0.18896408585194807,
    -1.3007978205391866, 1.3178854229828996, -0.17742818045333325,
    -1.5414548303614168, 3.4990753913887831, 0.736430516182038,
    -1.547064664126576, -0.39564631182206927, -0.90074606766912291,
    1.9801444542394835, -0.34221618452683988, -2.4539751941380978,
    3.2770515744737141, 1.1573817859016045, -2.2630366434169225,
    0.86066754916955412, 0.10736393782438061, 0.01140245992163645,
    -0.10865284122299239, -0.14649367541649297, -0.048343141312874524,
    0.13196520980945478, -0.10253991641833553, -0.2484611768972913,
    1.1666860145082445, -4.80770726487812, 2.4410656986365216,
    1.4760684079209629, -0.034656622816495827, -2.4599186423743356,
    0.46090572970068971, -1.959634136483541, -0.16378755597531139,
    -0.87449186744459673, 1.0537360459299836, -0.17615865852552698,
    -0.0538752379729844, -0.96118379135908538, -0.22421924375597774,
    -1.0853458829529981, 0.45952386572059273, -0.26435251194004766,
    0.31967188472434083, 0.59620316504338156, -0.26581062146654727,
    -0.50926898171511426, 0.65513751035249546, -0.43101349173666031,
    -0.61832460781117993, 1.06823904258369, -0.15099507230824613,
    -0.7383688536836871, 0.55127595137790109, 0.317856571592111,
    -0.64490302385358078, -0.32369003397368845, -0.061157611678667193,
    0.14561739679645069, 0.277630106421263, 0.29420133648001434,
    0.934329235501248, -0.27885347177744485, -0.035067191648449327,
    -1.0717415416525589, -0.13412892885085856, 0.31451253242911609,
    -0.10140637969328634, -0.10562378102289037, 0.015475238032569318,
    0.12859489183651432, -0.083651286898072275, -0.13936349977001447,
    0.90261850563337465, -1.203580409660876, 0.91127930903663334,
    0.778362952033977, -0.777227525409743, -0.93635691544673061,
    0.51228395878032473, -2.4069778086459435, 2.1273200360714171,
    -3.08351812039097, 1.6661077618028863, 2.5522445246390246,
    -0.59170709264117594, -1.7471908834720251, 0.04057872944602129,
    -2.2180169415997728, 0.25432714319257432, 0.43882180100684048,
    0.029628644628008832, -0.47369033105330932, 0.43887142231107057,
    -0.0043837188251823608, 0.019919644638033321, -0.72802465075228628,
    1.4158413404121586, -1.1331244764501247, 0.25613735015458139,
    0.79045540117727453, -0.078688633326299773, -0.46929067065558755,
    1.1074583801299289, 1.5586937017540237 };

  static const double weights_layer_3[16] = { 1.1244585306625683, 6.6502524293410294,
    0.83466997432988643, 3.1357642947377786, -3.8723989430019787,
    2.6895736638231669, -0.2533418798893467, 1.801942099242462,
    -2.18942031900452, -0.09338083899079859, -4.1737830182823359,
    3.1507813045844046, -1.0980672170047456, 0.16778083296857424,
    1.8564500087860198, -0.3300098884108103 };

  /*  ===== NORMALIZATION CONSTANTS ===== */
  /*  Input Normalization */
  inputs[0] = (AR - 1.0) * 0.33333333333333331 * 2.0 - 1.0;
  inputs[1] = (Re - 50.0) * 0.0066666666666666671 * 2.0 - 1.0;
  inputs[2] = (kappa - 0.1) * 5.0 * 2.0 - 1.0;
  inputs[3] = fabs(pz) * 1.25 * 2.0 - 1.0;
  inputs[4] = fabs(py) * 1.0526315789473679 * 2.0 - 1.0;

  /*  Layer 1 (tansig layer) */
  for (k = 0; k < 20; k++) {
    neuron_input_sum = 0.0;
    for (i = 0; i < 5; i++) {
      neuron_input_sum += weights_layer_1[k + 20 * i] * inputs[i];
    }

    outputs_layer_1[k] = tanh(offsets_layer_1[k] + neuron_input_sum);
  }

  /*  Layer 2 (tansig layer) */
  for (k = 0; k < 8; k++) {
    neuron_input_sum = 0.0;
    for (i = 0; i < 20; i++) {
      neuron_input_sum += weights_layer_2[k + (i << 3)] * outputs_layer_1[i];
    }

    outputs_layer_2[k] = tanh(offsets_layer_2[k] + neuron_input_sum);
  }

  /*  Layer 3 (linear layer) and Output Layer */
  neuron_input_sum = 0.0;
  for (k = 0; k < 8; k++) {
    neuron_input_sum += weights_layer_3[(index-1) + (k << 1)] * outputs_layer_2[k];
  }

  /* Output De-normalization */
  *CL = (((-1.3400483104042133 * (double)(index-1) + -2.9103877624314474) + neuron_input_sum) + 1.0)
    / 2.0 / (0.093013622399621965 * (double)(index-1) + 0.99451030312674) +
    (0.085999999999999965 * (double)(index-1) + -0.7776);
}



double CALCULATE_PRESSURE_GRADIENT(double W, double H, double mu, double Q)
{
  // DECLARATION OF VARIABLES
  double dPdx, s, k;
  int n;

  // calculate the pressure gradient at this point given the flow rate
  s = 0;
  for (n = 1; n <= 5; n++) {
      s += (1/pow(2*n-1,5))* (cosh(BETA(n,H) * W) - 1) / sinh(BETA(n,H) * W);
  }
  k = ((pow(H,3))*W)/(12*mu) - (16*pow(H,4))/(pow(M_PI,5)*mu) * s;
  dPdx = Q/k; // Pressure gradient

  return dPdx;
}


/***
 *  Solution to PDEs of 3D Poiseuille Flow in a Rectangular Channel
 *  Boussinesq derivation of flow velocity in rectangular channel
 *  
 * 	arguments:
 * 	double w: channel width in meters
 *  double h: channel height in meters
 *  double Q: volumetric flow rate in the channel
 */
double CALCULATE_MAX_VELOCITY(double H, double W, double mu, double dPdx)
{
	// DECLARATION OF VARIABLES
	double Umax, s;
	int n;

	s = 0;
	for (n = 1; n <= 5; n++) {
	    s += (1/pow(2*n-1,3)) * ((2*sinh(BETA(n,H) * W / 2)) / sinh(BETA(n,H) * W));
	}

	Umax = (dPdx/(2*mu))*(pow(H/2,2)) - ((4*dPdx*pow(H,2))/(mu*pow(M_PI,3))) * s;

	// return the velocity in the x direction
	return Umax;
}


double CALCULATE_VELOCITY(double H, double W, double y, double z, double mu, double dPdx)
{
	// DECLARATION OF VARIABLES
	double unm, ux;
	int n, m, maxiter;

	unm = 0;
	maxiter = 20;
	maxiter = (maxiter + (1 - maxiter % 2));
	for (n = 1; n <= maxiter; n += 2) {
		for (m = 1; m <= maxiter; m += 2) {
			unm += (1 / (n * m * (pow(n, 2) / pow(W, 2) + pow(m, 2) / pow(H, 2)))) * sin(n * (M_PI / W) * z) * sin(m * (M_PI / H) * y);
		}
	}
	ux = (16 / pow(M_PI, 4)) * (dPdx / mu) * unm;

	return ux;
}


void CALCULATE_CHANNEL_PARAMETERS(double particle_diameter, double mu, double rho, double *H, double *W, double *Re, double *kappa, double particle_position[3])
{
	double x = particle_position[0];
	double channel_center_offset_y = 0;
	double channel_center_offset_z = 0;
	double A, P, Q, Dh;

	Q = Get_Input_Parameter("volumetric_flow_rate_ul_per_min"); // uL/min
	Q = Q * 1e-9 / 60; // convert Q to m^3/s

	/* THIS IS HOW YOU DESCRIBE THE CHANNEL AR */
	double notch_length = Get_Input_Parameter("notch_length");
	double notch_midpoint = Get_Input_Parameter("notch_midpoint");

	if (x < (notch_midpoint - notch_length/2) || x > (notch_midpoint + notch_length/2)) { // we are in the main channel
		// no need to change the channel center location

		*W = 80e-6;
		*H = 40e-6;
	} else {
		channel_center_offset_y = 10e-6;
		channel_center_offset_z = 0;

		*W = 80e-6;
		*H = 20e-6;
	}

	/* normalize the particle position in y and z with respect to the center position of the channel */
	particle_position[1] = 2 * (particle_position[1] - channel_center_offset_y) / *H;
	particle_position[2] = 2 * (particle_position[2] - channel_center_offset_z) / *W;

	A = *W * *H; // channel area
	P = *W * 2 + *H * 2; // channel perimeter

	Dh = 4 * A / P; // hydraulic
	*Re = rho * Q * Dh / (mu * A);
	*kappa = particle_diameter / *H; // blockage ratio
}


DEFINE_DPM_BODY_FORCE(inertial_lift, p, i)
{
  if (i == 0) {
    return 0; // do not calculate lift force in the X direction
  }

	// DECLARATION OF VARIABLES
	double CL, FL, particle_diameter, H, W, Re, kappa, mu, rho, Umax;
	cell_t c = P_CELL(p); // the cell in which the particle is present
	Thread *t = P_CELL_THREAD(p); // thread initialization

	double *particle_position = P_POS(p);

	particle_diameter = P_DIAM(p); 		// Particle diameter
	mu = C_MU_L(c,t); 	// Cell dynamic viscosity
	rho = C_R(c,t); 	// Cell density

	CALCULATE_CHANNEL_PARAMETERS(particle_diameter, mu, rho, &H, &W, &Re, &kappa, particle_position);

	CALCULATE_LIFT_COEFFICIENTS(W/H, Re, kappa, particle_position[1], particle_position[2], i, &CL);

	// Maximum velocity for each cell should be precalculated prior to starting DPM calculation
	Umax = C_UDMI(c, t, 2);

	// Calculate the lift force based on Su et al 2021 Equation 7
	FL = CL * rho * pow(Umax, 2) * pow(particle_diameter, 4) / pow(H, 2);

	// An acceleration should be returned
	return (FL / P_MASS(p));
}


DEFINE_PROFILE(inlet_x_velocity, thread, position)
{
	// DECLARATION OF VARIABLES
	double face_centroid_pos[3]; // this will hold the position vector
	double dPdx, W, H, y, z, mu, Q;
	face_t f;

	W = Get_Input_Parameter("channel_width");
	H = Get_Input_Parameter("channel_height");

	// TODO: would be nice to not have to hard code this
	mu = 0.001003; // kg/(m*s)
	Q = Get_Input_Parameter("volumetric_flow_rate_ul_per_min"); // uL/min
	Q = Q * 1e-9 / 60; // convert Q to m^3/s

	dPdx = CALCULATE_PRESSURE_GRADIENT(H, W, mu, Q);

	begin_f_loop(f, thread)
	{
		F_CENTROID(face_centroid_pos, f, thread);
		
		y = face_centroid_pos[1] + 0.5 * H;
		z = face_centroid_pos[2] + 0.5 * W;
		F_PROFILE(f, thread, position) = CALCULATE_VELOCITY(H, W, y, z, mu, dPdx);
	}
	end_f_loop(f, thread)

	Message("Finished Loading Profile\n");
}


// TODO: not sure why this happens more than once
DEFINE_ON_DEMAND(on_demand_max_velocity_calculation)
{
	// DECLARATION OF VARIABLES
	double cell_centroid_position[ND_ND]; /* this will hold the position vector */
	double dPdx, W, H, x, y, z, mu, Q, Umax, particle_diameter, rho, Re, kappa, CL;
	Domain *domain;
	face_t f;
	cell_t c;
	Thread *t;

	Q = Get_Input_Parameter("volumetric_flow_rate_ul_per_min"); // uL/min
	Message("Volumetric Flow Rate (User Parameter): %.2f [uL/min]\n", Q);
	Q = Q * 1e-9 / 60; // convert Q to m^3/s

  particle_diameter = Get_Input_Parameter("particle_diameter");
  Message("Particle Diameter (User Parameter): %.2f [um]\n", particle_diameter * 1e6);

	domain=Get_Domain(1);

	thread_loop_c(t, domain)
	{
		begin_c_loop(c, t)
		{
			C_CENTROID(cell_centroid_position, c, t); // get the cell position
			mu = C_MU_L(c,t); 	// Cell dynamic viscosity
			rho = C_R(c,t); 	// Cell density

			CALCULATE_CHANNEL_PARAMETERS(particle_diameter, mu, rho, &H, &W, &Re, &kappa, cell_centroid_position);

      // calculate the coefficient of lift for a particle at the cell center
      CALCULATE_LIFT_COEFFICIENTS(W/H, Re, kappa, cell_centroid_position[1], cell_centroid_position[2], 1, &CL);
      C_UDMI(c, t, 0) = CL;
      CALCULATE_LIFT_COEFFICIENTS(W/H, Re, kappa, cell_centroid_position[1], cell_centroid_position[2], 2, &CL);
      C_UDMI(c, t, 1) = CL;

			// from the height and width of the channel, determine what the pressure graident and maximum velocity are
			dPdx = CALCULATE_PRESSURE_GRADIENT(H, W, mu, Q);
			Umax = CALCULATE_MAX_VELOCITY(H, W, mu, dPdx);

			// save the maximum velocity of the channel at this location
			C_UDMI(c, t, 2) = Umax;
		}
		end_c_loop(c, t)
	}

	Message("Finished Setting Umax for all Cells in Domain\n");
}