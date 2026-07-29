#include "choice_solution.hpp"

Choice_solution::Choice_solution(Parameters & par){

	/**
	 * @details
	 * @param[in] par contains all the values from the parameter
	 * @warning Error: the dimension is ***
	 * @warning This *** solution for L=*** does not exist!
	 * @warning *** solutions for the domain *** do not exist!
	 * @note If the solution does not exists, the code will exit with failure termination code.
	 * @todo Exceptions should be treated.
	 */

	dim2 = 2.*par.get_choicedim()-1.;
	if (abs(2*par.get_choicedim()-1 - dim2)>EPSILON){ // if 2*choicedim is not an integer
		cerr << "Error: the dimension is " << par.get_choicedim() << endl;
		exit(EXIT_FAILURE);
	}
	switch (dim2){
		case 1: // dimension 1
			switch (par.get_choicetype()){
				case 0:
					/******************************************************************************
					 * Solutions over an inclined plane
					 *  see inclined_plane* files for details (dimension, water heights, ...)
					 ******************************************************************************/
					switch (par.get_choicedomain()){
						case 1:
							/******************************************************************************
							 * L=10 m
							 * 1: supercritical flow
							 ******************************************************************************/
							switch (par.get_choice()){
								case 1:
									sol = new Inclined_plane(par);
									break;

								default:
									cerr<< "This solution over an inclined plane for L=10 m does not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
							break;

						case 2 :
							/********************************************************************************
							 * L=20 m
							 * 1: transient solution
							 ********************************************************************************/
							switch(par.get_choice()){
								case 1:
									sol = new Swash(par);
									break;

								case 2:
									sol=new Swash(par) ;
									break ;

								default:
									cerr<< "This solution over an inclined plane for L=20 m does not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
							break;

						default:
							cerr<< "Solutions over an inclined plane for the domain "<<par.get_choicedomain()<<" do not exist!"<< endl;
							exit(EXIT_FAILURE);

					}
					break;


				case 1:
					/******************************************************************************
					 * Bump solutions
					 *  see bump.* files for details (dimension, water heights, ...)
					 ******************************************************************************/
					switch (par.get_choicedomain()){
						case 1:
							/******************************************************************************
							 * L=25 m
							 * 1: subcritical flow
							 * 2: transcritical without shock (subcritical-supercritical)
							 * 3: transcritical with shock (subcritical-supercritical-subcritical)
							 * 4: lake at rest with an immersed bump
							 * 5: lake at rest with an emerged bump
							 ******************************************************************************/

							switch (par.get_choice()){
								case 1:
									sol = new Bump(par);
									break;

								case 2:
									sol = new Bump(par);
									break;

								case 3:
									sol = new Bump(par);
									break;

								case 4:
									sol = new Bump(par);
									break;

								case 5:
									sol = new Bump(par);
									break;

								default:
									cerr<< "This bump solution for L=25 m does not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
							break;

						default:
							cerr<< "Bump solutions for the domain "<<par.get_choicedomain()<<" do not exist!"<< endl;
							exit(EXIT_FAILURE);

					}
					break;

				case 2:
					/******************************************************************************
					 * MacDonald's like solutions
					 *  see MacDonald_like.* files for details (dimension, water heights, ...)
					 ******************************************************************************/
					switch (par.get_choicedomain()){
						case 1:
							/******************************************************************************
							 * Long channel : L=1000 m
							 * 1: subcritical flow with Darcy-Weisbach law
							 * 2: subcritical flow with Manning law
							 * 3: supercritical flow with Darcy-Weisbach law
							 * 4: supercritical flow with Manning law
							 * 5: sub- to super-critical flow with Darcy-Weisbach law
							 * 6: sub- to super-critical flow with Manning law
							 * 7: super- to sub-critical flow with Darcy-Weisbach law
							 * 8: super- to sub-critical flow with Manning law
							 ******************************************************************************/
							switch (par.get_choice()){
								case 1:
									sol = new MacDonald_like(par);
									break;

								case 2:
									sol = new MacDonald_like(par);
									break;

								case 3:
									sol = new MacDonald_like(par);
									break;

								case 4:
									sol = new MacDonald_like(par);
									break;

								case 5:
									sol = new MacDonald_like(par);
									break;

								case 6:
									sol = new MacDonald_like(par);
									break;

								case 7:
									sol = new MacDonald_like(par);
									break;

								case 8:
									sol = new MacDonald_like(par);
									break;

								default:
									cerr<< "This MacDonald solution for L=1000 m does not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
							break;

					 case 2:
							/******************************************************************************
							 * Short channel : L=100 m
							 * 2: supercritical-supercritical with Manning law
							 * 4: subcritical-subcritical with Manning law
							 * 6: subcritical-supercritical with Manning law
							 ******************************************************************************/
							switch (par.get_choice()){
								case 2:
									sol = new MacDonald_like(par);
									break;

								case 4:
									sol = new MacDonald_like(par);
									break;

								case 6:
									sol = new MacDonald_like(par);
									break;

								default:
									cerr<< "This MacDonald solution for L=100 m does not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
							break;

						case 3:
							/******************************************************************************
							 * Long, undulating, periodic channel : L=5000 m
							 * 2: subcritical-subcritical with Manning law
							 ******************************************************************************/
							switch (par.get_choice()){
								case 2:
									sol = new MacDonald_like(par);
									break;

								default:
									cerr<< "This MacDonald solution for an undulating channel with L=5000 m does not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
							break;

						case 4:
							/******************************************************************************
							 * Long channel : L=1000 m with rain
							 * 1: subcritical-subcritical with rain with Darcy-Weisbach law
							 * 2: subcritical-subcritical with rain with Manning law
							 * 3: supercritical-supercritical with rain with Darcy-Weisbach law
							 * 4: supercritical-supercritical with rain with Manning law
							 ******************************************************************************/
							switch (par.get_choice()){
								case 1:
									sol = new MacDonald_like(par);
									break;

								case 2:
									sol = new MacDonald_like(par);
									break;

								case 3:
									sol = new MacDonald_like(par);
									break;

								case 4:
									sol = new MacDonald_like(par);
									break;

								default:
									cerr<< "This MacDonald solution for L=1000 m with rain does not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
							break;

						case 5:
							/******************************************************************************
							 * Long channel : L=1000 m with diffusion
							 * 1: subcritical-subcritical with diffusion
							 * 2: supercritical-supercritical with diffusion
							 ******************************************************************************/
							switch (par.get_choice()){

								case 1:
									sol = new MacDonald_like_diffus(par);
									break;

								case 2:
									sol = new MacDonald_like_diffus(par);
									break;

								default:
									cerr<< "This MacDonald solution for L=1000 m with diffusion does not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
							break;

						default:
							cerr<< "MacDonald solutions for the domain "<< par.get_choicedomain() << " do not exist!"<< endl;
							exit(EXIT_FAILURE);
						}
					break;

				case 3:
					/******************************************************************************
					 * Dam break solutions
					 *  see *dam_break.*, Dressler_dam.* files for details (dimension, water heights, ...)
					 ******************************************************************************/
					switch (par.get_choicedomain()){
						case 1:
							/******************************************************************************
							 * L=10 m
							 * 1: Dam break on a wet domain without friction (Stoker's solution)
							 * 2: Dam break on a dry domain without friction (Ritter's solution)
							 * 3: Dam break on a dry domain with friction (Dressler's solution)
							 ******************************************************************************/

							switch (par.get_choice()){
								case 3:
									sol = new DresslerDam(par);
									break;

								default:
									cerr<< "This dam break solution for L=10 m does not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
							break;

						case 2:
							/******************************************************************************
							 * L=20 m
							 * 1: Self-similar dam break on a flat bottom with a laminar friction
							 * 2: Self-similar dam break on an inclined plane with a laminar friction
							 ******************************************************************************/

							switch (par.get_choice()){
								case 1:
									sol = new Selfsimilar_dam_break(par);
									break;

								case 2:
									sol = new Selfsimilar_dam_break(par);
									break;

								default:
									cerr<< "This dam break solution for L=20 m does not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
							break;


						default:
							cerr<< "Dam break solutions for the domain "<< par.get_choicedomain() << " do not exist!"<< endl;
							exit(EXIT_FAILURE);
						}
					break;

				case 4:
					/******************************************************************************
					 * Oscillations
					 *  see Thacker1D.* and Sampson.* files for details (dimension, water heights, ...)
					 ******************************************************************************/
					switch (par.get_choicedomain()){
						case 2:
							/******************************************************************************
							 * L=10000 m
							 * 1: Planar surface in a parabola with a linear friction (Sampson's solution)
							 ******************************************************************************/

							switch (par.get_choice()){
								case 1:
									sol = new Sampson(par);
									break;

								default:
									cerr<< "This oscillation solution for L=10000 m does not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
							break;

						default:
							cerr<< "Oscillation solutions for the domain "<< par.get_choicedomain() << " do not exist!"<< endl;
							exit(EXIT_FAILURE);
					}
					break;

				case 5:
					/******************************************************************************
					 * Bedload transport (Exner equation)
					 *  see bedload.* file for details (dimension, water heights, ...)
					 ******************************************************************************/
					switch (par.get_choicedomain()){
						case 1:

							/******************************************************************************
							 * L= 15 m
							 * 1: Grass equation
							 * 2: Meyer-Peter & Muler equation
							 ******************************************************************************/

							switch (par.get_choice()){
								case 1:
									sol = new Bedload(par);
									break;

								case 2:
									sol = new Bedload(par);
									break;

								default:
									cerr<< "This beadload solution for L=15 m does not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
							break;

						default:
							cerr<< "Bedload solutions for the domain "<< par.get_choicedomain() << " do not exist!"<< endl;
							exit(EXIT_FAILURE);
					}
					break;

				case 6:
					/******************************************************************************
					 * Sluice gate
					 *  see sluice_gate.* file for details (dimension, water heights, ...)
					 ******************************************************************************/
					switch (par.get_choicedomain()) {
					case 1:

						/******************************************************************************
						 * L= 10 m
						 * 1: Sluice gate opening on a dry domain 
						 * 2: Sluice gate opening on a wet domain with free flow and low h_right = 0.01 * gate_size
						 * 3: Sluice gate opening on a wet domain with free flow and h_right=gate_size
						 ******************************************************************************/

						switch (par.get_choice()) {
						case 1:
							sol = new Sluice_gate(par);
							break;

						case 2:
							sol = new Sluice_gate(par);
							break;

						case 3:
							sol = new Sluice_gate(par);
							break;

						default:
							cerr << "This sluice gate solution for L=10 m does not exist!" << endl;
							exit(EXIT_FAILURE);
						}
						break;

					default:
						cerr << "Sluice gates solutions for the domain " << par.get_choicedomain() << " do not exist!" << endl;
						exit(EXIT_FAILURE);
					}
					break;
					
				case 8:
					/******************************************************************************
					 * Solute problem 
					 * see solute.* file for details (dimension, configuration, ...)
					 ******************************************************************************/
					switch (par.get_choicedomain()) {
					case 1:

						/******************************************************************************
						 * L=1000 m
						 * 1: no degradation with initial dissolved concentration
						 * 2: no degradation with boundary dissolved concentration
						 * 3: degradation with initial dissolved concentration
						 * 4: degradation with boundary dissolved concentration
						 ******************************************************************************/

						switch (par.get_choice()) {
						case 1:
							sol = new Solute(par);
							break;
						case 2:
							sol = new Solute(par);
							break;
						case 3:
							sol = new Solute(par);
							break;
						case 4:
							sol = new Solute(par);
							break;
							
						default:
							cerr << "This solute solution for L=1000 m does not exist!" << endl;
							exit(EXIT_FAILURE);
						}
						break;

					default:
						cerr << "Solute solutions for the domain " << par.get_choicedomain() << " do not exist!" << endl;
						exit(EXIT_FAILURE);
					}
					break;
					
				case 9:
					/******************************************************************************
					 * Mobile rain problem
					 * See files rain.* for details (dimension, configuration...)
					 ******************************************************************************/
					switch (par.get_choicedomain()) {
					case 1:

						/******************************************************************************
						 * L=18000 m
						 * 1: rain with the same velocity as the flow
						 * 2: rain with a velocity smaller than the flow
						 * 3: rain with a velocity larger than the flow
						 ******************************************************************************/

						switch (par.get_choice()) {
						case 1:
							sol = new Rain(par);
							break;
						case 2:
							sol = new Rain(par);
							break;
						case 3:
							sol = new Rain(par);
							break; 
							
						default:
							cerr << "This mobile rain solution for L=18000 m does not exist!" << endl;
							exit(EXIT_FAILURE);
						}
						break;

					default:
						cerr << "Mobile rain solutions for the domain " << par.get_choicedomain() << " do not exist!" << endl;
						exit(EXIT_FAILURE);
					}
					break;

				default:
					cerr<< "This type of solutions in one dimension does not exist!"<< endl;
					exit(EXIT_FAILURE);

			}

		break;

		case 2: // pseudo2D
			switch (par.get_choicetype()){
					case 1:
							/******************************************************************************
							 * MacDonald PSEUDO 2D
							 *  see MacDonald.* files for details (dimension, water heights, ...)
							 ******************************************************************************/
							switch (par.get_choicedomain()){
								case 1:
									/******************************************************************************
									 * Rectangular short channel, shape B1: L=200 m
									 * 1: subcritical flow
									 * 2: supercritical flow
									 * 3: smooth transition
									 * 4: hydraulic jump
									 ******************************************************************************/

									switch (par.get_choice()){
										case 1:
											sol = new MacDonaldB1(par);
											break;
										case 2:
											sol = new MacDonaldB1(par);
											break;
										case 3:
											sol = new MacDonaldB1(par);
											break;
										case 4:
											sol = new MacDonaldB1(par);
											break;
										default:
											cerr<< "This MacDonald solution for a rectangular B1 channel L=200 m does not exist!"<< endl;
											exit(EXIT_FAILURE);
									}
									break;
								case 2 :
									/******************************************************************************
									 * Trapezoidal long channel, shape B2: L=400 m
									 * 1: subcritical flow
									 * 2: smooth transition and hydraulic jump
									  ******************************************************************************/
									switch (par.get_choice()){
										case 1:
											sol = new MacDonaldB2(par);
											break;
										case 2:
											sol = new MacDonaldB2(par);
											break;
										default:
											cerr<< "This MacDonald solution for a trapezoidal B2 channel L=400 m does not exist!"<< endl;
											exit(EXIT_FAILURE);
									}
									break;

								default:
									cerr<< "MacDonald pseudo-2D solutions for the domain "<< par.get_choicedomain() << " do not exist!"<< endl;
									exit(EXIT_FAILURE);
							}
					break;

				default:
					cerr<< "This type of pseudo2D solutions does not exist!"<< endl;
					exit(EXIT_FAILURE);

			}

			break;
		case 3: // 2D
			switch (par.get_choicetype()){
				case 2:
					/******************************************************************************
					 * Static dam
					 *  see dam_2D.* files for details (dimension, water heights, ...)
					 ******************************************************************************/
					switch (par.get_choicedomain()) {
					case 1:

						/******************************************************************************
						 * L=25 m     l=10 m    Parabolic dam
						 ******************************************************************************/

						switch (par.get_choice()) {
						case 1:
							sol = new Dam_2D(par);
							break;


						default:
							cerr << "This 2D Dam solution for L=25 m does not exist!" << endl;
							exit(EXIT_FAILURE);
						}
						break;

					case 2:

						/******************************************************************************
						 * L=l=10 m      Cross like dam
						 ******************************************************************************/
						switch (par.get_choice()) {
						case 1:
							sol = new Dam_2D(par);
							break;


						default:
							cerr << "This 2D Dam solution for L=10 m does not exist!" << endl;
							exit(EXIT_FAILURE);
						}
						break;

					default:
						cerr << "Dam solution for the domain " << par.get_choicedomain() << " do not exist!" << endl;
						exit(EXIT_FAILURE);
					}
					break;

				case 3:
					/******************************************************************************
					 * Spherical geometry
					 *  see spherical.* files for details (dimension, geometry, ...)
					 ******************************************************************************/
					switch (par.get_choicedomain()) {
					case 1:

						/******************************************************************************
						 * radius=6 371 220 m  omega=7.292*10^-5 s^-1   alpha=0 radiants
						 ******************************************************************************/

						switch (par.get_choice()) {
						case 1:
							sol = new Spherical(par);
							break;

						default:
							cerr << "This Spherical solution for this domain does not exist!" << endl;
							exit(EXIT_FAILURE);
						}
						break;

					case 2:

						/******************************************************************************
						 * radius=6 371 220 m  omega=7.292*10^-5 s^-1  alpha=0.406 radiants
						 ******************************************************************************/

						switch (par.get_choice()) {
						case 1:
							sol = new Spherical(par);
							break;

						default:
							cerr << "This Spherical solution for this domain does not exist!" << endl;
							exit(EXIT_FAILURE);
						}
						break;


					default:
						cerr << "Spherical solution for the domain " << par.get_choicedomain() << " do not exist!" << endl;
						exit(EXIT_FAILURE);
					}
					break;

				default:
					cerr << "This type of solutions in two dimensions does not exist!" << endl;
					exit(EXIT_FAILURE);

			}

			break;


		default:
		cerr<< "Dimension "<< par.get_choicedim()<<" impossible!!"<< endl;
		exit(EXIT_FAILURE);
	}
}


void Choice_solution::compute(){
	sol->compute();
}

Choice_solution::~Choice_solution(){
	if (sol != NULL){
		delete sol;
		sol = NULL;
	}
}
