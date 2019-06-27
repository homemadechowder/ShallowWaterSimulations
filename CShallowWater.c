#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#define RES 50
#define HALO 2

double gravity =  9.80616;

/*struct myArray{
	double *data;//data of array
	int *size;//amount of elements for each dimension size[i]
}*/

void GaussianFcn(double x, double y, double *h, double *u, double *v)
{
  double a;

  /*  Radius from the center of the hill */
  a = sqrt(x * x + y * y);

  /*  Half-width of the hill */
  /*  Fluid height perturbation */
  /*  Height of the hill at this point */
  *h = 500.0 + 5.0 * exp(-(a * a) / 1.6E+11);

  /*  Zero velocity */
  *u = 0.0;
  *v = 0.0;
}

void GeostrophicFcn(double x, double y, double *h, double *u, double *v, double omega)
{
  double h_tmp;
  (void)x;

  /*  Background height */
  /*  Background velocity */
  /*  Latitude (override) */
  /*  Coriolis parameter (override) */
  /*  Wavelength */
  /*  Height of the hill at this point */
  h_tmp = 6.2831853071795862 * y / 1.0E+6;
  *h = 500.0 + 30.0 * (2.0 * omega * 0.70710678118654746) / gravity * 1.0E+6 /
    6.2831853071795862 * cos(h_tmp);

  /*  Fixed zonal velocity */
  *u = 30.0 * sin(h_tmp);

  /*  Zero meridional velocity */
  *v = 0.0;
}

void GradientFcn(double x, double y, double *h, double *u, double *v)
{
  double r;
  double theta;
  double uphi;
  double h_tmp;

  /*  Radius from the center of the hill */
  r = sqrt(x * x + y * y);

  /*  Angle from the center of the hill */
  theta = atan2(y, x);

  /*  Half-width of the hill */
  /*  Fluid height perturbation */
  /*  Background fluid height */
  /*  Calculate fluid height */
  uphi = r * r;
  h_tmp = exp(-uphi / 4.0E+10);
  *h = 500.0 - 20.0 * h_tmp;

  /*  Rotational velocity */
  uphi = sqrt(gravity * 20.0 * (2.0 * uphi / 4.0E+10) * h_tmp);
  if (r == 0.0) {
    *u = 0.0;
    *v = 0.0;
  } else {
    *u = -sin(theta) * uphi;
    *v = cos(theta) * uphi;
  }
}

void ApplyBoundaryConditions(double h[RES + 2*HALO][RES + 2*HALO], double hu[RES + 2*HALO][RES + 2*HALO], double hv[RES + 2*HALO][RES + 2*HALO]){
	int i;
	int j;
	for(j = 0; j < RES + 2*HALO; j++){
		h[0][j] = h[RES][j];
		h[1][j] = h[RES + 1][j];
		h[RES + HALO][j] = h[HALO][j];
		h[(RES + 2 * HALO) - 1][j] = h[2 * HALO - 1][j];

		hu[0][j] = hu[RES][j];
		hu[1][j] = hu[RES + 1][j];
		hu[RES + HALO][j] = hu[HALO][j];
		hu[(RES + 2 * HALO) - 1][j] = hu[2 * HALO - 1][j];

		hv[0][j] = hv[RES][j];
		hv[1][j] = hv[RES + 1][j];
		hv[RES + HALO][j] = hv[HALO][j];
		hv[(RES + 2 * HALO) - 1][j] = hv[2 * HALO - 1][j];
	}
	for(i = 0; i < RES + 2*HALO; i++){
		h[i][0] = h[0][RES];
		h[i][1] = h[1][RES + 1];
		h[i][RES + HALO] = h[i][HALO];
		h[i][(RES + 2 * HALO) - 1] = h[i][2 * HALO - 1];

		hu[i][0] = hu[0][RES];
		hu[i][1] = hu[1][RES + 1];
		hu[i][RES + HALO] = hu[i][HALO];
		hu[i][(RES + 2 * HALO) - 1] = hu[i][2 * HALO - 1];

		hv[i][0] = hv[0][RES];
		hv[i][1] = hv[1][RES + 1];
		hv[i][RES + HALO] = hv[i][HALO];
		hv[i][(RES + 2 * HALO) - 1] = hv[i][2 * HALO - 1];
	}

}

void CalculateXFluxes(const double h[RES + 2*HALO][RES + 2*HALO], const double hu[RES + 2*HALO][RES + 2*HALO],
  const double hv[RES + 2*HALO][RES + 2*HALO], double F[RES + 1][RES][3]){
  	int i;
  	int j;
  	int io;
  	int jo;
  	double h_left;
  	double hu_left;
  	double hv_left;
  	double h_right;
  	double hu_right;
  	double hv_right;
  	double h_flux_left;
  	double hu_flux_left;
  	double hv_flux_left;
  	double h_flux_right;
  	double hu_flux_right;
  	double hv_flux_right;
  	double max_wave_speed;

	for(i = 0; i < RES + 1; i++){
		for(j = 0; j < RES; j++){
        	io = i + HALO;
        	jo = j + HALO;

        // Calculate left state (to 2nd order accuracy)
        	h_left  = 0.25 * h[io][jo]  + h[io-1][jo]  - 0.25 * h[io-2][jo];
        	hu_left = 0.25 * hu[io][jo] + hu[io-1][jo] - 0.25 * hu[io-2][jo];
        	hv_left = 0.25 * hv[io][jo] + hv[io-1][jo] - 0.25 * hv[io-2][jo];

        // Calculate right state (to 2nd order accuracy)
        	h_right  = 0.25 * h[io-1][jo]  + h[io][jo]  - 0.25 * h[io+1][jo];
        	hu_right  = 0.25 * hu[io-1][jo]  + hu[io][jo]  - 0.25 * hu[io+1][jo];
        	hv_right  = 0.25 * hv[io-1][jo]  + hv[io][jo]  - 0.25 * hv[io+1][jo];

        // Left flux
        	h_flux_left = hu_left;
        	hu_flux_left = hu_left * hu_left / h_left + 0.5 * gravity * h_left * h_left;
        	hv_flux_left = hu_left * hv_left / h_left;

        	h_flux_right = hu_right;
        	hu_flux_right = hu_right * hu_right / h_right + 0.5 * gravity * h_right * h_right;
        	hv_flux_right = hu_right * hv_right / h_right;

        // Max wave speed
        	max_wave_speed = abs(0.5 * (hu_left / h_left + hu_right / h_right)) + sqrt(gravity * 0.5 * (h_left + h_right));

        // Calculate flux
        	F[i][j][0] = 0.5 * ( h_flux_left +  h_flux_right) - 0.5 * max_wave_speed * (h_right - h_left);
        	F[i][j][1] = 0.5 * (hu_flux_left + hu_flux_right) - 0.5 * max_wave_speed * (hu_right - hu_left);
        	F[i][j][2] = 0.5 * (hv_flux_left + hv_flux_right) - 0.5 * max_wave_speed * (hv_right - hv_left);

		}
	}
}

void UpdateXFluxes(double h[RES + 2*HALO][RES + 2*HALO], double hu[RES + 2*HALO][RES + 2*HALO],
  double hv[RES + 2*HALO][RES + 2*HALO], const double XFlux[RES + 1][RES][3], double dx, double dt){
  	int io;
  	int jo;
  	int i,j;

	for(i = 0; i < RES + 1; i++){
		for(j = 0; j < RES; j++){
			io = i + HALO;
        	jo = j + HALO;

        	h[io][jo]   = h[io][jo]   + dt / dx * XFlux[i][j][0];
        	h[io-1][jo] = h[io-1][jo] - dt / dx * XFlux[i][j][0];

        	hu[io][jo]   = hu[io][jo]   + dt / dx * XFlux[i][j][1];
        	hu[io-1][jo] = hu[io-1][jo] - dt / dx * XFlux[i][j][1];

        	hv[io][jo]   = hv[io][jo]   + dt / dx * XFlux[i][j][2];
        	hv[io-1][jo] = hv[io-1][jo] - dt / dx * XFlux[i][j][2];

		}
	}
}

void CalculateYFluxes(const double h[RES + 2*HALO][RES + 2*HALO], const double hu[RES + 2*HALO][RES + 2*HALO],
  const double hv[RES + 2*HALO][RES + 2*HALO], double F[RES][RES + 1][3]){
  	int i;
  	int j;
  	int io;
  	int jo;
  	double h_left;
  	double hu_left;
  	double hv_left;
  	double h_right;
  	double hu_right;
  	double hv_right;
  	double h_flux_left;
  	double hu_flux_left;
  	double hv_flux_left;
  	double h_flux_right;
  	double hu_flux_right;
  	double hv_flux_right;
  	double max_wave_speed;

	for(i = 0; i < RES; i++){
		for(j = 0; j < RES + 1; j++){
        	io = i + HALO;
        	jo = j + HALO;

        // Calculate left state (to 2nd order accuracy)
        	h_left  = 0.25 * h[io][jo]  + h[io][jo-1]  - 0.25 * h[io][jo-2];
        	hu_left = 0.25 * hu[io][jo] + hu[io][jo-1] - 0.25 * hu[io][jo-2];
        	hv_left = 0.25 * hv[io][jo] + hv[io][jo-1] - 0.25 * hv[io][jo-2];

        // Calculate right state (to 2nd order accuracy)
        	h_right  = 0.25 * h[io][jo-1]  + h[io][jo]  - 0.25 * h[io][jo+1];
        	hu_right  = 0.25 * hu[io][jo-1]  + hu[io][jo]  - 0.25 * hu[io][jo+1];
        	hv_right  = 0.25 * hv[io][jo-1]  + hv[io][jo]  - 0.25 * hv[io][jo+1];

        // Left flux
        	h_flux_left = hu_left;
        	hv_flux_left = hu_left * hu_left / h_left + 0.5 * gravity * h_left * h_left;
        	hu_flux_left = hu_left * hv_left / h_left;

        	h_flux_right = hu_right;
        	hv_flux_right = hu_right * hu_right / h_right + 0.5 * gravity * h_right * h_right;
        	hu_flux_right = hu_right * hv_right / h_right;

        // Max wave speed
        	max_wave_speed = abs(0.5 * (hu_left / h_left + hu_right / h_right)) + sqrt(gravity * 0.5 * (h_left + h_right));

        // Calculate flux
        	F[i][j][0] = 0.5 * ( h_flux_left +  h_flux_right) - 0.5 * max_wave_speed * (h_right - h_left);
        	F[i][j][1] = 0.5 * (hu_flux_left + hu_flux_right) - 0.5 * max_wave_speed * (hu_right - hu_left);
        	F[i][j][2] = 0.5 * (hv_flux_left + hv_flux_right) - 0.5 * max_wave_speed * (hv_right - hv_left);

		}
	}
}

void UpdateYFluxes(double h[RES + 2*HALO][RES + 2*HALO], double hu[RES + 2*HALO][RES + 2*HALO],
  double hv[RES + 2*HALO][RES + 2*HALO], const double YFlux[RES][RES + 1][3], double dx, double dt){
  	int io;
  	int jo;
  	int i, j;

	for(i = 0; i < RES; i++){
		for(j = 0; j < RES + 1; j++){
			io = i + HALO;
        	jo = j + HALO;

        	h[io][jo]   = h[io][jo]   + dt / dx * YFlux[i][j][0];
        	h[io][jo-1] = h[io][jo-1] - dt / dx * YFlux[i][j][0];

        	hu[io][jo]   = hu[io][jo]   + dt / dx * YFlux[i][j][1];
        	hu[io][jo-1] = hu[io][jo-1] - dt / dx * YFlux[i][j][1];

        	hv[io][jo]   = hv[io][jo]   + dt / dx * YFlux[i][j][2];
        	hv[io][jo-1] = hv[io][jo-1] - dt / dx * YFlux[i][j][2];

		}
	}
}

int offset(int x, int y, int t){
	return (t * RES * RES) + (y * RES) + x;
}

double ShallowWaterModel(int fcnnum, double t_final, int n_res, const double x_domain[2], const double y_domain[2], double* X, double* Y, double **T, double Minith[RES][RES], double Minithu[RES][RES], double Minithv[RES][RES], double **Mh, double **Mhu, double **Mhv){
  double phi0 =  0.78539816339744828;//latitude
  double cfl = 0.5;
  double omega = pow(7.292, -5);//rotation rate

  double f = 2 * omega * sin(phi0);
  int halo = 2;

  double dx = (x_domain[1] - x_domain[0]) / n_res;
  double dy = (y_domain[1] - y_domain[0]) / n_res;

  int i;
  int j;
  double XEdge[n_res + 2];
  double YEdge[n_res + 2];
  double XFlux[RES + 1][RES][3];
  double YFlux[RES][RES + 1][3];


  for(i = 0; i <= n_res + 1; i++){
    XEdge[i] = i * dx + x_domain[0];
    YEdge[i] = i * dy + y_domain[0];
  }  

  // for(i = 0; i <= n_res + 1; i++){
  //   printf("XEdge: ");
  //   printf("%f ", XEdge[i]);
  // }  
  // printf("\n");

  // for(i = 0; i <= n_res + 1; i++){
  //   printf("YEdge: ");
  //   printf("%f ", YEdge[i]);
  // }  
  // printf("\n");
  //printf("Input number to select function: 1 - Gaussian, 2 - Geostrophic, 3 - Gradient");
  //scanf("%f",&fcnnum);
  for(i = 0; i < n_res; i++){
    X[i] = 0.5 * (XEdge[i] + XEdge[i]);
    Y[i] = 0.5 * (YEdge[i] + YEdge[i]);
  }

  /*for(i = 0; i < n_res; i++){
    printf("XCentroid: ");
    printf("%f ", X[i]);
  }  
  printf("\n");

  for(i = 0; i < n_res; i++){
    printf("YCentroid: ");
    printf("%f ", Y[i]);
  } 

  printf("\n"); */

  double h[RES + 2*halo][RES + 2*halo];//height and momentums
  double hu[RES + 2*halo][RES + 2*halo];
  double hv[RES + 2*halo][RES + 2*halo];
  double h_old[RES + 2*halo][RES + 2*halo];//height and momentums
  double hu_old[RES + 2*halo][RES + 2*halo];
  double hv_old[RES + 2*halo][RES + 2*halo];
  int io;
  int jo;

  for(i = 0; i < n_res; i++){
    for(j = 0; j < n_res; j++){
      io = i + halo;
      jo = j + halo;

      
	  //printf("Input function #: 1 - Gaussian, 2 - Geostrophic, 3 - Gradient");
	  //scanf("%d", &fcnnum);
	  
	  //switch (fcnnum){
		//case 1:
			GaussianFcn(X[i], Y[i], &(h[io][jo]), &(hu[io][jo]), &(hv[io][jo]));	
		//case 2:
		//	GeostrophicFcn(X[i], Y[i], &(h[io][jo]), &(hu[io][jo]), &(hv[io][jo]), omega);
	//	case 3:
	//		GradientFcn(X[i], Y[i], &(h[io][jo]), &(hu[io][jo]), &(hv[io][jo]));
	 // }
	  

      hu[io][jo] *= h[io][jo];
      hv[io][jo] *= h[io][jo];

      //doing Minit initial states with less loops
      if((io >= halo && io < n_res + halo) && (jo >= halo && jo < n_res + halo)){
        Minith[i][j] = h[io][jo];
        Minithu[i][j] = hu[io][jo] / Minith[i][j];
        Minithv[i][j] = hv[io][jo] / Minith[i][j];                 

      }      

    }
  }

  // for(i = 0; i < RES + 2*halo; i++){
  //   for(j = 0; j < RES + 2*halo; j++){
  //     printf("h: %f ", h[i][j]);
  //     if(j == (RES + 2*halo) - 1)
  //       printf("\n");
  //   }
  // }

  // for(i = 0; i < RES + 2*halo; i++){
  //   for(j = 0; j < RES + 2*halo; j++){
  //     printf("hu: %f ", hu[i][j]);
  //     if(j == (RES + 2*halo) - 1)
  //       printf("\n");
  //   }
  // }

  //calculate initial timestep
  double max_wave_speed = 0;
  double wave_speed;
  for(i = 0; i < n_res; i++){
  	for(j = 0; j < n_res; j++){
  		wave_speed = sqrt(hu[i][j] * hu[i][j]) / h[i][j] + sqrt(gravity * h[i][j]);

  		if(wave_speed > max_wave_speed)
  			max_wave_speed = wave_speed;
  	}
  }

  double dt = cfl * dx / max_wave_speed;

  double nt = ceil(t_final / dt);
  
  printf("Time step size: %f seconds\n", dt);
  printf("Number of time steps: %f\n", nt);

  *T = (double *)malloc(nt * sizeof(double));
  *Mh = (double *)malloc(nt * RES * RES * sizeof(double)); 
  *Mhu = (double *)malloc(nt * RES * RES * sizeof(double));  
  *Mhv = (double *)malloc(nt * RES * RES * sizeof(double));

  if(*T)
  	printf("T malloc success\n");
  if(*Mh)
  	printf("Mh malloc success\n");
  if(*Mhu)
  	printf("Mhu malloc success\n");
  if(*Mhv)
  	printf("Mhv malloc success\n");
  

  double t_current = 0;
  int t;

  for(t = 0; t < (int)nt; t++){
  	//special handling for final timestep
  	if(t == nt)
  		dt = t_final - t_current;

  	ApplyBoundaryConditions(h, hu, hv);

  	/*for(i = 0; i < RES + 2*halo; i++){
    	for(j = 0; j < RES + 2*halo; j++){
      	printf("h: %f ", h[i][j]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	}*/

  	//printf("\n");

  	//store old values
  	for(i = 0; i < RES + 2*halo; i++){
    	for(j = 0; j < RES + 2*halo; j++){
  			h_old[i][j] = h[i][j];
  			hu_old[i][j] = hu[i][j];
  			hv_old[i][j] = hv[i][j];
  		}
  	}

  	//predictor step
  	CalculateXFluxes(h, hu, hv, XFlux);
  	CalculateYFluxes(h, hu, hv, YFlux);

  	//printing inverse for some reason
  	/*for(i = 0; i < RES + 1; i++){
    	for(j = 0; j < RES; j++){
      	printf("XFlux1: %f ", XFlux[i][j][0]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	}

  	for(i = 0; i < RES + 1; i++){
    	for(j = 0; j < RES; j++){
      	printf("XFlux2: %f ", XFlux[i][j][1]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	}

	for(i = 0; i < RES + 1; i++){
    	for(j = 0; j < RES; j++){
      	printf("XFlux3: %f ", XFlux[i][j][2]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	}

  	for(i = 0; i < RES; i++){
    	for(j = 0; j < RES + 1; j++){
      	printf("YFlux1: %f ", YFlux[i][j][0]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	}

  	for(i = 0; i < RES; i++){
    	for(j = 0; j < RES + 1; j++){
      	printf("YFlux2: %f ", YFlux[i][j][1]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	}

	for(i = 0; i < RES; i++){
    	for(j = 0; j < RES + 1; j++){
      	printf("YFlux3: %f ", YFlux[i][j][2]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	} */ 	

  	//update half a timestep
  	UpdateXFluxes(h, hu, hv, XFlux, dx, 0.5 * dt);
  	UpdateYFluxes(h, hu, hv, YFlux, dx, 0.5 * dt);

  	/*for(i = 0; i < RES + 2*halo; i++){
    	for(j = 0; j < RES + 2*halo; j++){
      	printf("h: %f ", h[i][j]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	}

  	for(i = 0; i < RES + 2*halo; i++){
    	for(j = 0; j < RES + 2*halo; j++){
      	printf("hu: %f ", hu[i][j]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	}

  	for(i = 0; i < RES + 2*halo; i++){
    	for(j = 0; j < RES + 2*halo; j++){
      	printf("hv: %f ", hv[i][j]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	}*/

  	// Apply source terms
    for(i = 0; i < RES; i++){
    	for(j = 0; j < RES; j++){
        	io = i + halo;
        	jo = j + halo;

        	hu[io][jo] = hu[io][jo] + 0.5 * dt * f * hv_old[io][jo];
        	hv[io][jo] = hv[io][jo] - 0.5 * dt * f * hu_old[io][jo];
    	}
    }

    ApplyBoundaryConditions(h, hu, hv);

  	CalculateXFluxes(h, hu, hv, XFlux);
  	CalculateYFluxes(h, hu, hv, YFlux);

  	for(i = 0; i < RES + 2*halo; i++){
    	for(j = 0; j < RES + 2*halo; j++){
  			h[i][j] = h_old[i][j];
  			hu[i][j] = hu_old[i][j];
  			hv[i][j] = hv_old[i][j];
  		}
  	}

  	UpdateXFluxes(h, hu, hv, XFlux, dx, dt);
  	UpdateYFluxes(h, hu, hv, YFlux, dx, dt);

  	t_current = t_current + dt;

  	(*T)[t] = t_current;

  	/*printf("t_current: %f\n", t_current);
  	printf("T: %f\n", (*T)[t]);

  	printf("\n");

  	for(i = 0; i < RES + 2*halo; i++){
    	for(j = 0; j < RES + 2*halo; j++){
      	printf("h: %f ", h[i][j]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	}

  	printf("\n");

  	for(i = 0; i < RES + 2*halo; i++){
    	for(j = 0; j < RES + 2*halo; j++){
      	printf("hu: %f ", hu[i][j]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	}

  	printf("\n");

  	for(i = 0; i < RES + 2*halo; i++){
    	for(j = 0; j < RES + 2*halo; j++){
      	printf("hv: %f ", hv[i][j]);
      	if(j == (RES + 2*halo) - 1)
        	printf("\n");
    	}
  	}

  	  printf("\n");
*/
  	int off = 0;
  	for(i = 0; i < RES; i++){
  		for(j = 0; j < RES; j++){
  			off = offset(i, j, t);
  			(*Mh)[off] = h[HALO + i][HALO + j];
  			//printf("h: %f\n", h[HALO + i][HALO + j]);
  			(*Mhu)[off] = hu[HALO + i][HALO + j] / h[HALO + i][HALO + j];
  			//printf("hu: %f\n", hu[HALO + i][HALO + j]);
  			(*Mhv)[off] = hv[HALO + i][HALO + j] / h[HALO + i][HALO + j];
  			//printf("hv: %f\n", hv[HALO + i][HALO + j]);
  		}
  	}
   
  }

  return nt;
}

int main(){
  int fcnnum;
  int n_res = RES;
  double t_final = 7200;
  double x_domain[2] = {-1000000, 1000000};
  double y_domain[2] = {-1000000, 1000000};
  double X[RES];
  double Y[RES];
  double *T;
  double *Mh;
  double *Mhu;
  double *Mhv;
  double Minith[n_res][n_res];
  double Minithu[n_res][n_res];
  double Minithv[n_res][n_res];
  int i;
  int j;
  double nt;
  
  
  
  printf("Executing ShallowWaterModel");  
  nt = ShallowWaterModel(fcnnum, t_final, n_res, x_domain, y_domain, X, Y, &T, Minith, Minithu, Minithv, &Mh, &Mhu, &Mhv);

  /*for(i = 0; i < n_res; i++){
    printf("X: ");
    printf("%f ", X[i]);
  }  
  printf("\n");

  for(i = 0; i < n_res; i++){
    printf("Y: ");
    printf("%f ", Y[i]);
  } 

  printf("\n");*/

  
  for(i = 0; i < n_res; i++){
    for(j = 0; j < n_res; j++){
      printf("Minith: ");
      printf("%f\t", Minith[i][j]);
      if(j == n_res - 1)
        printf("\n");
    }
  } 

  printf("\n");

  for(i = 0; i < n_res; i++){
    for(j = 0; j < n_res; j++){
      printf("Minithu: ");
      printf("%f\t", Minithu[i][j]);
        if(j == n_res - 1)
        printf("\n");
    }
  } 

  printf("\n");

  for(i = 0; i < n_res; i++){
    for(j = 0; j < n_res; j++){
      printf("Minithv: ");
      printf("%f\t", Minithv[i][j]);
        if(j == n_res - 1)
        printf("\n");
    }
  } 

  printf("\n");

  printf("nt: %f\t", nt);

  for(i = 0; i < n_res; i++){
    for(j = 0; j < n_res; j++){
      int off = offset(i, j, nt - 1);
      printf("Mh: ");
      printf("%f\t", Mh[off]);
        if(j == n_res - 1)
        printf("\n");
    }
  }

  for(i = 0; i < n_res; i++){
    for(j = 0; j < n_res; j++){
      int off = offset(i, j, nt - 1);
      printf("Mhu: ");
      printf("%f\t", Mhu[off]);
        if(j == n_res - 1)
        printf("\n");
    }
  }

  for(i = 0; i < n_res; i++){
    for(j = 0; j < n_res; j++){
      int off = offset(i, j, nt - 1);
      printf("Mhv: ");
      printf("%f\t", Mhv[off]);
        if(j == n_res - 1)
        printf("\n");
    }
  }

  /*
  int t;
  for(t = 0; t < (int)nt; t++){
  	printf("T: %f\n", T[t]);

  	printf("\n");
  }*/   

  free(T);
  free(Mh);
  free(Mhv);
  free(Mhu);

  return 0;
}
