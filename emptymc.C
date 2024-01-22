
//****************************************************************************************
// CODE for MC exercises (NVT & NPT Hard Spheres)
// Debye Summer School 2014
//******************************************************************************************

//******************************************************************************************
// Load libraries
//******************************************************************************************

#include <stdio.h>
// Standard Input/output

#include <math.h>
// Mathematics library

#include <stdlib.h>
// Standard Library

#include <time.h>
// Use of the clock


//******************************************************************************************
// Random Number Generator
//******************************************************************************************

#include "mt19937ar.C"
// Random number generator

//******************************************************************************************
// FIXED PARAMETERS
//******************************************************************************************

#define N 500
// Number of particles

#define SimType 1
// Defines whether the simulation is NVT or NPT:
// 1 -> NVT
// 2 -> NPT

const double pressure = 5.0;
//Pressure in case of NPT simulations

const int cycles = 200000;
//The number of MC cycles per simulation

const int initcycles = 150000;
//The number of equilibration cycles per simulation

const double delta    = 0.05;
// Step size for particle moves

const double deltav   = 0.1;
//Volume step size

const double D = 0.010;
//Parameter to keep track of $\delta$ in the stretch

const double startingdensity  = 1.25;
//Starting density

double system_energy = 0.0;
// Global variable to keep track of system energy

// Lennard-Jones Parameters: 

const double sigma = 1.0;
//simulation reduced lennard-jones distance parameter

const double epsilon = 1.0; 
//simlation reduced lennard-jones energy parameter

const double kT = 1.0; 
//simulation variable corresponding to reduced lennard-jones parameters kT 

const double rcutoff = pow(2.0, 1.0/6.0); // Corrected cutoff distance for WCA


//******************************************************************************************
// GLOBAL VARIABLES
//******************************************************************************************

double boxx, boxy, boxz;
// The variables for the simulation box-size

double r[N][3];
//Variable for the coordinates

double tensor_product[3][3];
//Variable for the tensor product within the expression of virial pressure 
//******************************************************************************************
// FUNCTION DECLARATIONS
//******************************************************************************************
int move(int part);
int volumemove();

double wca_potential(double rij);
double potential_energy(int part);
double total_energy();

void setdensity();
void measure();
void readcoords();
void writecoords();
void initrandomseed();
void print_tensor_product();
int overlap(int part);

double genrand();                       //Defined in mt19937ar.c.  Generates a random number in [0,1).
void init_genrand(unsigned long s);     //Defined in mt19937ar.c.  Initializes the random number generator with seed s.

//******************************************************************************************
// MAIN CODE
//******************************************************************************************
int main()
{
  // ------------- Initialization ---------------- 
  double rcutoff2 = rcutoff * rcutoff; 
  // Square of cutoff distance for comparison

  initrandomseed();
  //Initialize random number generator

  readcoords();
  //Read coordinates from init.dat

  setdensity();
  //Set density to the correct value

  //writecoords();
  //Write starting coordinates

  // ---------------  MC local variables -------------------
  int cyc; // counter for number of cycles
  int step; // counter for the number of steps (particle or volume moves) per cycle
  int accepted; // dummy variable to check acceptance of moves
  int tvmoves = 0; //Counter for number of trial volume moves
  int tpmoves = 0; //Counter for number of trial particle moves
  int vmoves = 0;   //Counter for number of accepted volume moves
  int pmoves = 0;   //Counter for number of accepted particle moves


  // Open a text file for debugging
  // FILE *outputfile = fopen("Mcdebugging//.txt", "w");

  // Generating the filename using the variables
  //char filename1[50];
  //sprintf(filename1, "MCResults/file_%d_%.4f.txt", N, boxx);  

  // Open a text file for debugging
  //FILE *outputfile1 = fopen(filename1, "w");
  double potential_energy_value = total_energy();


  // ------------- MC simulation ----------------
  system_energy = potential_energy_value;
  printf("Total Energy before Simulation: %1f \n", potential_energy_value);

  
  // MC CYCLES (cyc)
  for (cyc = 0; cyc < cycles; cyc++)
  {
    // MC STEPS (VOLUME & PARTICLE MOVES)  (step)
    for ( step = 0; step < N; step++)
    {
      // Declare a variable for the index of the particle to move, if part = N, then that is a volume move.
      int part;

      //printf ("Step = %d\n", cyc);
      // Determine if the move will be a particle move or a volume move, and if a particle move, which particle is attempted to move.
      if (SimType == 1)  //NVT : Each cycle consists of N particle moves
      {
        part = (int) (genrand_real2() * (N));       //Generate a random integer in [0,...,N-1]
      }
      else if (SimType == 2) //NPT : Each cycle contains on average N particle moves and 1 volume move once every (N + 1) steps
      {
        part = (int) (genrand() * (N + 1 ));	//Generate a random integer in [0,...,N]
      }
      else
      {
         printf( "SimType must be either 1 (NVT) or 2 (NPT).\n");
         exit(1);
      }

      // VOLUME MOVE, tries to change the volume of the system
      if (part == N)
      {
        tvmoves ++; // Increase the number of trial volume moves
        accepted = volumemove();		//Change the volume
        if ( accepted ) vmoves ++; // If move was accepted, increase the number of accepted volume moves
      }
      // PARTICLE MOVE, tries to move particle "part"
      else
      {
        tpmoves ++; // Increase the number of trial particle moves
        accepted = move(part);			//Move the chosen particle
        if (accepted) pmoves++; // If a trial particle move was accepted, increase the number of accepted particle moves.
      }
    }

    

    // When cyc > initcycles, we can measure as the system has finished equilibrating, defined by initcycles
    if (cyc > initcycles) //&& cyc%20==0)
    {
      
      measure();				    //Measure density
      //writecoords();				//Write coordinates
      //potential_energy_value = total_energy();
      // Write potential energy and cycle number to the output file
      //fprintf(outputfile1, "%d %lf\n", cyc, potential_energy_value/N);
      
    }
  }

  //writecoords();					//Write last configuration
  printf("Test TPmoves %d Pmoves: %d\n", tpmoves, pmoves);
  printf ("Fraction of accepted particle moves = %lf\n",(double)pmoves/(double)tpmoves);
  printf("Total Energy after Simulation: %1f \n", total_energy());
  if (SimType==2) printf ("Fraction of accepted volume moves   = %lf\n",(double)vmoves/(double)tvmoves);

  return 0;
}

//***********************************************************************************
//***********************************************************************************
//
// MOVE : Tries to move a particle called "part"
// ---  To be filled in by student !!!
//
// Useful information:
// --- Note that a random number between 0 and 1 can be generated as follows:
//      x = genrand();
//     Note that x should be a "double".
// --- Note that the array containing the particle positions is called r, such that
//     r[25][2]  is the z component of the 26th particle.  Remember that in C, arrays start at 0, so that
//     r[0][0] is the x coordinate of the 1st particle in the array.
// --- You can use a function called "overlap" to determine if the new position of a particle overlaps
//     with any other particles in the simulations.  This function works as follows:
//     x = overlap(ipart)
//     x is 0 if there are no overlaps between particle "ipart" and any other particle in the system, and
//     x is 1 if there are overlaps
//     Note that x should be an integer
// --- Move should return 1 if the move was succesful, and 0 if the move was rejected.
//
//***********************************************************************************
//***********************************************************************************

int move (int part)
{
  
  // Generate random displacements in each direction
  double dx = (genrand_real2() - 0.5) * 2.0 * delta;
  double dy = (genrand_real2() - 0.5) * 2.0 * delta;
  double dz = (genrand_real2() - 0.5) * 2.0 * delta;

  // Backup the current coordinates of the particle
  double old_x = r[part][0];
  double old_y = r[part][1];
  double old_z = r[part][2];

  
  double E1 = potential_energy(part); // Current total energy  (Use potential energy)
  
  // Apply the proposed move
  r[part][0] += dx;
  r[part][1] += dy;
  r[part][2] += dz;

  // Apply periodic boundary conditions to bring particles back into the box
  while (r[part][0] >= boxx) r[part][0] -= boxx;
  while (r[part][0] < 0) r[part][0] += boxx;
  while (r[part][1] >= boxy) r[part][1] -= boxy;
  while (r[part][1] < 0) r[part][1] += boxy;
  while (r[part][2] >= boxz) r[part][2] -= boxz;
  while (r[part][2] < 0) r[part][2] += boxz;

  // Calculate the change in energy (dE) due to the move 
  double E2 = potential_energy(part);     // New total energy after move
  double dE = E2 - E1; //calculating the change in energy 
  //printf("%lf \n", dE);

  if (dE <= 0 || genrand_real2() < exp(-dE / (kT))) {
    system_energy += dE;
    return 1; // Return 1 to indicate the move was accepted

  }else {

    // Reject the move, restore particle coordinates
    r[part][0] = old_x;
    r[part][1] = old_y;
    r[part][2] = old_z;
    
    return 0; // Return 0 to indicate the move was rejected
  }


}

//***********************************************************************************
//***********************************************************************************
//
//VOLUME MOVE : Tries to change the volume during the simulation
// ---  To be filled in by student !!!
//
// Useful information:
// --- Note that a random number between 0 and 1 can be generated as follows:
//      x = genrand();
//     Note that x should be a "double".
// --- The simulation box is given by the three axis lengths: boxx, boxy, and boxz.
// --- Note that the array containing the particle positions is called r, such that
//     r[25][2]  is the z component of the 26th particle.  Remember that in C, arrays start at 0, so that
//     r[0][0] is the x coordinate of the 1st particle in the array.
// --- You can use a function called "overlap" to determine if the new position of a particle overlaps
//     with any other particles in the simulations.  This function works as follows:
//     x = overlap(ipart)
//     x is 0 if there are no overlaps between particle "ipart" and any other particle in the system, and
//     x is 1 if there are overlaps
//     Note that x should be an integer
// --- Volumemove should return 1 if the move was succesful, and 0 if the move was rejected.
// --- You can take a cube root using the function cbrt(x), a natural logarithm using log(x),
//     and an exponential using exp(x).
//
//***********************************************************************************
//***********************************************************************************
int volumemove()
{
  return 0;
}

double compute_force(double rij) {

    double r2, r6, r12;
    r2 = sigma * sigma / (rij * rij); // (sigma/rij)^2
    r6 = r2 * r2 * r2;  // (sigma/rij)^6
    r12 = r6 * r6;      // (sigma/rij)^12
    
    // Compute the force magnitude
    return 48.0 * epsilon * (r12 - 0.5*r6 ) / (rij);
}

// Function to calculate WCA potential energy given the distance rij
//use periodic boundary condition
double wca_potential(double rij) {

    if (rij > rcutoff) {
        return 0.0; // Potential is zero beyond the cutoff distance (is this even necessary)
    } else {
        
        double r2 = sigma * sigma / (rij * rij); // (sigma/rij)^2
        double r6 = r2 * r2 * r2;                // (sigma/rij)^6
        double r12 = r6 * r6;                    // (sigma/rij)^12
        return 4.0 * epsilon * (r12 - r6) + epsilon; // WCA potentiall
    }
}


double potential_energy(int part) {
  
  double U_part = 0.0;
  double dx, dy, dz, r1; // Note: using r2 for squared distance
  int i;

  //FILE *f = fopen("WCA.txt", "a"); // Open file in append mode

  for (i = 0; i < N; i++) { // Loop over all particles
    if (i != part) { // Do not check overlap with self
      
      dx = r[part][0] - r[i][0]; // Calculate distance vector
      dy = r[part][1] - r[i][1];
      dz = r[part][2] - r[i][2];

      // Periodic boundaries: Use nearest image convention
      while (dx >  boxx*0.5) dx -= boxx;
      while (dx < -boxx*0.5) dx += boxx;
      while (dy >  boxy*0.5) dy -= boxy; // Using boxy for y dimension
      while (dy < -boxy*0.5) dy += boxy; // Using boxy for y dimension
      while (dz >  boxz*0.5) dz -= boxz; // Using boxz for z dimension
      while (dz < -boxz*0.5) dz += boxz; // Using boxz for z dimension

      r1 = sqrt(dx*dx + dy*dy + dz*dz); // Squared distance 

      if (r1 <= rcutoff) { // Only consider interactions within squared rcutoff
        U_part += wca_potential(r1);
      }
    }
  }
  //fclose(f);
  return U_part;
}


// Function to calculate the total energy of the system
double total_energy() {
    double total_energy = 0.0;

    for (int i = 0; i < N; i++) {
        total_energy += potential_energy(i);
    }

    return total_energy * 0.5;
}


//***********************************************************************************
//***********************************************************************************
//
// MEASURE : Keeps track of the density during the simulation
// ---  To be filled in by student !!!
//
// Useful information:
// --- The simulation box is given by the three axis lengths: boxx, boxy, and boxz.
// --- Note that the array containing the particle positions is called r, such that
//     r[25][2]  is the z component of the 26th particle.  Remember that in C, arrays start at 0, so that
//     r[0][0] is the x coordinate of the 1st particle in the array.
// --- If you want to write to a file, you have to first declare a file pointer:
//         FILE* outpufile;
//     Then, you can open and empty it using:
//         outputfile = fopen ("filename.dat", "w");
//     Or, you can open and add to the end using:
//         outputfile = fopen ("filename.dat", "a");
//     To write to a file, use fprintf:
//         fprintf (outputfile, "%lf\n", varname)
//     This prints a real number (taken from the variable 'varname') to a line, and ends the line.
//     Note that integers are indicated with %d (instead of %lf) in this format.
//***********************************************************************************
//***********************************************************************************

void measure() {
  
  double dx, dy, dz, r1, fij;

  // Initialize the stress tensor to zero
  for (int a = 0; a < 3; a++)
      for (int b = 0; b < 3; b++)
          tensor_product[a][b] = 0.0;

  // Iterate over all particle pairs (i, j)
  for (int i = 0; i < N; i++) {
      for (int j = 0; j < i; j++) {
        dx = r[i][0] - r[j][0];
        dy = r[i][1] - r[j][1];
        dz = r[i][2] - r[j][2];
        
        while (dx >  boxx*0.5) dx -= boxx;		//Periodic boundaries: Use nearest image convention
        while (dx < -boxx*0.5) dx += boxx;
        while (dy >  boxy*0.5) dy -= boxy;
        while (dy < -boxy*0.5) dy += boxy;
        while (dz >  boxz*0.5) dz -= boxz;
        while (dz < -boxz*0.5) dz += boxz;

        r1 = sqrt(dx*dx + dy*dy + dz*dz); // Squared distance 

        if (r1 < rcutoff) {
          // Compute scalar force contribution
          fij = compute_force(r1); // Call the correct force function here
          tensor_product[0][0] += dx * (dx/r1) * fij;
          tensor_product[0][1] += dx * (dy/r1) * fij;
          tensor_product[0][2] += dx * (dz/r1) * fij;
          tensor_product[1][0] += dy * (dx/r1) * fij;
          tensor_product[1][1] += dy * (dy/r1) * fij;
          tensor_product[1][2] += dy * (dz/r1) * fij;
          tensor_product[2][0] += dz * (dx/r1) * fij;
          tensor_product[2][1] += dz * (dy/r1) * fij;
          tensor_product[2][2] += dz * (dz/r1) * fij;
          

      }
      

      
     }
  }
  
  print_tensor_product();
}

void print_tensor_product() {

  // Generating the filename for eq of state using the variables
  
  char filename3[50];
  sprintf(filename3, "MCResults/Betaepsilon1/pressureoutput_stretch_%0.4f.txt", D);
  //sprintf(filename3, "MCResults/Betaepsilon40/pressureoutput_%0.4f.txt", startingdensity);
  FILE *outputfile2 = fopen(filename3, "a");
  int i, j;
  //printf("Stress Tensor:\n");
  for (i = 0; i < 3; i++) {
      //printf("[");
      for (j = 0; j < 3; j++) {
          //printf("%lf\t ,", tensor_product[i][j]); // Print each element with a tab space
          fprintf(outputfile2, "%lf \t", tensor_product[i][j]);
      }
      fprintf(outputfile2, "\n");
      //printf("]\n"); // Move to the next line after printing each row
  }
  fclose(outputfile2);
}

//***********************************************************************************
// SETDENSITY: Set density to the correct starting density
//***********************************************************************************
void setdensity()
{
  int i;
  double volumenew = N / startingdensity;		//Calculate scaling factor
  double volumeold = boxx * boxy * boxz;
  double scalefactor = pow(volumenew / volumeold, 1.0/3.0);


  boxx *= scalefactor;				//Update box size
  boxy *= scalefactor;
  boxz *= scalefactor;
  for ( i = 0; i < N; i++)			//Update particle positions
  {
    r[i][0] *= scalefactor;
    r[i][1] *= scalefactor;
    r[i][2] *= scalefactor;
  }

  for (i=0; i<N; i++)
  {
    if (overlap(i)==1)
    {
      printf("Overlaps in beginning configuration.  Try lowering the starting density.\n");
      exit(1);
    }
  }
  printf("boxx = %1f \n", boxx);
  printf("boxy = %1f \n", boxy);
  printf("boxz = %1f \n", boxz);
  printf ("Scaled box: (%.2lf, %.2lf, %.2lf), Npart: %d\n", boxx, boxy, boxz, N);

  char filename2[50];
  sprintf(filename2, "MCResults/Betaepsilon1/pressureoutput_stretch_%0.4f.txt", D);
  FILE *outputfile2 = fopen(filename2, "a");
  fprintf(outputfile2, "%1f %1f %1f \n", boxx, boxy, boxz);
  fclose(outputfile2);

}


//***********************************************************************************
// WRITECOORDS: Write a snapshot to a file, following the LAMMPS trajectory format.
//              The name of the outputfile is "output.lammpstrj".
//***********************************************************************************
void writecoords()
{
  static int num = 0;         //Remember frame number
  int i;                      //Counter to loop over the particles

  char filename[100];
  sprintf (filename, "MCcoords/coords_%04d.dat", num);
  FILE* output;               //Create file pointer.
  output = fopen (filename, "w"); //Open file: create empty file

  fprintf (output, "%d\n", N);
  fprintf (output, "0.0  %lf\n", boxx);
  fprintf (output, "0.0  %lf\n", boxy);
  fprintf (output, "0.0  %lf\n", boxz);

  for (i = 0 ; i < N; i++)
  {
    fprintf (output, "%lf  %lf  %lf 1.0\n", r[i][0], r[i][1], r[i][2]);
  }

  fclose(output);

  num ++;
}


//***********************************************************************************
// READCOORDS:  Reads in the initial configuration, called "init.dat".
//***********************************************************************************
void readcoords()
{
  int i, npart;
  char filename2[50];
  sprintf(filename2, "MCinputlattice/updatedrun/xyz_stretched_%0.3f.dat", D);
  //sprintf(filename2, "MCinputlattice/updatedrun/.dat");
  FILE *initfile = fopen(filename2, "r");
  //FILE* initfile = fopen ("MCinputlattice/xyz6.dat", "r");				//Open the file

  int test = fscanf (initfile, "%d\n", &npart);					//Read number of particles
  if (test != 1) { printf ("Problem reading number of particles!\n"); exit(3);}
  if (npart != N)
  {
    
    printf ("Error! Number of particles in file does not match number at the top of the code!\n");
    exit(3);
  }

  test = fscanf (initfile, "%lf  %lf  %lf\n", &boxx, &boxy, &boxz);		//Read boxsize
  if (test != 3) { printf ("Problem reading box!\n"); exit(3);}
  
  printf ("Box: (%.2lf, %.2lf, %.2lf), Npart: %d\n", boxx, boxy, boxz, npart);

  double dummy;

  for ( i = 0; i < N; i++)
  {
    test = fscanf (initfile, "%lf %lf %lf %lf\n", &(r[i][0]), &(r[i][1]), &(r[i][2]), &dummy);
    if (test < 3) { printf ("Problem reading particles!\n"); exit(3);}
  }
  fclose ( initfile);

}



//***********************************************************************************
// OVERLAP: checks a particle (with index part) for overlaps with other particles.
//          This function returns 0 if there are no overlaps, and 1 if there is at least one.
//***********************************************************************************
int overlap(int part)
{
  double dx, dy, dz, r2;
  int i;
  int overl = 0;

  for (i = 0 ; i < N; i++)              //Loop over all particles
  {
    if ( i != part)                     //Do not check overlap with self
    {
      dx = r[part][0] - r[i][0];        //Calculate distance vector
      dy = r[part][1] - r[i][1];
      dz = r[part][2] - r[i][2];
      while (dx >  boxx*0.5) dx -= boxx;		//Periodic boundaries: Use nearest image convention
      while (dx < -boxx*0.5) dx += boxx;
      while (dy >  boxx*0.5) dy -= boxy;
      while (dy < -boxx*0.5) dy += boxy;
      while (dz >  boxx*0.5) dz -= boxz;
      while (dz < -boxx*0.5) dz += boxz;
      r2 = dx*dx+dy*dy+dz*dz;			//distance squared
      if (r2 < 1.0)
      {
        overl = 1;
        break;					//Stop checking
      }
    }
  }
  return overl;
}

//***********************************************************************************
// INITRANDOMSEED:  Initializes the random number generator.
//                    Note that the seed is taken from the current time.
//***********************************************************************************
void initrandomseed()
{
  unsigned long seed = time(NULL);
  init_genrand(seed);
}

