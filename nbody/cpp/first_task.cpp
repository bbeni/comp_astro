#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>

/*

Hermite 4th order

neeed a-dot

predictor step

xp = x0 + v0 dt + 1/2 a0 dt^2 + 1/6 a0. dt^3
xp. = x0. + a0 dt + 1/2 a0. dt^2

evaluator step
vector quantities: rij, vij, aij
xp. und xp -> rij und vij

ai = sum(mj rij / xij^3)
ai. = sum(mj vij / xij^3 - 3 alphaij aij)

corrector step (alles vector)
x1. = x0. + 1/2(a0+a1) dt + 1/12(a0. - a1.)dt^2
x1 = x0 + 1/2(x0.+x1.) dt + 1/12(a0 - a1)dt^2

trick is: first compute x1. then use it for evaluation corrected x1

best performing criterion for timestep dti
dti = (eta (|ai||ai..| + |ai.|^2) / (|ai.||ai...| + |ai..|))^1/2


plummer softening S(rij, epsilon) = -1/sqrt(rij^2 + epsilon^2)

*/


using namespace std;


#define NSTEPS 2000
#define eta 0.1 // accuracy
//#define SOFTENING 0.0004

typedef struct Particle
{
  double m;
	double x, y, z;
	double vx, vy, vz;
  double ax, ay, az;
  double jx, jy, jz; // jerk - a.

  // perdicted
  double xp, yp, zp;
  double vxp, vyp, vzp;

  // next - from evaluator step
  double axn, ayn, azn;
  double jxn, jyn, jzn;

  // calculated scalars
  double r;
  //double ekin;
  //double epot;

	double softening;
	double potential;

} Particle;

Particle new_particle()
{
  Particle p = {}; // all member variables to 0
  return p;
}

double cmpdt(std::vector<Particle> particles)
{
  double DT_min = 10.0;
  for(Particle& p : particles)
  {
    if (p.softening == 0) continue;
    double dt = eta * sqrt(p.softening/sqrt(p.ax*p.ax + p.ay*p.ay + p.az*p.az));
    if (dt < DT_min) DT_min = dt;
  }

  return DT_min;
}

void snapshot_to_csv(std::vector<Particle> particles, double DT, string filename)
{
  
  ofstream myfile(filename);
  if (myfile.is_open())
  {
    // header
    myfile << "m,x,y,z,vx,vy,vz,ax,ay,az,r,softening,potential,dt" << endl;
    for(auto const& p : particles)
    {
      myfile << p.m << ',' << p.x << ',' << p.y << ',' << p.z << ',' << p.vx << ',' << p.vy << ',' << p.vz << ',';
      myfile << p.ax << ',' << p.ay << ',' << p.az << ',' << p.r << ',' << p.softening << ',' << p.potential;
      myfile << DT;
      myfile << endl;

    }
    myfile.close();
  } 
  else cout << "Unable to open file" << endl;

}


vector<Particle> read_particle_data(string filename)
{
  std::vector<Particle> particles;

  string line;
  ifstream myfile (filename);
  if (myfile.is_open())
  {
    string header;
    getline(myfile, header);
    std::string header1 = header.substr(0, header.find(" "));
    
    int n = stoi(header1);
    cout << n << endl;

    // add n empty particles
    for (int i=0; i<n; i++)
    {
      Particle p = new_particle();
      particles.push_back(p);
    }

    std::vector<Particle>::size_type i;
    // masses
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      double m = stof(line);
      particles[i].m = m;
      i++;
    }

    // x
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      double x = stof(line);
      particles[i].x = x;
      i++;
    }

    // y
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      double x = stof(line);
      particles[i].y = x;
      i++;
    }
    
    // z
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      double x = stof(line);
      particles[i].z = x;
      i++;
    } 

    // vx
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      double x = stof(line);
      particles[i].vx = x;
      i++;
    }

    // vy
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      double x = stof(line);
      particles[i].vy = x;
      i++;
    }
    
    // vz
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      double x = stof(line);
      particles[i].vz = x;
      i++;
    }

    // softening
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      double x = stof(line);
      particles[i].softening = x;
#ifdef SOFTENING
      particles[i].softening = SOFTENING;
#endif
      i++;
    }
    // potential
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      double x = stof(line);
      particles[i].potential = x;
      i++;
    }
    cout << particles[10].m << endl;
    cout << particles[11].m << endl;


    myfile.close();
  }

  else cout << "Unable to open file";
  return particles;
}

void calc_direct_force(std::vector<Particle>& particles)
{
  double rx, ry, rz;
  double under;

  for(auto base=particles.begin(); base != particles.end(); ++base)
  {
    double ax=0, ay=0, az=0;
    for(auto it=particles.begin(); it != particles.end(); ++it)
    {
      if( base == it ) continue;
      rx = base->x - it->x;
      ry = base->y - it->y;
      rz = base->z - it->z;

      under = pow(rx*rx + ry*ry + rz*rz + base->softening*base->softening, 3.0/2.0);

      ax += it->m * rx / under;
      ay += it->m * ry / under;
      az += it->m * rz / under;
    }

    //int particle_nr = base - p.begin();
    //cout << "particle nr " << particle_nr << " " << ax << " " << ay << " " << az << endl;
    //cout << base - p.begin() << endl; 

    base->ax = ax;
    base->ay = ay;
    base->az = az;

    base->r = sqrt(pow(base->x, 2) + pow(base->y, 2) + pow(base->z, 2));

  }

}

void hermite_predict(Particle& p, double DT)
{
	// position
	p.xp = p.x + p.vx * DT + 0.5 * p.ax * DT*DT + 1.0/6.0 * p.jx * DT*DT*DT;
	p.yp = p.y + p.vy * DT + 0.5 * p.ay * DT*DT + 1.0/6.0 * p.jy * DT*DT*DT;
	p.zp = p.z + p.vz * DT + 0.5 * p.az * DT*DT + 1.0/6.0 * p.jz * DT*DT*DT;

	// velocity
	p.vxp = p.vx + p.ax * DT + 0.5 * p.jx * DT*DT;
	p.vyp = p.vy + p.ay * DT + 0.5 * p.jy * DT*DT;
	p.vzp = p.vz + p.az * DT + 0.5 * p.jz * DT*DT;
}

void hermite_correct(Particle& p, double DT)
{
	double vx1 = p.vx + 1.0/2.0*(p.ax+p.axn) * DT + 1.0/12.0*(p.jx-p.jxn) * DT*DT;
	double vy1 = p.vy + 1.0/2.0*(p.ay+p.ayn) * DT + 1.0/12.0*(p.jy-p.jyn) * DT*DT;
	double vz1 = p.vz + 1.0/2.0*(p.az+p.azn) * DT + 1.0/12.0*(p.jz-p.jzn) * DT*DT;

	// new quantities are set
	p.x = p.x + 1.0/2.0*(p.vx+vx1) * DT + 1.0/12.0*(p.ax-p.axn) * DT*DT;
	p.y = p.y + 1.0/2.0*(p.vy+vy1) * DT + 1.0/12.0*(p.ay-p.ayn) * DT*DT;
	p.z = p.z + 1.0/2.0*(p.vz+vz1) * DT + 1.0/12.0*(p.az-p.azn) * DT*DT;

	p.vx = vx1;
	p.vy = vy1;
	p.vz = vz1;

	p.ax = p.axn;
	p.ay = p.ayn;
	p.az = p.azn;
	
	p.jx = p.jxn;
	p.jy = p.jyn;
	p.jz = p.jzn;
}

void hermite_evaluate(std::vector<Particle>& particles, double DT)
{
	for(Particle& base : particles)
  {
    for(Particle& it : particles)
    {
      if( &base == &it ) continue;
      // use the predicted quantities xp, ...
      double rx = base.xp - it.xp;
      double ry = base.yp - it.yp;
      double rz = base.zp - it.zp;
      double vx = base.vxp - it.vxp;
      double vy = base.vyp - it.vyp;
      double vz = base.vzp - it.vzp;

      double rsquared = rx*rx + ry*ry + rz*rz;

      double alphaij = (rx*vx + ry*vy + rz*vz);

      double sq = rsquared + base.softening*base.softening;

      double under1 = pow(sq, 3.0/2.0);
      double under2 = pow(sq, 5.0/2.0);

      // G = 1
      base.axn -= it.m * rx / under1;
      base.ayn -= it.m * ry / under1;
      base.azn -= it.m * rz / under1;
      base.jxn -= it.m * (vx / under1 - 3 * alphaij * rx / under2);
			base.jyn -= it.m * (vy / under1 - 3 * alphaij * ry / under2);
			base.jzn -= it.m * (vz / under1 - 3 * alphaij * rz / under2);
    }
	}
}

void integrate(std::vector<Particle>& particles, double DT)
{
	for(Particle& p : particles)
	{
		hermite_predict(p, DT);
	}

	hermite_evaluate(particles, DT);

	for(Particle& p : particles)
	{
		hermite_correct(p, DT);

    // calculate scalars
		p.r = sqrt(pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2));

  }

}


int main(int argc, const char* argv[])
{
	std::string filename = "data_n2.dat";

	std::vector<Particle> particles = read_particle_data(filename);

  calc_direct_force(particles);

  for(int i = 0; i < NSTEPS; i++)
  {
  	cout << "step " << i << endl;

    double DT = cmpdt(particles);

  	integrate(particles, DT);
    
    cout << "DT is " << DT << endl;

    std::string fname = "snapshot" + std::string(6 - std::to_string(i).length(), '0') + std::to_string(i) + ".csv";
    cout << "saving snapshot " << fname << " ..." << endl;

  	snapshot_to_csv(particles, DT, fname);
  }

}