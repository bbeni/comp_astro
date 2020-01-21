#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

typedef struct Particle
{
  float m;
	float x, y, z;
	float vx, vy, vz;

  // acceleration
  float ax, ay, az;

  float r;

	float softening;
	float potential;



} Particle;

void particles_to_csv(std::vector<Particle> particles, string filename)
{
  
  ofstream myfile(filename);
  if (myfile.is_open())
  {
    // header
    myfile << "m,x,y,z,vx,vy,vz,ax,ay,az,r,softening,potential" << endl;
    for(auto const& p : particles)
    {
      myfile << p.m << ',' << p.x << ',' << p.y << ',' << p.z << ',' << p.vx << ',' << p.vy << ',' << p.vz << ',';
      myfile << p.ax << ',' << p.ay << ',' << p.az << ',' << p.r << ',' << p.softening << ',' << p.potential;
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
    //header.erase(0, pos + 1);
    int n = stoi(header1);
    cout << n << endl;

    // add n empty particles
    for (int i=0; i<n; i++)
    {
      Particle p;
      particles.push_back(p);
    }

    std::vector<Particle>::size_type i;
    // masses
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      float m = stof(line);
      particles[i].m = m;
      i++;
    }

    // x
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      float x = stof(line);
      particles[i].x = x;
      i++;
    }

    // y
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      float x = stof(line);
      particles[i].y = x;
      i++;
    }
    
    // z
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      float x = stof(line);
      particles[i].z = x;
      i++;
    } 

    // vx
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      float x = stof(line);
      particles[i].vx = x;
      i++;
    }

    // vy
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      float x = stof(line);
      particles[i].vy = x;
      i++;
    }
    
    // vz
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      float x = stof(line);
      particles[i].vz = x;
      i++;
    }

    // softening
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      float x = stof(line);
      particles[i].softening = x;
      i++;
    }
    // potential
    i = 0;
    while ( getline (myfile,line) && i != particles.size())
    {
      float x = stof(line);
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

void calc_direct_force(std::vector<Particle>& p)
{
  float rx, ry, rz;
  float under;

  for(auto base=p.begin(); base != p.end(); ++base)
  {
    float ax=0, ay=0, az=0;
    for(auto it=p.begin(); it != p.end(); ++it)
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

    int step = base - p.begin();
    cout << "step " << step << " " << ax << " " << ay << " " << az << endl;
    //cout << base - p.begin() << endl; 

    base->ax = ax;
    base->ay = ay;
    base->az = az;

    base->r = sqrt(pow(base->x, 2) + pow(base->y, 2) + pow(base->z, 2));

  }

}


int main(int argc, const char* argv[])
{
	std::string filename = "data.ascii";

	std::vector<Particle> particles = read_particle_data(filename);

  calc_direct_force(particles);
  particles_to_csv(particles, "test.csv");

}