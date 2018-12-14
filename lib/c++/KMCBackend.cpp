#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include <map>
#include <random>

#include "KMCBackend.hpp"

Simulator::Simulator(double *_cell, double *_basis, int _nbasis, \
  double _radius, double _jump, double _bond, double *_force, double _penalty) {
  for(int i=0;i<2;i++) cell[i] = _cell[i];
  for(int i=0;i<2;i++) force[i] = _force[i];
  for(int i=0;i<2*_nbasis;i++) basis.push_back(_basis[i]);
  nbasis = _nbasis;
  jumpE = _jump;
  bondE = _bond;
  radius = _radius;
  penaltyE = _penalty;
  built = false;
  occupied = false;
  build_neighbour_list();
  sim_time = 0.;
  sim_steps = 0;

  std::cout<<"cell: "<<cell[0]<<" "<<cell[1]<<std::endl;
  for(unsigned i=0;i<nbasis;i++) \
    std::cout<<"basis: "<<basis[2*i]<<" "<<basis[2*i+1]<<std::endl;
  std::cout<<"rad: "<<radius<<std::endl;
  std::cout<<"jump: "<<jumpE<<std::endl;
  std::cout<<"bond: "<<bondE<<std::endl;
  std::cout<<"penalty: "<<penaltyE<<std::endl;
  std::cout<<"force: "<<force[0]<<" "<<force[1]<<std::endl;
};

void Simulator::build_neighbour_list(){
  std::vector< neighbor_site > nl;
  int xi,yi;
  unsigned bi,si;
  double x,y,tx=0.,ty=0.;
  neighbor_site n_s;

  // look at self and 8 surrounding cells
  // 6 7 8     2  3  4
  // 3 4 5 -> -1  0  1
  // 0 1 2    -4 -3 -2

  for (si=0; si<nbasis; si++) {
    nl.clear();
    tx=0.; ty=0.;
    for(xi=-1;xi<2;xi++) for(yi=-1;yi<2;yi++) for(bi=0;bi<nbasis;bi++) {
      if(xi==0&&yi==0&&bi==si) continue;
      x = cell[0]*xi + basis[2*bi]-basis[2*si];
      y = cell[1]*yi + basis[2*bi+1]-basis[2*si+1];
      if (x*x+y*y <= radius*radius) {
        n_s.dx = xi;
        n_s.dy = yi;
        n_s.b = bi;
        n_s.r = sqrt(x*x+y*y);
        nl.push_back(n_s);
        tx += x;
        ty += y;
      }
    }
    nlist.push_back(nl);
    nlist_size.push_back(nl.size());
    std::cout<<"Basis atom "<<si<<" has "<<nl.size()<<" neighbors"<<std::endl;
  }
};

// Global index and neighbor_site -> global index of neighbor
unsigned Simulator::index(int si, neighbor_site n_s) {
  // i = b+nbasis*(x+ncells[0]*y)

  unsigned xi = (si / nbasis) % ncells[0];
  unsigned yi = (si / nbasis) / ncells[0];
  unsigned uxi = (xi+n_s.dx+ncells[0])%ncells[0];
  unsigned uyi = (yi+n_s.dy+ncells[1])%ncells[1];
  unsigned result = n_s.b + nbasis*(uxi+uyi*ncells[0]);
  return result;
};

void Simulator::build(int xsize,int ysize){
  int b,x,y;
  std::pair<double,double> apos;
  ncells[0] = int(xsize/cell[0]);
  ncells[1] = int(ysize/cell[1]);
  nsites = ncells[0] * ncells[1] * nbasis;

  // i = b+nbasis*(x+ncells[0]*y)
  for (unsigned i=0; i<nsites; i++) {
    b = i%nbasis;
    x = (i/nbasis)%ncells[0];
    y = (i/nbasis)/ncells[0];
    positions.push_back( cell[0] * x + basis[2*b]);
    positions.push_back( cell[1] * y + basis[2*b+1]);
  }
  occupation = std::vector<bool>(nsites,false);
  nnv = std::vector<unsigned>(nsites,0);
  built = true;
};

void Simulator::set_occupations(bool *data) {
  if(!built) {
    std::cout<<"Please build lattice"<<std::endl;
    return;
  }
  for(unsigned i=0;i<nsites;i++) occupation[i] = data[i];
  build_mobile_set();
  occupied = true;
};

unsigned Simulator::get_nsites(){
  return nsites;
};

double Simulator::get_time(){
  return sim_time;
};

double Simulator::barrier_function(unsigned nni, unsigned nnj) {
  return -2.*bondE*(double(nnj-1)- double(nni));
};

double Simulator::rate_function(unsigned i, unsigned j) {
  // will become more complicated
  double dx, dE = barrier_function(nnv[i],nnv[j]);
  for(int k=0;k<2;k++) {
    dx = positions[2*j+k]-positions[2*i+k];
    dx -= round(dx/cell[k]/ncells[k])*cell[k]*ncells[k];
    dE += -force[k]*dx;
  }
  double localjumpE = jumpE;
  if(nnv[i]==0 && nnv[j]>1) localjumpE += penaltyE;
  return exp( -( std::max(dE,0.) + localjumpE ) );
};

double Simulator::build_rates(std::vector< rate > &rates) {
  rate rate;
  double total_rate=0.;
  rates.clear();
  for (auto mi:mobile) {
    for (auto ni:nlist[mi%nbasis]) {
      rate.i = mi;
      rate.j = index(mi,ni);
      if (!occupation[rate.j]) {
        rate.k = rate_function(rate.i,rate.j);
        rates.push_back(rate);
        total_rate += rate.k;
      }
    }
  }
  return total_rate;
};

void Simulator::build_mobile_set(){
  // cycle through and build mobile list and nnv
  mobile.clear();
  for(unsigned i=0;i<nsites;i++) {
    nnv[i] = 0;
    for(auto ni: nlist[i%nbasis]) nnv[i] += unsigned(occupation[index(i,ni)]);
    if(nnv[i]<nlist_size[i%nbasis] && occupation[i]) mobile.insert(i);
  }
};

void Simulator::make_jump(rate jump) {
  // got through and recalculate nn of nlist for orig and dest

  // origin clearly not occupied
  occupation[jump.i] = false;
  mobile.erase(jump.i);

  // destination clearly occupied and mobile
  occupation[jump.j] = true;
  mobile.insert(jump.j);

  // newly mobile nn to origin ?
  nnv[jump.i]=0;
  for(auto ni: nlist[jump.i%nbasis]) {
    unsigned nind = index(jump.i,ni);
    mobile.erase(nind);
    if(occupation[nind]) nnv[jump.i]+=1;
    nnv[nind] = 0;
    for(auto nni: nlist[nind%nbasis])
      if(occupation[index(nind,nni)]) nnv[nind]+=1;
    if(nnv[nind]<nlist_size[nind%nbasis] && occupation[nind])
      mobile.insert(nind);
  }

  // newly immobile nn to destination?
  nnv[jump.j]=0;
  for(auto ni: nlist[jump.j%nbasis]) {
    unsigned nind = index(jump.j,ni);
    if(occupation[nind]) nnv[jump.j]+=1;
    mobile.erase(nind);
    nnv[nind] = 0;
    for(auto nni: nlist[nind%nbasis])
      if(occupation[index(nind,nni)]) nnv[nind]+=1;
    if(nnv[nind]<nlist_size[nind%nbasis] && occupation[nind])
      mobile.insert(nind);
  }
};

void Simulator::run(unsigned steps, bool restart, unsigned seed) {
  if(!built || !occupied) {
    std::cout<<"Please build / occupy lattice"<<std::endl;
    return;
  }
  if(restart) {
    sim_time = 0.;
    sim_steps = 0;
  }



  std::default_random_engine generator(seed);
  std::uniform_real_distribution<double> distribution(1.0e-10,1.0);
  std::vector<rate> rates;
  double total_rate,sel_rate,target;

  // One O(N) scan
  build_mobile_set();

  for(unsigned step=0; step<steps; step++,sim_steps++) {

    total_rate = build_rates(rates);
    target = distribution(generator) * total_rate;
    sel_rate=0.;
    for(auto r: rates) {
      sel_rate += r.k;
      if(sel_rate>target) {
        make_jump(r);
        break;
      }
    }
    sim_time += -log(distribution(generator)) / total_rate;
  }
  if(restart) std::cout<<"Ran "<<sim_steps<<" steps; KMC time="<<sim_time<<"ps; ";
};


extern "C" {

  void open_sim(double *lattice, double *basis, int nbasis, double radius,\
                double jump, double bond, double *force, double penalty,\
                void ** ptr) {

    Simulator *kmcsim = \
      new Simulator(lattice,basis,nbasis,radius,jump,bond,force,penalty);


    *ptr = (void *) kmcsim;
  };

  void build(void *ptr, int xsize, int ysize){
    Simulator *kmcsim = (Simulator *) ptr;
    kmcsim->build(xsize,ysize);
  };

  unsigned nsites(void *ptr) {
    Simulator *kmcsim = (Simulator *) ptr;
    unsigned result = kmcsim->get_nsites();
    return result;
  };

  void get_positions(void *ptr, double *data) {
    Simulator *kmcsim = (Simulator *) ptr;
    unsigned i;
    for (i=0;i<2*kmcsim->nsites;i++) data[i] = kmcsim->positions.at(i);
  };

  void set_occupations(void *ptr, bool *data) {
    Simulator *kmcsim = (Simulator *) ptr;
    kmcsim->set_occupations(data);
  };

  void get_occupations(void *ptr, bool *data) {
    Simulator *kmcsim = (Simulator *) ptr;
    unsigned i;
    for (i=0;i<kmcsim->nsites;i++) data[i] = kmcsim->occupation.at(i);
  };

  void get_neigh_count(void *ptr, unsigned *data) {
    Simulator *kmcsim = (Simulator *) ptr;
    unsigned i;
    for (i=0;i<kmcsim->nsites;i++) data[i] = kmcsim->nnv.at(i);
  };

  void run(void *ptr, unsigned steps, bool verbose, unsigned seed) {
    Simulator *kmcsim = (Simulator *) ptr;
    kmcsim->run(steps,verbose,seed);
  };

  double get_time(void *ptr) {
    Simulator *kmcsim = (Simulator *) ptr;
    double result = kmcsim->get_time();
    return result;
  };


  //void tester(void *ptr, double bb) {
    //Simulator *lmp = (Simulator *) ptr;
    //lmp->tester(bb);
  //};


};
