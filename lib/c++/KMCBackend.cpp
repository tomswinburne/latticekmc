#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include <map>
#include <random>

#include "KMCBackend.hpp"

Simulator::Simulator(double *_cell, double *_basis, unsigned _nbasis, \
                      double _radius, double *baren, \
                      unsigned nn, double *_force) {
  nbasis = _nbasis;
  nstates = nn;
  for(unsigned i=0; i<2; i++) {
    cell[i] = _cell[i];
    force[i] = _force[i];
  }

  basis = new double[2*nbasis];
  energies = new double[nstates];
  barriers = new double[nstates*nstates];

  unsigned i,j;
  for(i=0; i<2*nbasis; i++) basis[i] = _basis[i];
  for(i=0; i<nstates; i++) {
    energies[i] = baren[nstates*nstates+i];
    for(j=0; j<nstates; j++) barriers[i*nstates+j] = baren[i*nstates+j];
  }
  bondE =  baren[nstates*nstates+1];
  radius = _radius;
  built = false;
  occupied = false;
  build_neighbour_list();
  sim_time = 0.0;
  sim_steps = 0;
  std::cout<<"Unit Cell: "<<cell[0]<<" "<<cell[1]<<std::endl;
  for(unsigned i=0;i<nstates;i++)
    std::cout<< "E["<<i<<"]="<<energies[i]<<std::endl;
  std::cout<<"Barriers:\n";
  for(unsigned i=0;i<nstates;i++){
    std::cout<<"[ ";
    for(unsigned j=0;j<nstates;j++) std::cout<<barriers[i*nstates+j]<<" ";
    std::cout<<"]\n";
  }
  for(unsigned i=0;i<nbasis;i++) \
    std::cout<<"Basis: "<<basis[2*i]<<" "<<basis[2*i+1]<<std::endl;
  std::cout<<"Cutoff: "<<radius<<std::endl;
  std::cout<<"States: "<<nn<<std::endl;
  std::cout<<"Force: "<<force[0]<<" "<<force[1]<<std::endl;
};

void Simulator::build_neighbour_list(){
  std::vector< neighbor_site > nl;
  int xi,yi;
  unsigned bi,si;
  double x,y,tx=0.,ty=0.;
  neighbor_site n_s;

  nlist_size = new unsigned[nbasis];

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
    nlist_size[si] = unsigned(nl.size());
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
  unsigned i,j;
  ncells[0] = int(xsize/cell[0]);
  ncells[1] = int(ysize/cell[1]);
  nsites = ncells[0] * ncells[1] * nbasis;
  supercell[0] = 1.0*cell[0]*double(ncells[0]);
  supercell[1] = 1.0*cell[1]*double(ncells[1]);

  positions = new double[2*nsites];
  // MUCH faster to hardcode this list
  nnl = new unsigned[nstates*nsites];
  occupation = new bool[nsites];
  nncount = new unsigned[nsites];

  // i = b+nbasis*(x+ncells[0]*y)
  for (i=0; i<nsites; i++) {
    occupation[i] = false;
    nncount[i] = 0;
    b = i%nbasis;
    x = (i/nbasis)%ncells[0];
    y = (i/nbasis)/ncells[0];
    positions[2*i+0] = cell[0] * x + basis[2*b];
    positions[2*i+1] = cell[1] * y + basis[2*b+1];

    // nnl
    for(j=0; j<nlist_size[b];j++) nnl[nstates*i+j] = index(i,nlist[b][j]);
    for(;j<nstates;j++) nnl[nstates*i+j] = nsites;
  }


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

double Simulator::rate_function(unsigned i, unsigned j) {

  unsigned in,jn;
  double dx, barrier, dE; // dE = final-intial

  // Just the onsite considerations- destination has one less neighbor
  dE = energies[nncount[j]-1] - energies[nncount[i]] ;

  // all neighbors of i that are not neighbors of j lose a bond after the jump
  for (in=nstates*i; nnl[in]<nsites; in++) if(occupation[nnl[in]]) {
    for(jn=nstates*j; nnl[jn]<nsites; jn++) if(nnl[jn]==nnl[in]) break;
    if(nnl[jn]==nsites)
      dE += energies[nncount[nnl[in]]-1] - energies[nncount[nnl[in]]];
  }

  // all neighbors of j that are not neighbors of i gain a bond after the jump
  for (jn=nstates*j; nnl[jn]<nsites;jn++) if(occupation[nnl[jn]]) {
    for (in=nstates*i; nnl[in]<nsites; in++) if(nnl[jn]==nnl[in]) break;
    if(nnl[in]==nsites)
      dE += energies[nncount[nnl[jn]]+1] - energies[nncount[nnl[jn]]];
  }
  
  /*
  for(in=nstates*i;nnl[in]<nsites;in++)
    for(jn=nstates*j;nnl[jn]<nsites;jn++)
      if(nnl[jn]!=nnl[in]) {
        dE += double(occupation[nnl[in]]) * (energies[nncount[nnl[in]]-1] - energies[nncount[nnl[in]]]);
        dE += double(occupation[nnl[jn]]) * (energies[nncount[nnl[jn]]+1] - energies[nncount[nnl[jn]]]);
      }
  */
  // barrier matrix
  barrier = barriers[nncount[i]*nstates+nncount[j]-1];

  for(int k=0;k<2;k++) {
    dx = positions[2*j+k]-positions[2*i+k];
    dx -= round(dx/supercell[k])*supercell[k];
    dE += -force[k] * dx;
    barrier += -0.5 * force[k] * dx;
  }

  return exp( -1.0 * ( std::max(dE,0.) + barrier ) );
};

double Simulator::build_rates(std::vector< rate > &rates) {
  rate rate;
  double total_rate = 0.;
  unsigned ni,k;
  rates.clear();
  for (auto mi:mobile) {
    for(k=0;k<nstates-1;k++) {//ni=nstates*mi+k;nnl[ni]!=nsites;ni++) {
      ni=nstates*mi+k;
      if(occupation[nnl[ni]]) continue;
      rate.i = mi;
      rate.j = nnl[ni];
      rate.k = rate_function(mi,nnl[ni]);
      rates.push_back(rate);
      total_rate += rate.k;
    }
  }
  return total_rate;
};

void Simulator::build_mobile_set(){
  // cycle through and build mobile list and nncount
  mobile.clear();
  unsigned i,ni;
  for(i=0;i<nsites;i++) {
    nncount[i] = 0;
    for(ni=nstates*i;nnl[ni]!=nsites;ni++)
      nncount[i] += unsigned(occupation[nnl[ni]]);
    if ( (nncount[i]<nlist_size[i%nbasis]) && occupation[i] ) mobile.insert(i);
  }
};

void Simulator::make_jump(rate jump) {
  // got through and recalculate nn of nlist for orig and dest

  // origin clearly not occupied and not mobile
  occupation[jump.i] = false;
  mobile.erase(jump.i);

  // destination clearly occupied and mobile
  occupation[jump.j] = true;
  mobile.insert(jump.j);

  unsigned i,ni,nni,ii,k,kk;
  // recount neighbors for origin and neighbor sites

  for(ii=0;ii<2;ii++) {
    i = jump.j*ii + jump.i*(1-ii);
    nncount[i]=0;
    for(k=0;k<nstates-1;k++) {
      ni=nstates*i+k;
      nncount[i] += unsigned(occupation[nnl[ni]]);
      nncount[nnl[ni]] = 0;
      for(kk=0;kk<nstates-1;kk++) {
        nni=nstates*nnl[ni]+kk;
        nncount[nnl[ni]] += unsigned(occupation[nnl[nni]]);
      }
      if( (nncount[nnl[ni]]<nlist_size[nnl[ni]%nbasis]) && occupation[nnl[ni]])
        mobile.insert(nnl[ni]);
      else mobile.erase(nnl[ni]);
    }
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
  total_rate = 0.;
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

  void open_sim(double *lattice, double *basis, unsigned nbasis, double radius,\
                double *barriers, unsigned nn, double *force, void ** ptr) {
    Simulator *kmcsim = \
      new Simulator(lattice,basis,nbasis,radius,barriers,nn,force);
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
    for (i=0;i<2*kmcsim->nsites;i++) data[i] = kmcsim->positions[i];
  };

  void set_occupations(void *ptr, bool *data) {
    Simulator *kmcsim = (Simulator *) ptr;
    kmcsim->set_occupations(data);
  };

  void get_occupations(void *ptr, bool *data) {
    Simulator *kmcsim = (Simulator *) ptr;
    unsigned i;
    for (i=0;i<kmcsim->nsites;i++) data[i] = kmcsim->occupation[i];
  };

  void get_neigh_count(void *ptr, unsigned *data) {
    Simulator *kmcsim = (Simulator *) ptr;
    unsigned i;
    for (i=0;i<kmcsim->nsites;i++) data[i] = kmcsim->nncount[i];
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
