#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include <map>
#include <random>

struct neighbor_site {
  int dx,dy;
  unsigned b;
  double r;
};

struct ipair
{
    bool operator==( const ipair& other ) const {
        if((i==other.i && j==other.j)|| (i==other.j && j==other.i)) return true;
        else return false;
    }
    bool operator<( const ipair& other ) const {
        if(i!=other.i) return i<other.i;
        else return  j<other.j;
    }
    unsigned i, j;
};

struct rate
{
  unsigned i,j;
  double k;
};


class Simulator {
public:

  Simulator(double *_cell, double *_basis, int _nbasis, double _radius,\
            double _jump, double _bond, double *_force, double _penalty);

  void build_neighbour_list();

  // Global index and neighbor_site -> global index of neighbor
  unsigned index(int si, neighbor_site n_s);

  void build(int xsize,int ysize);

  void set_occupations(bool *data);

  unsigned get_nsites();

  double get_time();

  double barrier_function(unsigned nni, unsigned nnj);

  double rate_function(unsigned i, unsigned j);

  double build_rates(std::vector< rate > &rates);

  void build_mobile_set();

  void make_jump(rate jump);

  void run(unsigned steps, bool restart, unsigned seed);

  bool built,occupied;
  unsigned ncells[2], nsites, nbasis, sim_steps;
  double cell[2], force[2], penaltyE, jumpE, bondE, radius,sim_time;
  std::vector<double> basis, jvec, positions;
  std::vector<unsigned int> nnv,nlist_size;
  std::vector<bool> occupation;
  std::set<unsigned int> mobile;
  std::vector< std::vector< neighbor_site > > nlist;
};
