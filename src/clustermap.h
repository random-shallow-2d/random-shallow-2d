#ifndef CLUSTERMAP_H
#define CLUSTERMAP_H

#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>
#include "itensor/all.h"
using namespace itensor;

typedef std::array<double, 4> PauliDist;
typedef std::map<int, std::map<int, std::map<int, std::complex<double> > > >
    UnitaryData;

struct EntanglementData {
  std::map<int, double> s;
  std::map<int, std::vector<double> > svals;
};

class ClusterMap {
  std::map<int, ITensor> data;
  std::string name;
  int iteration_num;

 public:
  ClusterMap(std::string l, int N);
  ClusterMap(std::string l, int N, bool product);
  ClusterMap(MPS& M, std::string na, int it_num);


  int get_size();
  ITensor& get_tensor(int pos);
  std::string get_name();
  int get_iteration_num();
  EntanglementData get_entanglement();


  ClusterMap& sphere_channel(double alpha_lb, double alpha_ub);
  ClusterMap& effective_channel(double x_lb, double x_ub);

  ClusterMap& sphere_measure(double alpha_lb, double alpha_ub, int pos,
                             Index i);
  ClusterMap& weak_measure(double x_lb, double x_ub, int pos);
  ClusterMap& apply_sphere_U(double alpha_lb, double alpha_ub, int pos,
                             Index i);

  ClusterMap& measure(std::map<int, Index>& index_to_meas);
  ClusterMap& full_weak_measure(double x_lb, double x_ub);

  MPS to_mps();

  ClusterMap& brick_channel(UnitaryData& udat);
  ClusterMap& cz_brick_channel();
  ClusterMap& apply_brick_U(int pos, Index i);
};

#endif
