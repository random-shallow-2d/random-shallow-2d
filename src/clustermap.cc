#include "clustermap.h"

ClusterMap::ClusterMap(std::string l, int N) {
  iteration_num = 0;
  name = l;

  auto i = Index(name + " " + std::to_string(1), 2, Site);
  auto right = Index("bond " + std::to_string(1), 2, Link);
  data[1] = ITensor(i, right);
  data[1].set(i(1), right(1), 1 / sqrt(2));
  data[1].set(i(1), right(2), 0);
  data[1].set(i(2), right(1), 0);
  data[1].set(i(2), right(2), 1 / sqrt(2));

  for (int ind = 2; ind <= N - 1; ind++) {
    auto left = right;
    auto new_right = Index("bond " + std::to_string(ind), 2, Link);
    i = Index(name + " " + std::to_string(ind), 2, Site);

    data[ind] = ITensor(i, left, new_right);
    data[ind].set(i(1), left(1), new_right(1), 1 / sqrt(2));
    data[ind].set(i(1), left(1), new_right(2), 0);
    data[ind].set(i(1), left(2), new_right(1), 1 / sqrt(2));
    data[ind].set(i(1), left(2), new_right(2), 0);

    data[ind].set(i(2), left(1), new_right(1), 0);
    data[ind].set(i(2), left(1), new_right(2), 1 / sqrt(2));
    data[ind].set(i(2), left(2), new_right(1), 0);
    data[ind].set(i(2), left(2), new_right(2), -1 / sqrt(2)); right = new_right;
  }

  i = Index(name + " " + std::to_string(N), 2, Site);
  auto left = right;
  data[N] = ITensor(i, left);
  data[N].set(i(1), left(1), 1 / sqrt(2));
  data[N].set(i(1), left(2), 1 / sqrt(2));
  data[N].set(i(2), left(1), 1 / sqrt(2));
  data[N].set(i(2), left(2), -1 / sqrt(2));
}

ClusterMap::ClusterMap(std::string l, int N, bool product) {
  iteration_num = 0;
  name = l;

  if (product == false) {
    auto i = Index(name + " " + std::to_string(1), 2, Site);
    auto right = Index("bond " + std::to_string(1), 2, Link);
    data[1] = ITensor(i, right);
    data[1].set(i(1), right(1), 1 / sqrt(2));
    data[1].set(i(1), right(2), 0);
    data[1].set(i(2), right(1), 0);
    data[1].set(i(2), right(2), 1 / sqrt(2));

    for (int ind = 2; ind <= N - 1; ind++) {
      auto left = right;
      auto new_right = Index("bond " + std::to_string(ind), 2, Link);
      i = Index(name + " " + std::to_string(ind), 2, Site);

      data[ind] = ITensor(i, left, new_right);
      data[ind].set(i(1), left(1), new_right(1), 1 / sqrt(2));
      data[ind].set(i(1), left(1), new_right(2), 0);
      data[ind].set(i(1), left(2), new_right(1), 1 / sqrt(2));
      data[ind].set(i(1), left(2), new_right(2), 0);

      data[ind].set(i(2), left(1), new_right(1), 0);
      data[ind].set(i(2), left(1), new_right(2), 1 / sqrt(2));
      data[ind].set(i(2), left(2), new_right(1), 0);
      data[ind].set(i(2), left(2), new_right(2), -1 / sqrt(2));

      right = new_right;
    }

    i = Index(name + " " + std::to_string(N), 2, Site);
    auto left = right;
    data[N] = ITensor(i, left);
    data[N].set(i(1), left(1), 1 / sqrt(2));
    data[N].set(i(1), left(2), 1 / sqrt(2));
    data[N].set(i(2), left(1), 1 / sqrt(2));
    data[N].set(i(2), left(2), -1 / sqrt(2));
  }
  if (product == true) {
    for (int ind = 1; ind <= N; ind++) {
      auto i = Index(name + " " + std::to_string(ind), 2, Site);

      data[ind] = ITensor(i);
      data[ind].set(i(1), 1);
      data[ind].set(i(2), 0);
    }
  }
}

ClusterMap::ClusterMap(MPS& M, std::string na, int it_num) {
  iteration_num = it_num;
  name = na;
  for (int i = 1; i <= M.N(); i++) {
    data[i] = M.A(i);
  }
}

int ClusterMap::get_size() { return (this->data).size(); }

ITensor& ClusterMap::get_tensor(int pos) { return (this->data)[pos]; }

std::string ClusterMap::get_name() { return (this->name); }

int ClusterMap::get_iteration_num() { return (this->iteration_num); }

EntanglementData ClusterMap::get_entanglement() {
  int N = this->get_size();

  std::map<int, double> entanglement_data;
  std::map<int, std::vector<double> > svs;
  EntanglementData out;

  std::map<int, ITensor> lambdas;

  ITensor conj = data[1];
  conj.conj();
  conj.prime(Link);

  lambdas[1] = data[1] * conj;

  for (int m = 2; m <= N; m++) {
    conj = data[m];
    conj.conj();
    conj.prime(Link);
    lambdas[m] = lambdas[m - 1] * conj * data[m];
  }

  for (int m = 1; m <= N - 1; m++) {
    int D = lambdas[m].inds().index(1).m();
    entanglement_data[m] = 0;
    for (int p = 1; p <= D; p++) {
      double val = lambdas[m].real(p, p);
      svs[m].push_back(val);
      entanglement_data[m] += -val * log2(val);
    }
  }

  out.s = entanglement_data;
  out.svals = svs;
  return out;
}

ClusterMap& ClusterMap::sphere_channel(double alpha_lb, double alpha_ub) {
  iteration_num = iteration_num + 1;
  int N = this->get_size();

  std::map<int, Index> index_to_meas;
  std::map<int, Index> new_phys_index;

  ClusterMap new_cluster("it:" + std::to_string(iteration_num) + " ", N);

  for (int pos = 1; pos <= N; pos++) {
    ITensor T0 = data[pos];
    ITensor T1 = new_cluster.get_tensor(pos);
    Index i;
    auto k = Index("inter bond " + std::to_string(pos), 2, Link);
    std::complex<double> im(0, 1);

    for (auto& I : T0.inds()) {
      if (I.type() == Site) {
        i = I;
      }
    }
    auto ii = Index(name + " " + std::to_string(pos), 2, Site);

    auto CZ_T0 = ITensor(k, i, ii);

    CZ_T0.set(k(1), i(1), ii(1), im);
    CZ_T0.set(k(1), i(2), ii(2), -im);
    CZ_T0.set(k(1), i(1), ii(2), 0);
    CZ_T0.set(k(1), i(2), ii(1), 0);

    CZ_T0.set(k(2), i(1), ii(1), sqrt(2));
    CZ_T0.set(k(2), i(2), ii(2), 0);
    CZ_T0.set(k(2), i(1), ii(2), 0);
    CZ_T0.set(k(2), i(2), ii(1), 0);

    T0 = T0 * CZ_T0;

    index_to_meas[pos] = ii;

    for (auto& I : T1.inds()) {
      if (I.type() == Site) i = I;
    }
    ii =
        Index("it:" + std::to_string(iteration_num) + " " + std::to_string(pos),
              2, Site);
    auto CZ_T1 = ITensor(k, i, ii);

    CZ_T1.set(k(1), i(1), ii(1), im);
    CZ_T1.set(k(1), i(2), ii(2), -im);
    CZ_T1.set(k(1), i(1), ii(2), 0);
    CZ_T1.set(k(1), i(2), ii(1), 0);

    CZ_T1.set(k(2), i(1), ii(1), sqrt(2));
    CZ_T1.set(k(2), i(2), ii(2), 0);
    CZ_T1.set(k(2), i(1), ii(2), 0);
    CZ_T1.set(k(2), i(2), ii(1), 0);

    new_phys_index[pos] = ii;
    T1 = T1 * CZ_T1;

    data[pos] = T0 * T1;
  }

  for (int pos = 1; pos <= (N - 1); pos++) {
    auto T0 = data[pos];
    auto T1 = data[pos + 1];

    std::vector<Index> common_inds;
    for (auto& I : T0.inds()) {
      if (hasindex(T1, I)) {
        common_inds.push_back(I);
      }
    }
    auto C = combiner(common_inds);
    data[pos] = C * data[pos];
    data[pos + 1] = C * data[pos + 1];
  }

  for (int pos = 1; pos <= N; pos++) {
    *this = this->apply_sphere_U(alpha_lb, alpha_ub, pos, index_to_meas[pos]);
  }
  for (int pos = 1; pos <= N; pos++) {
    for (auto& I : data[pos].inds()) {
      if ((I.type() == Site) && (I != new_phys_index[pos])) {
        index_to_meas[pos] = I;
      }
    }
  }

  *this = this->measure(index_to_meas);

  return *this;
}

ClusterMap& ClusterMap::effective_channel(double x_lb, double x_ub) {
  int N = this->get_size();
  this->full_weak_measure(x_lb, x_ub);
  ITensor CZ0;
  ITensor CZ1;
  Index i;
  Index ii;
  std::complex<double> im(0.0, 1.0);

  std::map<int, Index> phys_inds;
  for (int pos = 1; pos <= N; pos++) {
    for (auto& I : data[pos].inds()) {
      if (I.type() == Site) {
        phys_inds[pos] = I;
      }
    }
  }

  for (int pos = 1; pos <= N - 1; pos++) {
    auto ii = Index("Physical " + std::to_string(pos), 2, Site);
    auto k = Index("CZ bond " + std::to_string(pos), 2, Link);
    i = phys_inds[pos];
    CZ0 = ITensor(phys_inds[pos], ii, k);
    CZ0.set(k(1), i(1), ii(1), im);
    CZ0.set(k(1), i(2), ii(2), -im);
    CZ0.set(k(1), i(1), ii(2), 0);
    CZ0.set(k(1), i(2), ii(1), 0);

    CZ0.set(k(2), i(1), ii(1), sqrt(2));
    CZ0.set(k(2), i(2), ii(2), 0);
    CZ0.set(k(2), i(1), ii(2), 0);
    CZ0.set(k(2), i(2), ii(1), 0);

    data[pos] = data[pos] * CZ0;
    phys_inds[pos] = ii;

    i = phys_inds[pos + 1];
    ii = Index("Physical " + std::to_string(pos + 1), 2, Site);
    auto CZ1 = ITensor(k, i, ii);

    CZ1.set(k(1), i(1), ii(1), im);
    CZ1.set(k(1), i(2), ii(2), -im);
    CZ1.set(k(1), i(1), ii(2), 0);
    CZ1.set(k(1), i(2), ii(1), 0);

    CZ1.set(k(2), i(1), ii(1), sqrt(2));
    CZ1.set(k(2), i(2), ii(2), 0);
    CZ1.set(k(2), i(1), ii(2), 0);
    CZ1.set(k(2), i(2), ii(1), 0);
    data[pos + 1] = data[pos + 1] * CZ1;
    phys_inds[pos + 1] = ii;

    std::vector<Index> common_inds;
    for (auto& I : (data[pos]).inds()) {
      if (hasindex(data[pos + 1], I)) {
        common_inds.push_back(I);
      }
    }
    auto C = combiner(common_inds);
    data[pos] = C * data[pos];
    data[pos + 1] = C * data[pos + 1];
  }

  return *this;
}

ClusterMap& ClusterMap::sphere_measure(double alpha_lb, double alpha_ub,
                                       int pos, Index i) {
  int N = this->get_size();
  std::map<int, ITensor> copy;
  for (int k = 1; k <= N; k++) {
    copy[k] = data[k];
  }

  double glob, psi, phi, chi;
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> distribution1(0.0, 6.28318530718);
  std::uniform_real_distribution<double> distribution2(alpha_lb, alpha_ub);
  std::uniform_real_distribution<double> sample(0.0, 1.0);
  glob = distribution1(generator);
  psi = distribution1(generator);
  chi = distribution1(generator);
  auto X = sqrt(distribution2(generator));
  phi = asin(X);

  auto newI = Index(this->get_name() + " " + std::to_string(pos), 2, Site);

  std::complex<double> im(0.0, 1.0);
  auto U = ITensor(i, newI);
  std::complex<double> val =
      std::exp(im * (glob + psi)) * std::complex<double>(cos(phi), 0.0);
  U.set(i(1), newI(1), val);

  val = std::exp(im * (glob + chi)) * std::complex<double>(sin(phi), 0.0);
  U.set(i(2), newI(1), val);

  val = std::exp(im * (glob - chi)) * std::complex<double>(-sin(phi), 0.0);
  U.set(i(1), newI(2), val);

  val = std::exp(im * (glob - psi)) * std::complex<double>(cos(phi), 0.0);
  U.set(i(2), newI(2), val);

  copy[pos] = copy[pos] * U;

  auto fake_I = Index("Fake", 2, Site);
  auto proj0 = ITensor(newI, fake_I);
  proj0.set(newI(1), fake_I(1), 1);
  proj0.set(newI(2), fake_I(1), 0);
  proj0.set(newI(1), fake_I(2), 0);
  proj0.set(newI(2), fake_I(2), 0);
  auto proj1 = ITensor(newI, fake_I);
  proj1.set(newI(1), fake_I(1), 0);
  proj1.set(newI(2), fake_I(1), 0);
  proj1.set(newI(1), fake_I(2), 0);
  proj1.set(newI(2), fake_I(2), 1);

  copy[pos] = copy[pos] * proj0;

  for (int k = 1; k <= N; k++) {
    std::vector<Index> phys_inds;
    for (auto& I : copy[k].inds()) {
      if (I.type() == Site) {
        phys_inds.push_back(I);
      }
    }
    auto C = combiner(phys_inds, {"IndexType", Site});
    copy[k] = copy[k] * C;
  }

  auto sites = SiteSet(N, 4);
  auto copy_mps = MPS(sites);
  int k = 1;
  for (int k = 1; k <= N; k++) {
    copy_mps.setA(k, copy[k]);
  }

  copy_mps.orthogonalize();
  double prob = std::pow(norm(copy_mps), 2);
  proj0 = ITensor(newI);
  proj0.set(newI(1), 1);
  proj0.set(newI(2), 0);
  proj1 = ITensor(newI);
  proj1.set(newI(1), 0);
  proj1.set(newI(2), 1);
  if (sample(generator) < prob) {
    proj0.set(newI(1), 1 / sqrt(prob));
    data[pos] = data[pos] * U * proj0;
  } else {
    proj1.set(newI(2), 1 / sqrt(1 - prob));
    data[pos] = data[pos] * U * proj1;
  }

  return *this;
}

ClusterMap& ClusterMap::weak_measure(double x_lb, double x_ub, int pos) {
  int N = this->get_size();
  std::map<int, ITensor> copy;
  for (int k = 1; k <= N; k++) {
    copy[k] = data[k];
  }

  double X, phi;
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> distribution1(0.0, 6.28318530718);
  std::uniform_real_distribution<double> distribution2(x_lb, x_ub);
  std::uniform_real_distribution<double> sample(0.0, 1.0);
  phi = distribution1(generator);
  X = distribution2(generator);

  auto newI = Index(this->get_name() + " " + std::to_string(pos), 2, Site);
  auto newI_H = Index(this->get_name() + " " + std::to_string(pos), 2, Site);
  auto newI_phase =
      Index(this->get_name() + " " + std::to_string(pos), 2, Site);

  Index i;
  for (auto& I : data[pos].inds()) {
    if (I.type() == Site) {
      i = I;
    }
  }

  std::complex<double> im(0.0, 1.0);
  auto U = ITensor(i, newI);
  auto H = ITensor(newI_phase, newI_H);
  auto phase = ITensor(newI, newI_phase);
  std::complex<double> val;
  U.set(i(1), newI(1), sqrt((1 + X) / 2.0));
  U.set(i(2), newI(1), 0);
  U.set(i(1), newI(2), 0);
  U.set(i(2), newI(2), sqrt((1 - X) / 2.0));

  H.set(newI_phase(1), newI_H(1), 1 / sqrt(2));
  H.set(newI_phase(1), newI_H(2), 1 / sqrt(2));
  H.set(newI_phase(2), newI_H(1), 1 / sqrt(2));
  H.set(newI_phase(2), newI_H(2), -1 / sqrt(2));

  phase.set(newI_phase(1), newI(1), 1.0);
  phase.set(newI_phase(1), newI(2), 0.0);
  phase.set(newI_phase(2), newI(1), 0.0);
  val = std::exp(im * phi);
  phase.set(newI_phase(2), newI(2), val);

  copy[pos] = copy[pos] * U;

  auto sites = SiteSet(N, 2);
  auto copy_mps = MPS(sites);
  int k = 1;
  for (int k = 1; k <= N; k++) {
    copy_mps.setA(k, copy[k]);
  }

  copy_mps.orthogonalize();
  double prob = std::pow(norm(copy_mps), 2);

  if (sample(generator) < prob) {
    U.set(i(1), newI(1), sqrt((1 + X) / 2.0) / sqrt(prob));
    U.set(i(2), newI(2), sqrt((1 - X) / 2.0) / sqrt(prob));
    data[pos] = data[pos] * U * phase * H;
  } else {
    U.set(i(1), newI(1), sqrt((1 - X) / 2.0) / sqrt(1 - prob));
    U.set(i(2), newI(2), sqrt((1 + X) / 2.0) / sqrt(1 - prob));
    data[pos] = data[pos] * U * phase * H;
  }

  return *this;
}

MPS ClusterMap::to_mps() {
  int N = this->get_size();
  auto sites = SiteSet(N, 2);
  auto psi = MPS(sites);
  for (int pos = 1; pos <= N; pos++) {
    psi.setA(pos, data[pos]);
  }
  return psi;
}

ClusterMap& ClusterMap::measure(std::map<int, Index>& index_to_meas) {
  int N = (this->get_size());
  std::map<int, ITensor> aux_map = data;
  std::map<int, ITensor> aux_combiners;
  std::map<int, ITensor> aux_projs0;
  std::map<int, ITensor> aux_projs1;
  std::chrono::steady_clock::time_point begin, end;

  for (int k = 1; k <= N; k++) {
    std::vector<Index> phys_inds;
    for (auto& I : aux_map[k].inds()) {
      if (I.type() == Site) {
        phys_inds.push_back(I);
        if (I != index_to_meas[k]) {
          Index ii = Index("new", 2, Site);
          aux_projs0[k] = ITensor(index_to_meas[k], I, ii);
          aux_projs1[k] = ITensor(index_to_meas[k], I, ii);

          aux_projs0[k].set(index_to_meas[k](1), I(1), ii(1), 1);
          aux_projs0[k].set(index_to_meas[k](1), I(2), ii(1), 0);
          aux_projs0[k].set(index_to_meas[k](1), I(1), ii(2), 0);
          aux_projs0[k].set(index_to_meas[k](1), I(2), ii(2), 1);
          aux_projs0[k].set(index_to_meas[k](2), I(1), ii(1), 0);
          aux_projs0[k].set(index_to_meas[k](2), I(2), ii(1), 0);
          aux_projs0[k].set(index_to_meas[k](2), I(1), ii(2), 0);
          aux_projs0[k].set(index_to_meas[k](2), I(2), ii(2), 0);

          aux_projs1[k].set(index_to_meas[k](1), I(1), ii(1), 0);
          aux_projs1[k].set(index_to_meas[k](1), I(2), ii(1), 0);
          aux_projs1[k].set(index_to_meas[k](1), I(1), ii(2), 0);
          aux_projs1[k].set(index_to_meas[k](1), I(2), ii(2), 0);
          aux_projs1[k].set(index_to_meas[k](2), I(1), ii(1), 1);
          aux_projs1[k].set(index_to_meas[k](2), I(2), ii(1), 0);
          aux_projs1[k].set(index_to_meas[k](2), I(1), ii(2), 0);
          aux_projs1[k].set(index_to_meas[k](2), I(2), ii(2), 1);
        }
      }
    }
    auto C = combiner(phys_inds, {"IndexType", Site});
    aux_map[k] = aux_map[k] * C;
    aux_projs0[k] = aux_projs0[k] * C;
    aux_projs1[k] = aux_projs1[k] * C;
  }

  auto sites = SiteSet(N, 4);
  auto copy_mps = MPS(sites);
  for (int k = N; k >= 1; k--) {
    copy_mps.setA(k, aux_map[k]);
  }
  for (int pos = 1; pos <= N; pos++) {
    aux_map[pos] = copy_mps.A(pos);
  }

  ITensor copy, copyconj;
  ITensor edge;


  for (int pos = 1; pos <= N; pos++) {
    double prob;
    if (pos == 1) {
      copy = aux_map[pos];
      copy = copy * aux_projs0[pos];
      prob = std::pow(norm(copy), 2);
    } else {
      copy = aux_map[pos];
      copy = copy * aux_projs0[pos];

      copyconj = copy;
      copyconj.conj();

      for (auto& I : copyconj.inds()) {
        if (hasindex(edge, I)) {
          copyconj.prime(I);
        }
      }

      auto T = copyconj * copy * edge;
      prob = T.real();
    }

    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> sample(0.0, 1.0);
    auto proj0 = ITensor(index_to_meas[pos]);
    auto proj1 = ITensor(index_to_meas[pos]);
    proj0.set(index_to_meas[pos](1), 1);
    proj0.set(index_to_meas[pos](2), 0);
    proj1.set(index_to_meas[pos](1), 0);
    proj1.set(index_to_meas[pos](2), 1);

    if (sample(generator) < prob) {
        aux_map[pos] = aux_map[pos] * (aux_projs0[pos] * (1 / sqrt(prob)));
        data[pos] = data[pos] * (proj0 * (1 / sqrt(prob)));
    } else {
        aux_map[pos] = aux_map[pos] * (aux_projs1[pos] * (1 / sqrt(1 - prob)));
        data[pos] = data[pos] * (proj1 * (1 / sqrt(1 - prob)));
    }

    if (pos == 1) {
      copy = aux_map[pos];
      copyconj = copy;
      copyconj.prime(Link);
      copyconj.conj();

      edge = copy * copyconj;
    } else {
      copy = aux_map[pos];
      copyconj = copy;
      copyconj.prime(Link);
      copyconj.conj();

      begin = std::chrono::steady_clock::now();
      begin = std::chrono::steady_clock::now();
      edge = edge * copy * copyconj;
      end = std::chrono::steady_clock::now();
    }
  }
  return *this;
}

ClusterMap& ClusterMap::apply_sphere_U(double alpha_lb, double alpha_ub,
                                       int pos, Index i) {
  int N = this->get_size();

  double glob, psi, phi, chi;
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> distribution1(0.0, 6.28318530718);
  std::uniform_real_distribution<double> distribution2(alpha_lb, alpha_ub);
  std::uniform_real_distribution<double> sample(0.0, 1.0);
  glob = distribution1(generator);
  psi = distribution1(generator);
  chi = distribution1(generator);
  auto X = sqrt(distribution2(generator));
  phi = asin(X);

  auto newI = Index(this->get_name() + " " + std::to_string(pos), 2, Site);

  std::complex<double> im(0.0, 1.0);
  auto U = ITensor(i, newI);
  std::complex<double> val =
      std::exp(im * (glob + psi)) * std::complex<double>(cos(phi), 0.0);
  U.set(i(1), newI(1), val);

  val = std::exp(im * (glob + chi)) * std::complex<double>(sin(phi), 0.0);
  U.set(i(2), newI(1), val);

  val = std::exp(im * (glob - chi)) * std::complex<double>(-sin(phi), 0.0);
  U.set(i(1), newI(2), val);

  val = std::exp(im * (glob - psi)) * std::complex<double>(cos(phi), 0.0);
  U.set(i(2), newI(2), val);

  data[pos] = data[pos] * U;

  return *this;
}

ClusterMap& ClusterMap::full_weak_measure(double x_lb, double x_ub) {
  int N = this->get_size();
  std::map<int, ITensor> aux_map;

  auto sites = SiteSet(N, 2);
  auto copy_mps = MPS(sites);
  for (int k = N; k >= 1; k--) {
    copy_mps.setA(k, data[k]);
  }

  for (int k = 1; k <= N; k++) {
    aux_map[k] = copy_mps.A(k);
  }

  ITensor copy, copyconj;
  ITensor edge;
  for (int pos = 1; pos <= N; pos++) {
    double prob;
    double X, phi;
    std::random_device rand_dev;
    std::mt19937 generator(rand_dev());
    std::uniform_real_distribution<double> distribution1(0.0, 6.28318530718);
    std::uniform_real_distribution<double> distribution2(x_lb, x_ub);
    std::uniform_real_distribution<double> sample(0.0, 1.0);
    phi = distribution1(generator);
    X = distribution2(generator);

    auto newI = Index(this->get_name() + " " + std::to_string(pos), 2, Site);
    auto newI_H = Index(this->get_name() + " " + std::to_string(pos), 2, Site);
    auto newI_phase =
        Index(this->get_name() + " " + std::to_string(pos), 2, Site);

    Index i;
    for (auto& I : data[pos].inds()) {
      if (I.type() == Site) {
        i = I;
      }
    }

    auto U = ITensor(i, newI);
    auto H = ITensor(newI_phase, newI_H);
    auto phase = ITensor(newI, newI_phase);

    std::complex<double> im(0.0, 1.0);
    std::complex<double> val;
    U.set(i(1), newI(1), sqrt((1 + X) / 2.0));
    U.set(i(2), newI(1), 0);
    U.set(i(1), newI(2), 0);
    U.set(i(2), newI(2), sqrt((1 - X) / 2.0));

    H.set(newI_phase(1), newI_H(1), 1 / sqrt(2));
    H.set(newI_phase(1), newI_H(2), 1 / sqrt(2));
    H.set(newI_phase(2), newI_H(1), 1 / sqrt(2));
    H.set(newI_phase(2), newI_H(2), -1 / sqrt(2));

    phase.set(newI_phase(1), newI(1), 1.0);
    phase.set(newI_phase(1), newI(2), 0.0);
    phase.set(newI_phase(2), newI(1), 0.0);
    val = std::exp(im * phi);
    phase.set(newI_phase(2), newI(2), val);

    if (pos == 1) {
      copy = aux_map[pos];
      copy = copy * U;
      prob = std::pow(norm(copy), 2);
    } else {
      copy = aux_map[pos];
      copy = copy * U;
      copyconj = copy;
      copyconj.conj();
      for (auto& I : copyconj.inds()) {
        if (hasindex(edge, I)) {
          copyconj.prime(I);
        }
      }
      auto T = edge * copy * copyconj;
      prob = T.real();
    }

    if (sample(generator) < prob) {
      U.set(i(1), newI(1), sqrt((1 + X) / 2.0) / sqrt(prob));
      U.set(i(2), newI(2), sqrt((1 - X) / 2.0) / sqrt(prob));
      aux_map[pos] = aux_map[pos] * U * phase * H;
      data[pos] = data[pos] * U * phase * H;
    } else {
      U.set(i(1), newI(1), sqrt((1 - X) / 2.0) / sqrt(1 - prob));
      U.set(i(2), newI(2), sqrt((1 + X) / 2.0) / sqrt(1 - prob));
      aux_map[pos] = aux_map[pos] * U * phase * H;
      data[pos] = data[pos] * U * phase * H;
    }

    if (pos == 1) {
      copy = aux_map[pos];
      copyconj = copy;
      copyconj.prime(Link);
      copyconj.conj();
      edge = copy * copyconj;
    } else {
      copy = aux_map[pos];
      copyconj = copy;
      copyconj.prime(Link);
      copyconj.conj();
      edge = edge * copy * copyconj;
    }
  }
  return *this;
}

ClusterMap& ClusterMap::brick_channel(UnitaryData& udat) {
  iteration_num = iteration_num + 1;
  int N = this->get_size();
  int M = udat.size();

  if (M < (N - 1) + 2 * ((N - 1) / 8.0)) {
    std::cout << "ERROR: Not enough unitary data for applying channel."
              << std::endl;
  }

  ClusterMap new_cluster("it:" + std::to_string(iteration_num) + " ", N, true);

  auto psi = new_cluster.to_mps();

  for (int pos = 1; pos <= N - 1; pos++) {
    psi.position(pos);
    ITensor T0 = psi.A(pos);
    ITensor T1 = psi.A(pos + 1);
    Index i, j;

    for (auto& I : T0.inds()) {
      if (I.type() == Site) {
        i = I;
      }
    }
    auto ii = prime(i);

    for (auto& I : T1.inds()) {
      if (I.type() == Site) j = I;
    }
    auto jj = prime(j);
    auto U = ITensor(i, j, ii, jj);

    U.set(i(1), j(1), ii(1), jj(1), udat[pos][0][0]);
    U.set(i(1), j(2), ii(1), jj(1), udat[pos][0][1]);
    U.set(i(2), j(1), ii(1), jj(1), udat[pos][0][2]);
    U.set(i(2), j(2), ii(1), jj(1), udat[pos][0][3]);

    U.set(i(1), j(1), ii(1), jj(2), udat[pos][1][0]);
    U.set(i(1), j(2), ii(1), jj(2), udat[pos][1][1]);
    U.set(i(2), j(1), ii(1), jj(2), udat[pos][1][2]);
    U.set(i(2), j(2), ii(1), jj(2), udat[pos][1][3]);

    U.set(i(1), j(1), ii(2), jj(1), udat[pos][2][0]);
    U.set(i(1), j(2), ii(2), jj(1), udat[pos][2][1]);
    U.set(i(2), j(1), ii(2), jj(1), udat[pos][2][2]);
    U.set(i(2), j(2), ii(2), jj(1), udat[pos][2][3]);

    U.set(i(1), j(1), ii(2), jj(2), udat[pos][3][0]);
    U.set(i(1), j(2), ii(2), jj(2), udat[pos][3][1]);
    U.set(i(2), j(1), ii(2), jj(2), udat[pos][3][2]);
    U.set(i(2), j(2), ii(2), jj(2), udat[pos][3][3]);

    auto wf = T0 * T1;
    wf *= U;
    wf.noprime();

    ITensor S, V;
    ITensor W = psi.A(pos);
    svd(wf, W, S, V, {"Cutoff=", 1E-18});
    psi.setA(pos, W);
    psi.setA(pos + 1, S * V);
  }
  psi.orthogonalize();
  normalize(psi);

  new_cluster = ClusterMap(psi, new_cluster.get_name(), iteration_num);
  for (int pos = 1; pos <= N; pos++) {
  }

  std::map<int, Index> index_to_meas;
  std::map<int, Index> new_phys_index;
  int counter = 0;

  for (int pos = 1; pos <= N; pos++) {
    ITensor T0 = data[pos];
    ITensor T1 = new_cluster.get_tensor(pos);
    Index i, j;
    for (auto& I : T0.inds()) {
      if (I.type() == Site) {
        i = I;
      }
    }
    for (auto& I : T1.inds()) {
      if (I.type() == Site) {
        j = I;
      }
    }


    if (iteration_num % 2 == 1) {
      if (((pos - 1) % 8 == 1) || ((pos - 1) % 8 == 3)) {

        auto ii = prime(i);
        auto jj = prime(j);
        auto U = ITensor(i, j, ii, jj);


        U.set(i(1), j(1), ii(1), jj(1), udat[N + counter][0][0]);
        U.set(i(1), j(2), ii(1), jj(1), udat[N + counter][0][1]);
        U.set(i(2), j(1), ii(1), jj(1), udat[N + counter][0][2]);
        U.set(i(2), j(2), ii(1), jj(1), udat[N + counter][0][3]);

        U.set(i(1), j(1), ii(1), jj(2), udat[N + counter][1][0]);
        U.set(i(1), j(2), ii(1), jj(2), udat[N + counter][1][1]);
        U.set(i(2), j(1), ii(1), jj(2), udat[N + counter][1][2]);
        U.set(i(2), j(2), ii(1), jj(2), udat[N + counter][1][3]);

        U.set(i(1), j(1), ii(2), jj(1), udat[N + counter][2][0]);
        U.set(i(1), j(2), ii(2), jj(1), udat[N + counter][2][1]);
        U.set(i(2), j(1), ii(2), jj(1), udat[N + counter][2][2]);
        U.set(i(2), j(2), ii(2), jj(1), udat[N + counter][2][3]);

        U.set(i(1), j(1), ii(2), jj(2), udat[N + counter][3][0]);
        U.set(i(1), j(2), ii(2), jj(2), udat[N + counter][3][1]);
        U.set(i(2), j(1), ii(2), jj(2), udat[N + counter][3][2]);
        U.set(i(2), j(2), ii(2), jj(2), udat[N + counter][3][3]);

        data[pos] = U * (T0 * T1);
        data[pos].noprime();
        counter++;
      } else {
        data[pos] = T0 * T1;
      }
    } else {
      if (((pos - 1) % 8 == 5) || ((pos - 1) % 8 == 7)) {

        auto ii = prime(i);
        auto jj = prime(j);
        auto U = ITensor(i, j, ii, jj);
        U.set(i(1), j(1), ii(1), jj(1), udat[N + counter][0][0]);
        U.set(i(1), j(2), ii(1), jj(1), udat[N + counter][0][1]);
        U.set(i(2), j(1), ii(1), jj(1), udat[N + counter][0][2]);
        U.set(i(2), j(2), ii(1), jj(1), udat[N + counter][0][3]);

        U.set(i(1), j(1), ii(1), jj(2), udat[N + counter][1][0]);
        U.set(i(1), j(2), ii(1), jj(2), udat[N + counter][1][1]);
        U.set(i(2), j(1), ii(1), jj(2), udat[N + counter][1][2]);
        U.set(i(2), j(2), ii(1), jj(2), udat[N + counter][1][3]);

        U.set(i(1), j(1), ii(2), jj(1), udat[N + counter][2][0]);
        U.set(i(1), j(2), ii(2), jj(1), udat[N + counter][2][1]);
        U.set(i(2), j(1), ii(2), jj(1), udat[N + counter][2][2]);
        U.set(i(2), j(2), ii(2), jj(1), udat[N + counter][2][3]);

        U.set(i(1), j(1), ii(2), jj(2), udat[N + counter][3][0]);
        U.set(i(1), j(2), ii(2), jj(2), udat[N + counter][3][1]);
        U.set(i(2), j(1), ii(2), jj(2), udat[N + counter][3][2]);
        U.set(i(2), j(2), ii(2), jj(2), udat[N + counter][3][3]);

        data[pos] = U * (T0 * T1);
        data[pos].noprime();
        counter++;
      } else {
        data[pos] = T0 * T1;
      }
    }
    new_phys_index[pos] = j;
  }

  for (int pos = 1; pos <= (N - 1); pos++) {
    auto T0 = data[pos];
    auto T1 = data[pos + 1];

    std::vector<Index> common_inds;
    for (auto& I : T0.inds()) {
      if (hasindex(T1, I)) {
        common_inds.push_back(I);
      }
    }
    auto C = combiner(common_inds);
    data[pos] = C * data[pos];
    data[pos + 1] = C * data[pos + 1];
  }

  for (int pos = 1; pos <= N; pos++) {
    for (auto& I : data[pos].inds()) {
      if ((I.type() == Site) && (I != new_phys_index[pos])) {
        index_to_meas[pos] = I;
      }
    }
  }

  *this = this->measure(index_to_meas);

  return *this;
}

ClusterMap& ClusterMap::cz_brick_channel() {
  iteration_num = iteration_num + 1;
  int N = this->get_size();
  ClusterMap new_cluster("it:" + std::to_string(iteration_num) + " ", N);

  auto psi = new_cluster.to_mps();

  for (int pos = 1; pos <= N - 1; pos++) {
    psi.position(pos);
    ITensor T0 = psi.A(pos);
    ITensor T1 = psi.A(pos + 1);
    Index i, j;

    for (auto& I : T0.inds()) {
      if (I.type() == Site) {
        i = I;
      }
    }
    auto ii = prime(i);

    for (auto& I : T1.inds()) {
      if (I.type() == Site) j = I;
    }
    auto jj = prime(j);
    auto U = ITensor(i, j, ii, jj);
      
    U.set(i(1), j(1), ii(1), jj(1), 1);
    U.set(i(1), j(2), ii(1), jj(1), 0);
    U.set(i(2), j(1), ii(1), jj(1), 0);
    U.set(i(2), j(2), ii(1), jj(1), 0);

    U.set(i(1), j(1), ii(1), jj(2), 0);
    U.set(i(1), j(2), ii(1), jj(2), 1);
    U.set(i(2), j(1), ii(1), jj(2), 0);
    U.set(i(2), j(2), ii(1), jj(2), 0);

    U.set(i(1), j(1), ii(2), jj(1), 0);
    U.set(i(1), j(2), ii(2), jj(1), 0);
    U.set(i(2), j(1), ii(2), jj(1), 1);
    U.set(i(2), j(2), ii(2), jj(1), 0);

    U.set(i(1), j(1), ii(2), jj(2), 0);
    U.set(i(1), j(2), ii(2), jj(2), 0);
    U.set(i(2), j(1), ii(2), jj(2), 0);
    U.set(i(2), j(2), ii(2), jj(2), 1);

    auto wf = T0 * T1;
    wf *= U;
    wf.noprime();

    ITensor S, V;
    ITensor W = psi.A(pos);
    svd(wf, W, S, V, {"Cutoff=", 1E-18});
    psi.setA(pos, W);
    psi.setA(pos + 1, S * V);
  }
  psi.orthogonalize();
  normalize(psi);

  new_cluster = ClusterMap(psi, new_cluster.get_name(), iteration_num);

  std::map<int, Index> index_to_meas;
  std::map<int, Index> new_phys_index;
  int counter = 0;

  for (int pos = 1; pos <= N; pos++) {
    ITensor T0 = data[pos];
    ITensor T1 = new_cluster.get_tensor(pos);
    Index i, j;
    for (auto& I : T0.inds()) {
      if (I.type() == Site) {
        i = I;
      }
    }
    for (auto& I : T1.inds()) {
      if (I.type() == Site) {
        j = I;
      }
    }


    if (iteration_num % 2 == 1) {
      if (((pos - 1) % 8 == 1) || ((pos - 1) % 8 == 3)) {

        auto ii = prime(i);
        auto jj = prime(j);
        auto U = ITensor(i, j, ii, jj);

        U.set(i(1), j(1), ii(1), jj(1), 1.0);
        U.set(i(1), j(2), ii(1), jj(1), 0.0);
        U.set(i(2), j(1), ii(1), jj(1), 0.0);
        U.set(i(2), j(2), ii(1), jj(1), 0.0);

        U.set(i(1), j(1), ii(1), jj(2), 0.0);
        U.set(i(1), j(2), ii(1), jj(2), 1.0);
        U.set(i(2), j(1), ii(1), jj(2), 0.0);
        U.set(i(2), j(2), ii(1), jj(2), 0.0);

        U.set(i(1), j(1), ii(2), jj(1), 0.0);
        U.set(i(1), j(2), ii(2), jj(1), 0.0);
        U.set(i(2), j(1), ii(2), jj(1), 1.0);
        U.set(i(2), j(2), ii(2), jj(1), 0.0);

        U.set(i(1), j(1), ii(2), jj(2), 0.0);
        U.set(i(1), j(2), ii(2), jj(2), 0.0);
        U.set(i(2), j(1), ii(2), jj(2), 0.0);
        U.set(i(2), j(2), ii(2), jj(2), -1.0);

        data[pos] = U * (T0 * T1);
        data[pos].noprime();
        counter++;

      } else {
        data[pos] = T0 * T1;
      }
    } else {
      if (((pos - 1) % 8 == 5) || ((pos - 1) % 8 == 7)) {

        auto ii = prime(i);
        auto jj = prime(j);
        auto U = ITensor(i, j, ii, jj);

        U.set(i(1), j(1), ii(1), jj(1), 1.0);
        U.set(i(1), j(2), ii(1), jj(1), 0.0);
        U.set(i(2), j(1), ii(1), jj(1), 0.0);
        U.set(i(2), j(2), ii(1), jj(1), 0.0);

        U.set(i(1), j(1), ii(1), jj(2), 0.0);
        U.set(i(1), j(2), ii(1), jj(2), 1.0);
        U.set(i(2), j(1), ii(1), jj(2), 0.0);
        U.set(i(2), j(2), ii(1), jj(2), 0.0);

        U.set(i(1), j(1), ii(2), jj(1), 0.0);
        U.set(i(1), j(2), ii(2), jj(1), 0.0);
        U.set(i(2), j(1), ii(2), jj(1), 1.0);
        U.set(i(2), j(2), ii(2), jj(1), 0.0);

        U.set(i(1), j(1), ii(2), jj(2), 0.0);
        U.set(i(1), j(2), ii(2), jj(2), 0.0);
        U.set(i(2), j(1), ii(2), jj(2), 0.0);
        U.set(i(2), j(2), ii(2), jj(2), -1.0);

        data[pos] = U * (T0 * T1);

        data[pos].noprime();
        counter++;
      } else {
        data[pos] = T0 * T1;
      }
    }
    new_phys_index[pos] = j;
  }

  for (int pos = 1; pos <= (N - 1); pos++) {
    auto T0 = data[pos];
    auto T1 = data[pos + 1];

    std::vector<Index> common_inds;
    for (auto& I : T0.inds()) {
      if (hasindex(T1, I)) {
        common_inds.push_back(I);
      }
    }
    auto C = combiner(common_inds);
    data[pos] = C * data[pos];
    data[pos + 1] = C * data[pos + 1];
  }

  for (int pos = 1; pos <= N; pos++) {
    for (auto& I : data[pos].inds()) {
      if ((I.type() == Site) && (I != new_phys_index[pos])) {
        index_to_meas[pos] = I;
      }
    }
  }
  for (int pos = 1; pos <= N; pos++) {
    *this = this->apply_brick_U(pos, index_to_meas[pos]);
  }
  for (int pos = 1; pos <= N; pos++) {
    for (auto& I : data[pos].inds()) {
      if ((I.type() == Site) && (I != new_phys_index[pos])) {
        index_to_meas[pos] = I;
      }
    }
  }

  *this = this->measure(index_to_meas);

  return *this;
}

ClusterMap& ClusterMap::apply_brick_U(int pos, Index i) {
  int N = this->get_size();

  double glob, psi, phi, chi;
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<double> distribution1(0.0, 6.28318530718);
  std::uniform_real_distribution<double> distribution2(0.5, 0.5);
  std::uniform_real_distribution<double> sample(0.0, 1.0);
  glob = 0;
  chi = 0;
  phi =  3.14159265359 / 4.0;

  auto X = sample(generator);

  if (X<=0.5){
    psi =  3.14159265359 / 2.0;
  } else {
    psi =  -3.14159265359 / 2.0;
  }

  auto newI = Index(this->get_name() + " " + std::to_string(pos), 2, Site);

  std::complex<double> im(0.0, 1.0);
  auto U = ITensor(i, newI);
  U.set(i(1), newI(1), 1.0/sqrt(2));

  auto val = std::exp(im * (psi))/sqrt(2.0);
  U.set(i(2), newI(1), val);
  U.set(i(1), newI(2), 1.0/sqrt(2.0));
  U.set(i(2), newI(2), -val);

  data[pos] = data[pos] * U;

  return *this;
}

