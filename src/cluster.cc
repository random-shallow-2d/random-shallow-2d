#include "clustermap.h"
#include "itensor/all.h"

using namespace itensor;

int main(int argc, char *argv[]) {
  int method = strtol(argv[1],NULL,10);
  float err = atof(argv[2]);
  int maxD = strtol(argv[3],NULL,10);
  int N = strtol(argv[4],NULL,10);
  int M = strtol(argv[5],NULL,10);
  std::string out_filename = argv[6];

  float x_lb = 0.0;
  float x_ub = 1.0;

  std::ofstream ofile;
  ofile.open(out_filename);
  if (method == 2) {
    ofile << "Method: truncation with error " << err << std::endl;
  } else if (method == 1) {
    ofile << "Method: truncation with max bond " << maxD << std::endl;
  }

  ofile << "Iterations:" << M << std::endl;

  std::string name = "Test";

  ClusterMap cluster(name, N);

  for (int m = 1; m <= M; m++) {
    ofile << "***************************************************************"
        << std::endl;
    ofile << "***************************************************************"
        << std::endl;
    ofile << "***************************************************************"
        << std::endl;
    ofile << "It:" << m << std::endl;
    std::cout << "Current It: " << m << " ";
    cluster.sphere_channel(x_lb, x_ub);
    auto psi = cluster.to_mps();
    auto psi2 = psi;

    ofile << "Before ortho max:" << maxM(psi) << std::endl;
    psi2.orthogonalize();
    if (method == 2) {
      psi.orthogonalize({"Cutoff", err});
    }
    if (method == 1) {
      psi.orthogonalize({"Maxm", maxD});
    }
    normalize(psi);
    normalize(psi2);
    ofile << "After ortho max:" << maxM(psi) << std::endl;
    cluster = ClusterMap(psi, name, m);

    EntanglementData edata = cluster.get_entanglement();
    ofile << "Entanglement:" << std::endl;
    for (int j = 1; j<= N-1 ;j++) {
        ofile << edata.s[j] << ",";
    }
    ofile << std::endl;

    auto svs = edata.svals;
    ofile << "Singular Values at N/2:" << std::endl;
    for (int i = 0; i < svs[N/2].size(); i++) {
        ofile << svs[N/2][i] << ",";
    }
    ofile << std::endl;

    auto inner = overlap(psi, psi2);
    ofile << "Inner Product between iterations:" << inner << std::endl;
  }
  ofile.close();
  return 0;
}
