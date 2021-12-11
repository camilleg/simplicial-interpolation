// Exercise the R^d to R^e mapping,
// for various d, e, and number of point-pairs n.

#include "si.h"

int main(int argc, char** argv) {
  if (argc != 1) {
    printf("usage: %s\n", argv[0]);
    return 1;
  }
  for (auto d=5; d>=2; --d) {
  for (auto e=d+1; e<d+10; e+=3) {
  for (auto n=d+10; n<d+35; n+=10) {
    cout << "d=" << d << ", e=" << e << ", n=" << n << std::endl;
    if (!init(d, e, n, qi_kind::geneticAlgorithm))
      return -1;

    std::vector<vertex> qtest;
    randomSites(qtest, d, 5);
    for (const auto& q: qtest) {
      dump_v("   query ", q);
      const auto p = eval(q);
      dump_v("  result ", p);
    }

    terminate();
  }}}
  return 0;
}
