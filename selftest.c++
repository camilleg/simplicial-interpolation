// Exercise the R^d to R^e mapping,
// for various d, e, and number of point-pairs n.

#include "si.h"

bool test(int d, int e, int n, qi_kind k) {
  if (!init(d, e, n, k))
    return false;

  vector<vertex> qtest;
  randomSites(qtest, d, 5, 1e6);
  for (const auto& q: qtest) {
    dump_v("   query ", q);
    const auto p = eval(q);
    dump_v("  result ", p);
  }
  return true;
}

int main(int argc, char** argv) {
  if (argc != 1) {
    printf("usage: %s\n", argv[0]);
    return 1;
  }
  for (auto d=5; d>=2; --d) {
  for (auto e=d+1; e<d+10; e+=3) {
  for (auto n=d+10; n<d+35; n+=10) {
    cout << "d=" << d << ", e=" << e << ", n=" << n << std::endl;
    for (auto k: {qi_kind::random, qi_kind::geneticAlgorithm})
      if (!test(d, e, n, k))
	return -1;

    terminate();
  }}}
  return 0;
}
