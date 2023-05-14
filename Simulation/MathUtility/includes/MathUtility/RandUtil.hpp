/* Random number utility functions */

#ifndef RANDUTIL_HPP
#define RANDUTIL_HPP

#include <random>
#include <chrono>

class RandUtil
{
private:
  ulong mSeed { constants::SEED };              //!< PRNG seed
  std::mt19937 mGenerator { mSeed };            //!< PRNG generator

public:
  RandUtil()
  {}
  RandUtil( ulong _seed ) : mSeed{ _seed }, mGenerator{ _seed }
  {}

  ulong getSeed() { return mSeed; }
  void setSeed(ulong _newSeed) { mSeed=_newSeed; mGenerator.seed(mSeed); }
  void setRandomSeed() {
    setSeed(
      std::chrono::high_resolution_clock::now().time_since_epoch().count()
    );
  }

  double getUniformRand(double a=0, double b=1)
  {
    std::uniform_real_distribution<double> uniform(a,b);
    return uniform(mGenerator);
  }

  double getNormalRand(double mu=0, double sig=1)
  {
    std::normal_distribution<double> normal(mu,sig);
    return normal(mGenerator);
  }

  double getGammaRand(double alpha, double beta)
  {
    std::gamma_distribution<double> gamma(alpha,beta);
    return gamma(mGenerator);
  }

  int getBinomialRand(int trials, double prob_success)
  {
    std::binomial_distribution<int> binomial(trials,prob_success);
    return binomial(mGenerator);
  }
};

extern RandUtil gen_rand;

#endif
