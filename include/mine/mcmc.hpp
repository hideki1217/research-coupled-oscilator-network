#include <algorithm>
#include <memory>
#include <random>

namespace research {

template <typename State>
struct Metropolice_ {
 public:
  using rng_t = std::mt19937;
  using state_t = State;
  using hamiltonian_t = std::function<double(const state_t&, rng_t&)>;

 public:
  enum class Result {
    Accepted,
    Rejected,
  };

 public:
  const double beta;

 private:
  std::unique_ptr<state_t> _state;
  double _E;
  rng_t rng;
  hamiltonian_t H;

 public:
  Metropolice_(hamiltonian_t H, const state_t& initial, double beta, int seed)
      : rng(seed),
        _state(std::make_unique<state_t>(initial)),
        beta(beta),
        H(H) {
    _E = H(*_state, rng);
  }

  const state_t& state() { return _state; }
  double energy() { return _E; }

  void swap(Metropolice_<state_t>& other) {
    assert(H == other.H);
    _state.swap(other._state);
    std::swap(_E, other._E);
  }

  Result update() {
    auto rebert = fluctuate(*_state, rng);
    const auto E = H(*_state, rng);

    const auto p = std::exp(-beta * (E - _E));
    if (1 < p || random_p(rng) < p) {
      _E = E;

      return Result::Accepted;
    } else {
      rebert();

      return Result::Rejected;
    }
  }

 private:
  double random_p(rng_t& rng) {
    std::uniform_real_distribution unif(0., 1.);
    return unif(rng);
  }

  auto fluctuate(state_t& state, rng_t& rng) {
    std::normal_distribution normal(0., 1.);
    std::uniform_int_distribution indexer(0, int(state.size() - 1));

    const auto index = indexer(rng);
    const auto noise = normal(rng);

    const auto prev = state[index];
    state[index] += noise;
    auto rebert = [=]() mutable { state[index] = prev; };
    return rebert;
  }
};

template <typename State>
Metropolice_<State> Metropolice(typename Metropolice_<State>::hamiltonian_t H,
                                const State& initial, double beta,
                                int seed = 20) {
  return Metropolice_(H, initial, beta);
}
}  // namespace research