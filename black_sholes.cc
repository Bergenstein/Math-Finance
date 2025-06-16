#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <stdexcept>
#include <algorithm>

class BlackScholes {
private:
    double S0;
    double K;
    double sigma;
    double T;
    int M;
    int Ite;
    double r;
    double dt;

    double norm_cdf(double x) const {
        return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
    }

public:
    BlackScholes(double S0, double K, double sigma, double T, int M, int Ite, double r)
        : S0(S0), K(K), sigma(sigma), T(T), M(M), Ite(Ite), r(r) 
    {
        if (sigma <= 0)
            throw std::invalid_argument("sigma must be a positive number");
        if (T <= 0)
            throw std::invalid_argument("T must be a positive number");
        if (M <= 0)
            throw std::invalid_argument("M must be a positive integer");
        if (Ite <= 0)
            throw std::invalid_argument("Ite must be a positive integer");

        dt = T / M;
        if (dt <= 0)
            throw std::invalid_argument("dt must be positive");
    }

    std::vector<std::vector<double>> BS_paths() const {
        std::vector<std::vector<double>> S(M + 1, std::vector<double>(Ite, 0.0));
        for (int j = 0; j < Ite; ++j) {
            S[0][j] = S0;
        }
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> dist(0.0, 1.0);
        double sqrt_dt = std::sqrt(dt);
        for (int t = 1; t <= M; ++t) {
            for (int j = 0; j < Ite; ++j) {
                double ran = dist(gen);
                S[t][j] = S[t - 1][j] * std::exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrt_dt * ran);
            }
        }
        return S;
    }

    double bs_option_mc(char optype) const {
        auto paths = BS_paths();
        double sum_payoff = 0.0;
        for (int j = 0; j < Ite; ++j) {
            double payoff = 0.0;
            if (optype == 'C') {
                payoff = std::max(0.0, paths[M][j] - K);
            } else if (optype == 'P') {
                payoff = std::max(0.0, K - paths[M][j]);
            } else {
                throw std::invalid_argument("Invalid option type. Use 'C' for Call or 'P' for Put.");
            }
            sum_payoff += payoff;
        }
        double average = sum_payoff / Ite;
        return std::exp(-r * T) * average;
    }

    double bs_option_cf(char optype) const {
        double d1 = (std::log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
        double d2 = d1 - sigma * std::sqrt(T);
        double price = 0.0;
        if (optype == 'C') {
            price = S0 * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
        } else if (optype == 'P') {
            price = K * std::exp(-r * T) * norm_cdf(-d2) - S0 * norm_cdf(-d1);
        } else {
            throw std::invalid_argument("Invalid option type. Use 'C' for Call or 'P' for Put.");
        }
        return price;
    }

    void print_paths(int n) const {
        if (n <= 0 || n > Ite)
            throw std::invalid_argument("n must be a positive integer not greater than Ite");
        auto paths = BS_paths();
        std::cout << "Simulated Price Paths (first " << n << " iterations):\n";
        for (int t = 0; t <= M; ++t) {
            std::cout << "Step " << t << ": ";
            for (int j = 0; j < n; ++j) {
                std::cout << paths[t][j] << " ";
            }
            std::cout << "\n";
        }
    }
};

int main() {
    double S0 = 80.0;
    double K = 80.0;
    double sigma = 0.35;
    double T = 0.25;
    int M = 100;
    int Ite = 10000;
    double r = 0.055;

    try {
        BlackScholes bs(S0, K, sigma, T, M, Ite, r);
        double call_mc = bs.bs_option_mc('C');
        double call_cf = bs.bs_option_cf('C');
        double put_mc = bs.bs_option_mc('P');
        double put_cf = bs.bs_option_cf('P');
        std::cout << "Call option price (Monte Carlo): " << call_mc << "\n";
        std::cout << "Call option price (Closed Form): " << call_cf << "\n";
        std::cout << "Put option price (Monte Carlo): " << put_mc << "\n";
        std::cout << "Put option price (Closed Form): " << put_cf << "\n";
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
    return 0;
}
