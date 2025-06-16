#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>

std::vector<std::complex<double>> fft(const std::vector<std::complex<double>> &a) {
    int n = a.size();
    if(n == 1) return a;
    std::vector<std::complex<double>> even(n/2), odd(n/2);
    for (int i = 0; i < n/2; i++) {
        even[i] = a[2*i];
        odd[i] = a[2*i+1];
    }
    std::vector<std::complex<double>> Feven = fft(even);
    std::vector<std::complex<double>> Fodd = fft(odd);
    std::vector<std::complex<double>> y(n);
    for (int k = 0; k < n/2; k++) {
        std::complex<double> t = std::polar(1.0, -2 * M_PI * k / n) * Fodd[k];
        y[k] = Feven[k] + t;
        y[k+n/2] = Feven[k] - t;
    }
    return y;
}

double simpson_integration(const std::function<double(double)> &f, double a, double b, int n) {
    double h = (b - a) / n;
    double s = f(a) + f(b);
    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        s += (i % 2 == 0 ? 2.0 : 4.0) * f(x);
    }
    return s * h / 3.0;
}

class SemiCompleteOptionPricingLibrary {
    double _r, _S0, _K, _T;
public:
    SemiCompleteOptionPricingLibrary(double r, double S0, double K, double T) : _r(r), _S0(S0), _K(K), _T(T) {}
    std::complex<double> heston_characteristic_function(std::complex<double> u, double kappa_v, double theta_v, double sigma_v, double rho, double v0) {
        std::complex<double> I(0, 1);
        double c1 = kappa_v * theta_v;
        std::complex<double> c2 = -std::sqrt(std::pow(rho * sigma_v * u * I - kappa_v, 2) - sigma_v * sigma_v * (-u * I - u * u));
        std::complex<double> c3 = (kappa_v - rho * sigma_v * u * I + c2) / (kappa_v - rho * sigma_v * u * I - c2);
        std::complex<double> H1 = _r * u * I * _T + (c1 / (sigma_v * sigma_v)) * ((kappa_v - rho * sigma_v * u * I + c2) * _T - 2.0 * std::log((1.0 - c3 * std::exp(c2 * _T)) / (1.0 - c3)));
        std::complex<double> H2 = ((kappa_v - rho * sigma_v * u * I + c2) / (sigma_v * sigma_v)) * ((1.0 - std::exp(c2 * _T)) / (1.0 - c3 * std::exp(c2 * _T)));
        return std::exp(H1 + H2 * v0);
    }
    double heston_integration_function(double u, double kappa_v, double theta_v, double sigma_v, double rho, double v0) {
        std::complex<double> I(0, 1);
        std::complex<double> char_func = heston_characteristic_function(u - 0.5 * I, kappa_v, theta_v, sigma_v, rho, v0);
        std::complex<double> integrand = std::exp(I * u * std::log(_S0 / _K)) * char_func;
        return (1.0 / (u * u + 0.25)) * integrand.real();
    }
    double heston_call_value(double kappa_v, double theta_v, double sigma_v, double rho, double v0) {
        auto f = [=](double u) { return heston_integration_function(u, kappa_v, theta_v, sigma_v, rho, v0); };
        double int_value = simpson_integration(f, 0.0, 200.0, 10000);
        double call_value = std::max(0.0, _S0 - std::exp(-_r * _T) * std::sqrt(_S0 * _K) / M_PI * int_value);
        return call_value;
    }
    std::complex<double> merton_characteristic_function_jump(std::complex<double> u, double lamb, double mu, double delta) {
        std::complex<double> I(0, 1);
        double omega = -lamb * (std::exp(mu + 0.5 * delta * delta) - 1.0);
        return std::exp((I * u * omega + lamb * (std::exp(I * u * mu - 0.5 * u * u * delta * delta) - 1.0)) * _T);
    }
    double merton_integration_function(double u, double dummy, double lamb, double mu, double delta) {
        std::complex<double> I(0, 1);
        std::complex<double> char_func = merton_characteristic_function_jump(u - 0.5 * I, lamb, mu, delta);
        std::complex<double> integrand = std::exp(I * u * std::log(_S0 / _K)) * char_func;
        return (1.0 / (u * u + 0.25)) * integrand.real();
    }
    double merton_call_value(double dummy, double lamb, double mu, double delta) {
        auto f = [=](double u) { return merton_integration_function(u, dummy, lamb, mu, delta); };
        double int_value = simpson_integration(f, 0.0, 50.0, 10000);
        double call_value = std::max(0.0, _S0 - std::exp(-_r * _T) * std::sqrt(_S0 * _K) / M_PI * int_value);
        return call_value;
    }
    double calibration_error_function(const std::vector<double> &p0, const std::vector<double> &options, bool isMerton) {
        double sum = 0.0;
        for (double market_price : options) {
            double model_price = 0.0;
            if (isMerton) model_price = merton_call_value(p0[0], p0[1], p0[2], p0[3]);
            else model_price = heston_call_value(p0[0], p0[1], p0[2], p0[3], p0[4]);
            double diff = model_price - market_price;
            sum += diff * diff;
        }
        double rmse = std::sqrt(sum / options.size());
        return rmse;
    }
    std::vector<double> calibrate_merton(const std::vector<double> &options) {
        std::vector<double> best(4, 0.0);
        double best_err = 1e12;
        for (double p0 = 0.075; p0 <= 0.2001; p0 += 0.025) {
            for (double p1 = 0.10; p1 <= 0.4001; p1 += 0.1) {
                for (double p2 = -0.5; p2 <= 0.001; p2 += 0.1) {
                    for (double p3 = 0.10; p3 <= 0.3001; p3 += 0.1) {
                        std::vector<double> params = {p0, p1, p2, p3};
                        double err = calibration_error_function(params, options, true);
                        if (err < best_err) { best_err = err; best = params; }
                    }
                }
            }
        }
        for (int iter = 0; iter < 100; iter++) {
            std::vector<double> new_best = best;
            for (int i = 0; i < 4; i++) {
                double step = 0.001;
                for (int d = -1; d <= 1; d++) {
                    std::vector<double> temp = best;
                    temp[i] += d * step;
                    double err = calibration_error_function(temp, options, true);
                    if (err < best_err) { best_err = err; new_best = temp; }
                }
            }
            best = new_best;
        }
        return best;
    }
    std::complex<double> bates_characteristic_function(std::complex<double> u, double kappa_v, double theta_v, double sigma_v, double rho, double v0, double lamb, double mu, double delta) {
        return heston_characteristic_function(u, kappa_v, theta_v, sigma_v, rho, v0) * merton_characteristic_function_jump(u, lamb, mu, delta);
    }
    double bates_integration_function(double u, double kappa_v, double theta_v, double sigma_v, double rho, double v0, double lamb, double mu, double delta) {
        std::complex<double> I(0, 1);
        std::complex<double> char_func_value = bates_characteristic_function(u - 0.5 * I, kappa_v, theta_v, sigma_v, rho, v0, lamb, mu, delta);
        std::complex<double> integrand = std::exp(I * u * std::log(_S0 / _K)) * char_func_value;
        return (1.0 / (u * u + 0.25)) * integrand.real();
    }
    double bates_call_value(double kappa_v, double theta_v, double sigma_v, double rho, double v0, double lamb, double mu, double delta) {
        auto f = [=](double u) { return bates_integration_function(u, kappa_v, theta_v, sigma_v, rho, v0, lamb, mu, delta); };
        double int_value = simpson_integration(f, 0.0, 200.0, 10000);
        double call_value = std::max(0.0, _S0 - std::exp(-_r * _T) * std::sqrt(_S0 * _K) / M_PI * int_value);
        return call_value;
    }
    double bates_call_fft(double kappa_v, double theta_v, double sigma_v, double rho, double v0, double lamb, double mu, double delta) {
        double k = std::log(_K / _S0);
        double g = 1.0;
        int N = static_cast<int>(g * 4096);
        double eps = 1.0 / (g * 150.0);
        double eta = 2 * M_PI / (N * eps);
        double b_val = 0.5 * N * eps - k;
        std::vector<double> u(N), vo(N);
        for (int i = 0; i < N; i++) {
            u[i] = i + 1;
            vo[i] = eta * i;
        }
        std::vector<std::complex<double>> modcharFunc;
        std::vector<std::complex<double>> modcharFunc1, modcharFunc2;
        double alpha;
        std::complex<double> I(0, 1);
        if (_S0 >= 0.95 * _K) {
            alpha = 1.5;
            modcharFunc.resize(N);
            for (int i = 0; i < N; i++) {
                std::complex<double> v_arr = vo[i] - (alpha + 1.0) * I;
                std::complex<double> denom = (alpha * alpha + alpha - vo[i] * vo[i] + I * (2 * alpha + 1) * vo[i]);
                modcharFunc[i] = std::exp(-_r * _T) * (bates_characteristic_function(v_arr, kappa_v, theta_v, sigma_v, rho, v0, lamb, mu, delta) / denom);
            }
        } else {
            alpha = 1.1;
            modcharFunc1.resize(N);
            modcharFunc2.resize(N);
            for (int i = 0; i < N; i++) {
                std::complex<double> v_arr1 = (vo[i] - I * alpha) - I;
                std::complex<double> denom1 = ((vo[i] - I * alpha) * (vo[i] - I * alpha) - I * (vo[i] - I * alpha));
                modcharFunc1[i] = std::exp(-_r * _T) * (1.0 / (1.0 + I * (vo[i] - I * alpha)) - std::exp(_r * _T) / (I * (vo[i] - I * alpha)) - bates_characteristic_function(v_arr1, kappa_v, theta_v, sigma_v, rho, v0, lamb, mu, delta) / denom1);
                std::complex<double> v_arr2 = (vo[i] + I * alpha) - I;
                std::complex<double> denom2 = ((vo[i] + I * alpha) * (vo[i] + I * alpha) - I * (vo[i] + I * alpha));
                modcharFunc2[i] = std::exp(-_r * _T) * (1.0 / (1.0 + I * (vo[i] + I * alpha)) - std::exp(_r * _T) / (I * (vo[i] + I * alpha)) - bates_characteristic_function(v_arr2, kappa_v, theta_v, sigma_v, rho, v0, lamb, mu, delta) / denom2);
            }
        }
        std::vector<double> delt(N, 0.0);
        delt[0] = 1.0;
        std::vector<double> SimpsonW(N, 0.0);
        for (int i = 0; i < N; i++) {
            int j = i + 1;
            SimpsonW[i] = (3.0 + std::pow(-1.0, j) - delt[i]) / 3.0;
        }
        std::vector<std::complex<double>> FFTFunc(N);
        std::vector<std::complex<double>> payoff;
        if (_S0 >= 0.95 * _K) {
            for (int i = 0; i < N; i++) {
                FFTFunc[i] = std::exp(I * b_val * vo[i]) * modcharFunc[i] * eta * SimpsonW[i];
            }
            payoff = fft(FFTFunc);
            std::vector<double> CallValueM(N);
            for (int i = 0; i < N; i++) {
                CallValueM[i] = std::exp(-alpha * k) / M_PI * payoff[i].real();
            }
            int pos = static_cast<int>((k + b_val) / eps);
            return CallValueM[pos] * _S0;
        } else {
            for (int i = 0; i < N; i++) {
                FFTFunc[i] = std::exp(I * b_val * vo[i]) * (modcharFunc1[i] - modcharFunc2[i]) * 0.5 * eta * SimpsonW[i];
            }
            payoff = fft(FFTFunc);
            std::vector<double> CallValueM(N);
            for (int i = 0; i < N; i++) {
                CallValueM[i] = payoff[i].real() / (std::sinh(alpha * k) * M_PI);
            }
            int pos = static_cast<int>((k + b_val) / eps);
            return CallValueM[pos] * _S0;
        }
    }
};

int main(){
    double S0 = 100.0, K = 100.0, T = 1.0, r = 0.05;
    double kappa_v = 1.5, theta_v = 0.02, sigma_v = 0.15, rho = 0.1, v0 = 0.01;
    double lamb = 0.25, mu = -0.2, delta = 0.1;
    double sigma = std::sqrt(v0);
    SemiCompleteOptionPricingLibrary pricing_lib(r, S0, K, T);
    double bates_price_integration = pricing_lib.bates_call_value(kappa_v, theta_v, sigma_v, rho, v0, lamb, mu, delta);
    std::cout << "Bates Call Price (Integration): " << bates_price_integration << std::endl;
    double bates_price_fft = pricing_lib.bates_call_fft(kappa_v, theta_v, sigma_v, rho, v0, lamb, mu, delta);
    std::cout << "Bates Call Price (FFT): " << bates_price_fft << std::endl;
    return 0;
}
