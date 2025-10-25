#include <iostream>
#include <algorithm>

// Simulate DINA probability computation
double compute_dina_prob(double guess, double slip, bool has_all_skills, double eps = 1e-12) {
    double eta = has_all_skills ? 1.0 : 0.0;
    double p = eta * (1.0 - slip) + (1.0 - eta) * guess;
    
    // Clamp to [eps, 1-eps]
    p = std::max(eps, std::min(1.0 - eps, p));
    
    return p;
}

void test_params(double guess, double slip) {
    std::cout << "guess=" << guess << ", slip=" << slip << "\n";
    std::cout << "  P(Y=1|no skills) = " << compute_dina_prob(guess, slip, false, 1e-12) << "\n";
    std::cout << "  P(Y=1|has skills) = " << compute_dina_prob(guess, slip, true, 1e-12) << "\n";
    
    // Check if parameters make sense
    if (guess < 0 || guess > 1) std::cout << "  WARNING: guess outside [0,1]\n";
    if (slip < 0 || slip > 1) std::cout << "  WARNING: slip outside [0,1]\n";
    
    double p0 = compute_dina_prob(guess, slip, false, 1e-12);
    double p1 = compute_dina_prob(guess, slip, true, 1e-12);
    if (p0 > p1) {
        std::cout << "  ERROR: P(Y=1|no skills) > P(Y=1|has skills) - model is non-monotonic!\n";
    }
    std::cout << "\n";
}

int main() {
    std::cout << "Test 1: Valid parameters\n";
    test_params(0.2, 0.1);
    
    std::cout << "Test 2: guess > 1\n";
    test_params(1.5, 0.1);
    
    std::cout << "Test 3: slip > 1\n";
    test_params(0.2, 1.5);
    
    std::cout << "Test 4: Both negative\n";
    test_params(-0.3, -0.2);
    
    std::cout << "Test 5: guess + slip > 1 (non-monotonic!)\n";
    test_params(0.8, 0.3);
    
    return 0;
}
