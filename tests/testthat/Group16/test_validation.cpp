#include <iostream>
#include <vector>
#include <cmath>

// Simulate the validation logic from cdm.h lines 1381-1391
bool validate_simplex(const std::vector<double>& w) {
    if (w.empty()) return false;
    
    // Check finite
    for (auto x : w) {
        if (!std::isfinite(x)) {
            std::cout << "  Failed: not finite\n";
            return false;
        }
    }
    
    // Check bounds [-1e-6, 1+1e-6]
    for (auto x : w) {
        if (x < -1e-6 || x > 1.0 + 1e-6) {
            std::cout << "  Failed: value " << x << " outside bounds\n";
            return false;
        }
    }
    
    // Check sum
    double sum = 0;
    for (auto x : w) sum += x;
    if (std::abs(sum - 1.0) > 1e-6) {
        std::cout << "  Failed: sum = " << sum << ", not close to 1\n";
        return false;
    }
    
    return true;
}

int main() {
    std::cout << "Test 1: Valid simplex\n";
    std::vector<double> v1 = {0.2, 0.3, 0.1, 0.4};
    std::cout << "  Result: " << (validate_simplex(v1) ? "PASS" : "FAIL") << "\n\n";
    
    std::cout << "Test 2: Negative probability\n";
    std::vector<double> v2 = {0.4, 0.5, 0.3, -0.2};
    std::cout << "  Result: " << (validate_simplex(v2) ? "PASS" : "FAIL") << "\n\n";
    
    std::cout << "Test 3: Sum > 1\n";
    std::vector<double> v3 = {0.4, 0.4, 0.4, 0.4};
    std::cout << "  Result: " << (validate_simplex(v3) ? "PASS" : "FAIL") << "\n\n";
    
    std::cout << "Test 4: Small negative (within tolerance)\n";
    std::vector<double> v4 = {0.3, 0.3, 0.4, -0.000001};
    std::cout << "  Result: " << (validate_simplex(v4) ? "PASS" : "FAIL") << "\n\n";
    
    std::cout << "Test 5: Zero probability at boundary\n";
    std::vector<double> v5 = {0.333, 0.333, 0.334, 0.0};
    std::cout << "  Result: " << (validate_simplex(v5) ? "PASS" : "FAIL") << "\n\n";
    
    return 0;
}
