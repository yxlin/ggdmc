#include <iostream>
#include <vector>

// Simulate the update_profile_probabilities logic
void test_profile_sum(const std::vector<double>& pi_compact, size_t n_profile) {
    std::vector<double> complete(n_profile);
    
    double sum = 0.0;
    for (size_t i = 0; i < n_profile - 1; ++i) {
        complete[i] = pi_compact[i];
        sum += pi_compact[i];
    }
    complete[n_profile - 1] = 1.0 - sum;
    
    // Print results
    std::cout << "Input compact probs: ";
    for (auto p : pi_compact) std::cout << p << " ";
    std::cout << "\n";
    
    std::cout << "Complete probs: ";
    for (auto p : complete) std::cout << p << " ";
    std::cout << "\n";
    
    std::cout << "Sum: " << sum << ", Last prob: " << complete[n_profile-1] << "\n";
    
    // Check for issues
    if (complete[n_profile-1] < 0.0) {
        std::cout << "ERROR: Last probability is negative!\n";
    }
    if (complete[n_profile-1] > 1.0) {
        std::cout << "ERROR: Last probability > 1!\n";
    }
    std::cout << "\n";
}

int main() {
    std::cout << "Test 1: Valid probabilities\n";
    test_profile_sum({0.2, 0.3, 0.1}, 4);
    
    std::cout << "Test 2: Probabilities sum > 1\n";
    test_profile_sum({0.4, 0.5, 0.3}, 4);
    
    std::cout << "Test 3: Probabilities sum to exactly 1\n";
    test_profile_sum({0.333, 0.333, 0.334}, 4);
    
    std::cout << "Test 4: Near-boundary case\n";
    test_profile_sum({0.01, 0.01, 0.97}, 4);
    
    return 0;
}
