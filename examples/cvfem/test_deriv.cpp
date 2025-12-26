#include <iostream>
#include <iomanip>
#include "external/master_element/Hex8CVFEM.h"

int main() {
    accel::HexSCS meSCS;
    const double* intgLocShift = meSCS.integration_location_shift();
    
    std::cout << std::scientific << std::setprecision(10);
    
    for (int ip = 0; ip < 12; ++ip) {
        double s1 = intgLocShift[ip*3 + 0];
        double s2 = intgLocShift[ip*3 + 1];
        double s3 = intgLocShift[ip*3 + 2];
        
        std::cout << "// ip " << ip << ": s1=" << s1 << ", s2=" << s2 << ", s3=" << s3 << std::endl;
        std::cout << "{{";
        
        // Compute derivatives manually
        double half = 0.5;
        double one4th = 0.25;
        double s1s2 = s1 * s2;
        double s2s3 = s2 * s3;
        double s1s3 = s1 * s3;
        
        double deriv[8][3];
        
        // d/dksi
        deriv[0][0] = half * (s3 + s2) - s2s3 - one4th;
        deriv[1][0] = half * (-s3 - s2) + s2s3 + one4th;
        deriv[2][0] = half * (-s3 + s2) - s2s3 + one4th;
        deriv[3][0] = half * (+s3 - s2) + s2s3 - one4th;
        deriv[4][0] = half * (-s3 + s2) + s2s3 - one4th;
        deriv[5][0] = half * (+s3 - s2) - s2s3 + one4th;
        deriv[6][0] = half * (+s3 + s2) + s2s3 + one4th;
        deriv[7][0] = half * (-s3 - s2) - s2s3 - one4th;
        
        // d/deta
        deriv[0][1] = half * (s3 + s1) - s1s3 - one4th;
        deriv[1][1] = half * (s3 - s1) + s1s3 - one4th;
        deriv[2][1] = half * (-s3 + s1) - s1s3 + one4th;
        deriv[3][1] = half * (-s3 - s1) + s1s3 + one4th;
        deriv[4][1] = half * (-s3 + s1) + s1s3 - one4th;
        deriv[5][1] = half * (-s3 - s1) - s1s3 - one4th;
        deriv[6][1] = half * (s3 + s1) + s1s3 + one4th;
        deriv[7][1] = half * (s3 - s1) - s1s3 + one4th;
        
        // d/dzeta
        deriv[0][2] = half * (s2 + s1) - s1s2 - one4th;
        deriv[1][2] = half * (s2 - s1) + s1s2 - one4th;
        deriv[2][2] = half * (-s2 - s1) - s1s2 - one4th;
        deriv[3][2] = half * (-s2 + s1) + s1s2 - one4th;
        deriv[4][2] = half * (-s2 - s1) + s1s2 + one4th;
        deriv[5][2] = half * (-s2 + s1) - s1s2 + one4th;
        deriv[6][2] = half * (s2 + s1) + s1s2 + one4th;
        deriv[7][2] = half * (s2 - s1) - s1s2 + one4th;
        
        for (int n = 0; n < 8; ++n) {
            std::cout << "{" << deriv[n][0] << ", " << deriv[n][1] << ", " << deriv[n][2] << "}";
            if (n < 7) std::cout << ", ";
        }
        std::cout << "}";
        if (ip < 11) std::cout << ",";
        std::cout << std::endl;
    }
    
    return 0;
}
