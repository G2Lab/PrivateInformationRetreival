#pragma once

namespace constants
{
    // DEBUG PARAMETERS

    const int DEBUG = 0;
    
    // BGV PARAMETERS

    // Plaintext prime modulus
    const unsigned long P = 131;
    // Cyclotomic polynomial - defines phi(m)
    const unsigned long M = 17293;
    // Hensel lifting (default = 1)
    const unsigned long R = 1;
    // Number of bits of the modulus chain
    const unsigned long BITS = 431;
    // Number of columns of Key-Switching matrix (default = 2 or 3)
    const unsigned long C = 3;

}
