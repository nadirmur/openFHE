#include "openfhe.h"
#include <iostream>
#include <algorithm>
#include <cmath> 
#include <vector>
#include <tuple>
#include <random> 
#include <fstream>
#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <NTL/ZZ.h>


using namespace lbcrypto;
    

// AUXILIARY FUNCTIONS 

// Multiplication between two polynomials with large integer coefficients (type NTL::ZZ)
std::vector<NTL::ZZ> multiplyPolynomials(const std::vector<NTL::ZZ>& p1, const std::vector<NTL::ZZ>& p2) {
    int64_t n1 = p1.size(); 
    int64_t n2 = p2.size(); 

    std::vector<NTL::ZZ> result(n1 + n2 - 1, NTL::ZZ(0));

    for (int64_t i = 0; i < n1; ++i) {
        for (int64_t j = 0; j < n2; ++j) {
            result[i + j] = (result[i + j] + (p1[i] * p2[j]));
        }
    }

    return result;
}

// Reduction of a polynomial with large integer coefficients (type NTL::ZZ) modulo X^n+1 
std::vector<NTL::ZZ> reducePolynomial(const std::vector<NTL::ZZ>& poly, int64_t n) {
    std::vector<NTL::ZZ> reduced(n, NTL::ZZ(0));

    for (int64_t i = 0; i < static_cast<int64_t>(poly.size()); ++i) {
        if (i < n) {
            reduced[i] += poly[i];
        } else {
            int64_t index = i % n;
            reduced[index] -= poly[i];
        }
    }

    return reduced;
}

// Reduction modulo q with centered representatives of large integers (type NTL::ZZ)
NTL::ZZ mod_reduce(NTL::ZZ num, NTL::ZZ q) {
    NTL::ZZ reduced = num % q;
    if (reduced < -q / 2) {
        reduced += q;
    } else if (reduced >= q / 2) {
        reduced -= q;
    }
    return reduced;
}


// Reduction mod q with centered representatives of a vector of large integers (type NTL::ZZ)
std::vector<NTL::ZZ> reduceVectorModulo(const std::vector<NTL::ZZ>& vec, NTL::ZZ q) {
    std::vector<NTL::ZZ> reduced_vec(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        reduced_vec[i] = mod_reduce(vec[i], q);
    }
    return reduced_vec;
}


// Sum of vectors of large integers (type NTL::ZZ)
std::vector<NTL::ZZ> addVectors(const std::vector<NTL::ZZ>& v1, const std::vector<NTL::ZZ>& v2) {
    size_t size = std::max(v1.size(), v2.size());
    std::vector<NTL::ZZ> result(size, NTL::ZZ(0));

    for (size_t i = 0; i < size; ++i) {
        NTL::ZZ val1 = (i < v1.size()) ? v1[i] : NTL::ZZ(0); 
        NTL::ZZ val2 = (i < v2.size()) ? v2[i] : NTL::ZZ(0); 
        result[i] = val1 + val2;
    }

    return result;
}


// Generate a random vector of length "size" of integers between 0 and "max_val" (set to 5 in this case)
std::vector<int64_t> generateRandomVector(size_t size, int64_t max_val = 5) {
    std::vector<int64_t> vec(size);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, max_val);
    for (auto& v : vec) {
        v = dis(gen);
    }
    return vec;
}

// Return true if exponent is a power of two.
bool isPowerOfTwo(uint64_t exponent) {
    return exponent != 0 && (exponent & (exponent - 1)) == 0;
}

// Multiplicative depth required by the square-and-multiply routine below.
uint32_t multiplicativeDepthSquareAndMultiply(uint64_t exponent) {
    if (exponent == 0) {
        throw std::invalid_argument("The exponent e must be at least 1.");
    }

    uint32_t floorLog2 = 0;
    uint64_t temp = exponent;
    while (temp >>= 1) {
        ++floorLog2;
    }

    return isPowerOfTwo(exponent) ? floorLog2 : floorLog2 + 1;
}

// Compute base^exponent with the square-and-multiply algorithm.
Ciphertext<DCRTPoly> squareAndMultiply(
    const CryptoContext<DCRTPoly>& cryptoContext,
    const Ciphertext<DCRTPoly>& base,
    uint64_t exponent) {

    if (exponent == 0) {
        throw std::invalid_argument("The exponent e must be at least 1.");
    }

    Ciphertext<DCRTPoly> result;
    Ciphertext<DCRTPoly> currentPower = base;
    bool resultInitialized = false;

    while (exponent > 0) {
        if ((exponent & 1ULL) != 0) {
            if (!resultInitialized) {
                result = currentPower;
                resultInitialized = true;
            } else {
                result = cryptoContext->EvalMult(result, currentPower);
            }
        }

        exponent >>= 1;
        if (exponent > 0) {
            currentPower = cryptoContext->EvalMult(currentPower, currentPower);
        }
    }

    return result;
}


///////// MAIN PROGRAM ////////////

int main() {
	auto start = std::chrono::high_resolution_clock::now();	
	
	int64_t n = 8192; // Degree of the polynomial X^n+1, i.e., Ring Dimension
	uint64_t samples = 100; // Number of samples
	int aux_size = 12; // Length of vectors used for creating the plaintexts
	uint64_t e = 6; // Exponent in C^e. Change this value to test another exponent.
	uint32_t M = std::max<uint32_t>(1, multiplicativeDepthSquareAndMultiply(e)); // Multiplicative depth required to compute C^e with square-and-multiply
	
	//Set cryptographic context
    CCParams<CryptoContextBGVRNS> parameters;
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetStandardDeviation(3.2);
    parameters.SetRingDim(n);
    parameters.SetPlaintextModulus(65537);
    parameters.SetMultiplicativeDepth(M);
		
	CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
	
    cryptoContext->Enable(PKE);
	cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
		
	// Output file for the first coefficient of the critical value of C^e
	std::string filename = "critical_values_e" + std::to_string(e) + ".txt";
	std::ofstream outfile(filename, std::ios::app);
	if (!outfile.is_open()) {
		std::cerr << "Impossible open file " << filename << " to write." << std::endl;
		return 1;
	}

	//Start loop
    for (size_t i = 0; i < samples; ++i) {

        //Key generation
        KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
		cryptoContext->EvalMultKeyGen(keyPair.secretKey);
        
        
        //Generation of plaintext
		auto rand_vec = generateRandomVector(aux_size);
		Plaintext pt = cryptoContext->MakePackedPlaintext(rand_vec);
		
       		
		// Compute C^e with square-and-multiply
		auto ciphertext = cryptoContext->Encrypt(keyPair.publicKey, pt);
		auto cPow = squareAndMultiply(cryptoContext, ciphertext, e);
		
		
		// Start computation critical value of the final ciphertext
		
		DCRTPoly s = keyPair.secretKey->GetPrivateElement();
		s.SetFormat(Format::COEFFICIENT);
		auto sk = s.CRTInterpolate(); //extraction of the secret key
		
        auto ciphertext00 = cPow->GetElements()[0];
        auto ciphertext11 = cPow->GetElements()[1];
		ciphertext00.SetFormat(Format::COEFFICIENT);
        ciphertext11.SetFormat(Format::COEFFICIENT);
		auto c0 = ciphertext00.CRTInterpolate(); //extraction of the first polynomial c0 from the ciphertext
        auto c1 = ciphertext11.CRTInterpolate(); //extraction of the second polynimial c1 from the ciphertext
		
		// Start transform type of coefficients of the polynomial sk (secret key) c0 and c1 (polynomials of the ciphertext) in NTL::ZZ in order to manage large integers
		std::vector<NTL::ZZ> c0int;
		std::vector<NTL::ZZ> c1int;
		std::vector<NTL::ZZ> skint;
		
		
		for (size_t j = 0; j < c0.GetLength(); ++j) {
			std::string c0_str = c0[j].ToString();  
			c0int.push_back(NTL::conv<NTL::ZZ>(c0_str.c_str()));
			std::string c1_str = c1[j].ToString();  
			c1int.push_back(NTL::conv<NTL::ZZ>(c1_str.c_str()));
			std::string sk_str = sk[j].ToString();  
			skint.push_back(NTL::conv<NTL::ZZ>(sk_str.c_str())); 
		}
		// End transformation
		
		auto qtemp = c0.GetModulus();
		std::string qtemp2 = qtemp.ToString();
		NTL::ZZ q = NTL::conv<NTL::ZZ>(qtemp2.c_str()); //extraction of ciphertext modulo and conversion in type NTL::ZZ
		
		// Computation of the critical value: criticalValue = c0 + c1 sk
		std::vector<NTL::ZZ> c1s = multiplyPolynomials(c1int, skint);
		std::vector<NTL::ZZ> CV = addVectors(c0int,c1s);
		std::vector<NTL::ZZ> CVModPol = reducePolynomial(CV, n);
		std::vector<NTL::ZZ> criticalValue = reduceVectorModulo(CVModPol, q);		
		
		
		// Save the first coefficient of the critical value of C^e.
		outfile << criticalValue[1] << std::endl;
			
				
	}
	//End loop
	
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    
	std::cout << "Exponent e: " << e << std::endl;
	std::cout << "Output file: " << filename << std::endl;
    std::cout << "Tempo di esecuzione: " << duration.count() << " secondi" << std::endl;

    return 0;
}

