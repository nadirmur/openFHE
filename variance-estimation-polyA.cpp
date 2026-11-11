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
#include <sstream>
#include <NTL/ZZ.h>


using namespace lbcrypto;


////////////////////////////////////////////////////
//                                                //
//   GIVEN A CIPHERTEXT C, THIS PROGRAM EVALUATE  //
//   THE VARIANCE OF THE FIRST COEFFICIENT OF THE //
//	 CRITICAL VALUE OF f(C), FOR A GIVEN          // 
//	 POLYNOMIAL f, IN THE BGV SCHEME.             //
//                                                //
////////////////////////////////////////////////////    

// AUXILIARY FUNCTIONS 

//Auxiliary function used in evalAddCountForIntegerCoefficient
uint64_t popcountUint64(uint64_t x) {
    uint64_t count = 0;
    while (x > 0) {
        count += (x & 1ULL);
        x >>= 1;
    }
    return count;
}

//Auxiliary function used in evalAddCountForIntegerCoefficient
uint64_t floorLog2Uint64(uint64_t x) {
    uint64_t result = 0;
    while (x > 1) {
        x >>= 1;
        ++result;
    }
    return result;
}

//Auxiliary function used in evalAddCountForIntegerCoefficient
uint64_t absInt64AsUint64(int64_t x) {
    if (x < 0) {
        return static_cast<uint64_t>(-(x + 1)) + 1ULL;
    }
    return static_cast<uint64_t>(x);
}

//Auxiliary function used in estimateEvalAddCountPolynomial
uint64_t evalAddCountForIntegerCoefficient(int64_t coefficient) {
    if (coefficient == 0 || coefficient == 1) {
        return 0;
    }

    uint64_t absCoefficient = absInt64AsUint64(coefficient);

    uint64_t count =
        floorLog2Uint64(absCoefficient)
        + popcountUint64(absCoefficient)
        - 1;

    if (coefficient < 0) {
        count += 2;
    }

    return count;
}

// This function provide the number of addition required in the BGV scheme for evaluating the polynomial f, with coefficients {a0, a1, ..., ad}
uint64_t estimateEvalAddCountPolynomial(const std::vector<int64_t>& coefficients) {
    uint64_t count = 0;
    uint64_t nonzeroNonconstantTerms = 0;

    for (size_t k = 1; k < coefficients.size(); ++k) {
        int64_t coefficient = coefficients[k];
        if (coefficient == 0) {
            continue;
        }

        ++nonzeroNonconstantTerms;
        count += evalAddCountForIntegerCoefficient(coefficient);
    }

    if (nonzeroNonconstantTerms > 1) {
        count += nonzeroNonconstantTerms - 1;
    }

    if (!coefficients.empty() && coefficients[0] != 0) {
        if (nonzeroNonconstantTerms > 0) {
            count += 1;
        } else {
            count += 2;
        }
    }

    if (count == 0) {
        count = 1;
    }

    return count;
}

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

// Return true if exponent is a power of two
bool isPowerOfTwo(uint64_t exponent) {
    return exponent != 0 && (exponent & (exponent - 1)) == 0;
}

// Multiplicative depth required by the square-and-multiply routine.
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


// Degree of a polynomial represented by coefficients {a0, a1, ..., ad},
int64_t polynomialDegree(const std::vector<int64_t>& coefficients) {
    if (coefficients.empty()) {
        throw std::invalid_argument("The polynomial must contain at least one coefficient.");
    }

    for (int64_t i = static_cast<int64_t>(coefficients.size()) - 1; i >= 0; --i) {
        if (coefficients[static_cast<size_t>(i)] != 0) {
            return i;
        }
    }

    return -1;
}

// Multiplicative depth required to evaluate f(C)
uint32_t multiplicativeDepthPolynomial(const std::vector<int64_t>& coefficients) {
    int64_t degree = polynomialDegree(coefficients);
    if (degree <= 1) {
        return 1;
    }

    return multiplicativeDepthSquareAndMultiply(static_cast<uint64_t>(degree));
}

// Print the multiplicative depth required by square-and-multiply for each
// nonzero power C^k appearing in the polynomial.
void printMultiplicativeDepthPerPower(const std::vector<int64_t>& coefficients) {
    int64_t degree = polynomialDegree(coefficients);

    bool foundNonconstantPower = false;
    uint32_t maxPowerDepth = 0;

    for (int64_t k = 1; k <= degree; ++k) {
        int64_t coefficient = coefficients[static_cast<size_t>(k)];
        if (coefficient == 0) {
            continue;
        }

        uint32_t powerDepth = multiplicativeDepthSquareAndMultiply(static_cast<uint64_t>(k));
        maxPowerDepth = std::max(maxPowerDepth, powerDepth);
        foundNonconstantPower = true;

        std::cout << "  C";
        if (k >= 2) {
            std::cout << "^" << k;
        }
        std::cout << ": Multiplicative depth for the power = " << powerDepth;
                  

        if (k == 1) {
            std::cout << " (no EvalMult for C)";
        }

        std::cout << std::endl;
    }

    if (!foundNonconstantPower) {
        std::cout << "  No non-constant power in the polynomial." << std::endl;
    }

}

// Plaintext encoding of a constant term.
Plaintext makeConstantPlaintext(
    const CryptoContext<DCRTPoly>& cryptoContext,
    int64_t value,
    size_t plaintextSize) {

    std::vector<int64_t> constantVector(plaintextSize, value);
    return cryptoContext->MakePackedPlaintext(constantVector);
}

// Encrypted zero with the same basic parameters as a reference ciphertext.
Ciphertext<DCRTPoly> makeEncryptedZeroLike(
    const CryptoContext<DCRTPoly>& cryptoContext,
    const Ciphertext<DCRTPoly>& reference) {

    return cryptoContext->EvalSub(reference, reference);
}

// Multiplication of a ciphertext by an integer coefficient using double-and-add.
Ciphertext<DCRTPoly> multiplyCiphertextByIntegerCoefficient(
    const CryptoContext<DCRTPoly>& cryptoContext,
    const Ciphertext<DCRTPoly>& ciphertext,
    int64_t coefficient) {

    if (coefficient == 0) {
        return makeEncryptedZeroLike(cryptoContext, ciphertext);
    }

    if (coefficient == 1) {
        return ciphertext;
    }

    uint64_t absCoefficient;
    if (coefficient < 0) {
        absCoefficient = static_cast<uint64_t>(-(coefficient + 1)) + 1ULL;
    } else {
        absCoefficient = static_cast<uint64_t>(coefficient);
    }

    Ciphertext<DCRTPoly> result;
    Ciphertext<DCRTPoly> addend = ciphertext;
    bool resultInitialized = false;

    while (absCoefficient > 0) {
        if ((absCoefficient & 1ULL) != 0) {
            if (!resultInitialized) {
                result = addend;
                resultInitialized = true;
            } else {
                result = cryptoContext->EvalAdd(result, addend);
            }
        }

        absCoefficient >>= 1;
        if (absCoefficient > 0) {
            addend = cryptoContext->EvalAdd(addend, addend);
        }
    }

    if (coefficient < 0) {
        return cryptoContext->EvalSub(makeEncryptedZeroLike(cryptoContext, ciphertext), result);
    }

    return result;
}

// Evaluate f(C) = a0 + a1*C + ... + ad*C^d.
Ciphertext<DCRTPoly> evaluatePolynomial(
    const CryptoContext<DCRTPoly>& cryptoContext,
    const Ciphertext<DCRTPoly>& base,
    const std::vector<int64_t>& coefficients,
    size_t plaintextSize) {

    int64_t degree = polynomialDegree(coefficients);

    Ciphertext<DCRTPoly> result;
    bool resultInitialized = false;

    for (int64_t k = 1; k <= degree; ++k) {
        int64_t coefficient = coefficients[static_cast<size_t>(k)];
        if (coefficient == 0) {
            continue;
        }

        Ciphertext<DCRTPoly> power = (k == 1)
            ? base
            : squareAndMultiply(cryptoContext, base, static_cast<uint64_t>(k));

        Ciphertext<DCRTPoly> term = multiplyCiphertextByIntegerCoefficient(
            cryptoContext,
            power,
            coefficient);

        if (!resultInitialized) {
            result = term;
            resultInitialized = true;
        } else {
            result = cryptoContext->EvalAdd(result, term);
        }
    }


    int64_t constantCoefficient = coefficients[0];
    if (constantCoefficient != 0) {
        Plaintext constantPlaintext = makeConstantPlaintext(
            cryptoContext,
            constantCoefficient,
            plaintextSize);

        if (!resultInitialized) {
            result = cryptoContext->EvalAdd(
                makeEncryptedZeroLike(cryptoContext, base),
                constantPlaintext);
            resultInitialized = true;
        } else {
            result = cryptoContext->EvalAdd(result, constantPlaintext);
        }
    }

    if (!resultInitialized) {
        result = makeEncryptedZeroLike(cryptoContext, base);
    }

    return result;
}

std::string polynomialToString(const std::vector<int64_t>& coefficients) {
    int64_t degree = polynomialDegree(coefficients);
    if (degree == -1) {
        return "0";
    }

    std::ostringstream oss;
    bool firstTerm = true;

    for (int64_t k = 0; k <= degree; ++k) {
        int64_t coefficient = coefficients[static_cast<size_t>(k)];
        if (coefficient == 0) {
            continue;
        }

        if (!firstTerm) {
            oss << (coefficient > 0 ? " + " : " - ");
        } else if (coefficient < 0) {
            oss << "-";
        }

        uint64_t absCoefficient = coefficient < 0
            ? static_cast<uint64_t>(-(coefficient + 1)) + 1ULL
            : static_cast<uint64_t>(coefficient);
        bool showCoefficient = (k == 0 || absCoefficient != 1);

        if (showCoefficient) {
            oss << absCoefficient;
        }

        if (k >= 1) {
            if (showCoefficient) {
                oss << "*";
            }
            oss << "C";
            if (k >= 2) {
                oss << "^" << k;
            }
        }

        firstTerm = false;
    }

    return oss.str();
}


///////// MAIN PROGRAM ////////////

int main() {
	auto start = std::chrono::high_resolution_clock::now();	
	
	int64_t n = 8192; // Degree of the polynomial X^n+1, i.e., Ring Dimension
	uint64_t samples = 1000; // Number of samples
	int aux_size = 12; // Length of vectors used for creating the plaintexts

	// Definition of the polynomial f(C) = a0 + a1*C + ... + ad*C^d.
	// Insert the coefficients manually in increasing degree order: {a0, a1, ..., ad}.
	std::vector<int64_t> polynomialCoefficients = {0, 128, 1, 1};

	uint32_t M = std::max<uint32_t>(1, multiplicativeDepthPolynomial(polynomialCoefficients)); // Multiplicative depth required to evaluate f(C)
	
	std::cout << "Polynomial f(C): " << polynomialToString(polynomialCoefficients) << std::endl;
	printMultiplicativeDepthPerPower(polynomialCoefficients);
	std::cout << "Multiplicative depth in CryptoContext: " << M << std::endl;
	std::cout << std::endl;
	
	NTL::ZZ sum1; // Variable where we store the sum of squares of the first coefficient of the critical value for all sampled critical values
	
	//Set cryptographic context
    CCParams<CryptoContextBGVRNS> parameters;
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetStandardDeviation(3.2);
    parameters.SetRingDim(n);
    parameters.SetPlaintextModulus(65537);
    parameters.SetMultiplicativeDepth(M);
	
	uint32_t A = static_cast<uint32_t>(
    std::max<uint64_t>(1, estimateEvalAddCountPolynomial(polynomialCoefficients))
						);
	std::cout << "Add Count in Cryptocontext: " << A << std::endl;
	parameters.SetEvalAddCount(A); //Number of addition required in Cryptocontext
		
	CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
	
    cryptoContext->Enable(PKE);
	cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
		

	//Start loop
    for (size_t i = 0; i < samples; ++i) {

        //Key generation
        KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
		cryptoContext->EvalMultKeyGen(keyPair.secretKey);
        
        
        //Generation of plaintext
		auto rand_vec = generateRandomVector(aux_size);
		Plaintext pt = cryptoContext->MakePackedPlaintext(rand_vec);
		
       		
		// Compute f(C), where C is the ciphertext.
		auto ciphertext = cryptoContext->Encrypt(keyPair.publicKey, pt);
		// Modulo extraction
		auto f00 = ciphertext->GetElements()[0];
		f00.SetFormat(Format::COEFFICIENT);
		auto f0 = f00.CRTInterpolate();
		auto qf = f0.GetModulus();
		//std::cout << "Modulo: " << qf << std::endl;
		
		auto fCiphertext = evaluatePolynomial(cryptoContext, ciphertext, polynomialCoefficients, aux_size);
		
		// Start computation critical value of f(C)
		
		DCRTPoly s = keyPair.secretKey->GetPrivateElement();
		s.SetFormat(Format::COEFFICIENT);
		auto sk = s.CRTInterpolate(); //extraction of the secret key
		
        auto ciphertext00 = fCiphertext->GetElements()[0];
        auto ciphertext11 = fCiphertext->GetElements()[1];
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
		
		
		// Evaluate sum of the squares of the first coefficient of the critical value for all the sampled critical values
        sum1 += criticalValue[1] * criticalValue[1];
			
				
	}
	//End loop
	
	double sum1Double = to_double(sum1);
    double var1 = sum1Double / static_cast<double>(samples); //variance of the first coefficient of the critical value
    double logvar1 = log2(var1); //log of the variance of the first coefficient of the critical value
	
	
	//write the result in the file expPolynomial.txt
	std::string filename = "expPolynomial.txt";
            std::ofstream outfile(filename, std::ios::app);
            if (outfile.is_open()) {
                outfile << logvar1 << std::endl;
                outfile.close();
            } else {
                std::cerr << "Impossible open file " << filename << " to write." << std::endl;
			}

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    
	std::cout << "Var exp: " << logvar1 << std::endl;
    std::cout << "Execution time: " << duration.count() << " seconds" << std::endl;

    return 0;
}

