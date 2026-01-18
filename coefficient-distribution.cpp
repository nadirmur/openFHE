#include "openfhe.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <random>
#include <fstream>
#include <NTL/ZZ.h>


using namespace lbcrypto;


////////////////////////////////////////////////////
//                                                //
//   THIS PROGRAM COLLECT THE FIRST COEFFICIENT   //
//   OF THE CRITICAL VALUE OF A CIPHERTEXT        //
//   OBTAINED AFTER L MULTIPLICATION FOR L FROM   //
//   1 TO M (MULTIPLICATIVE DEPTH) FOR SEVERAL    //
//   SAMPLES IN ORDER TO CHECK THE DISTRIBUTION   //
//   OF SUCH COEFFICIENTS FOR THE BGV SCHEME      //
//                                                //
//////////////////////////////////////////////////// 


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

// Perform multiplication between 2^L plaintexts as in the circuit of Fig.4 of the paper and output the ciphertexts for each level
// "Accurate BGV Parameters Selection: Accounting for Secret and Public Key Dependencies in Average-Case Analysis"
// B. Biasioli, C. Marcolla, N. Murru, M. Urani
std::vector<std::vector<Ciphertext<DCRTPoly>>> buildMultiplicationLevels(
        CryptoContext<DCRTPoly> cryptoContext,
        const PublicKey<DCRTPoly>& pk,
        const std::vector<Plaintext>& plain_vec) {
    
    std::vector<Ciphertext<DCRTPoly>> ciphers;
    for (const auto& pt : plain_vec) {
        ciphers.push_back(cryptoContext->Encrypt(pk, pt));
    }

    std::vector<std::vector<Ciphertext<DCRTPoly>>> levels;
    levels.push_back(ciphers); // livello 0 (fresh)

    while (ciphers.size() > 1) {
        std::vector<Ciphertext<DCRTPoly>> next_level;
        for (size_t i = 0; i < ciphers.size(); i += 2) {
            auto mult = cryptoContext->EvalMult(ciphers[i], ciphers[i + 1]);
            next_level.push_back(mult);
        }
        levels.push_back(next_level);
        ciphers = next_level;
    }

    return levels;
}


///////// MAIN PROGRAM ////////////

int main() {

	int64_t n = 8192; // Degree of the polynomial X^n+1, i.e., Ring Dimension
	int M =6; // Multiplicative depth
	uint64_t samples = 50000; // Number of samples
	int num_inputs = 1 << M; //2^M initial plaintexts
	int aux_size = 12; // Length of vectors used for creating the plaintexts
    
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
        
	//Start loop
    for (size_t i = 0; i < samples; ++i) {
        
		//Key generation
		KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
        cryptoContext->EvalMultKeyGen(keyPair.secretKey);
        
        //Generation of 2^M random plaintexts
        std::vector<Plaintext> plaintexts;
        for (int j = 0; j < num_inputs; ++j) {
            auto rand_vec = generateRandomVector(aux_size);
            Plaintext pt = cryptoContext->MakePackedPlaintext(rand_vec);
            plaintexts.push_back(pt);
        }
        
        // Perform the multiplication between all the plaintexts through the levels and output the ciphertexts of each level
		// see fig. 4 paper "Accurate BGV Parameters Selection: Accounting for Secret and Public Key Dependencies in Average-Case Analysis"
		// B. Biasioli, C. Marcolla, N. Murru, M. Urani
        auto levels = buildMultiplicationLevels(cryptoContext, keyPair.publicKey, plaintexts);

		//Computation of the critical value of ciphertexts in each level
        for (size_t lvl = 0; lvl < levels.size(); ++lvl) {
            auto ciphertext = levels[lvl][0]; //take first ciphertext of the conisdered level

            DCRTPoly s = keyPair.secretKey->GetPrivateElement();
            s.SetFormat(Format::COEFFICIENT);
            auto sk = s.CRTInterpolate(); //extraction of the secret key

            auto ciphertext00 = ciphertext->GetElements()[0];
            auto ciphertext11 = ciphertext->GetElements()[1];
            ciphertext00.SetFormat(Format::COEFFICIENT);
            ciphertext11.SetFormat(Format::COEFFICIENT);
            auto c0 = ciphertext00.CRTInterpolate(); //extraction of the first polynomial c0 from the ciphertext
            auto c1 = ciphertext11.CRTInterpolate(); //extraction of the second polynimial c1 from the ciphertext

			// Start transform type of coefficients of the polynomial sk (secret key) c0 and c1 (polynomials of the ciphertext) in NTL::ZZ in order to manage large integers
            std::vector<NTL::ZZ> c0int, c1int, skint;
            for (size_t j = 0; j < c0.GetLength(); ++j) {
                c0int.push_back(NTL::conv<NTL::ZZ>(c0[j].ToString().c_str()));
                c1int.push_back(NTL::conv<NTL::ZZ>(c1[j].ToString().c_str()));
                skint.push_back(NTL::conv<NTL::ZZ>(sk[j].ToString().c_str()));
            }
			// End transformation

            auto qtemp = c0.GetModulus();
            NTL::ZZ q = NTL::conv<NTL::ZZ>(qtemp.ToString().c_str()); //extraction of ciphertext modulo and conversion in type NTL::ZZ

			// Computation of the critical value: criticalValue = c0 + c1 sk
            std::vector<NTL::ZZ> c1s = multiplyPolynomials(c1int, skint);
            std::vector<NTL::ZZ> CV = addVectors(c0int, c1s);
            std::vector<NTL::ZZ> CVModPol = reducePolynomial(CV, n);
            std::vector<NTL::ZZ> criticalValue = reduceVectorModulo(CVModPol, q);

			// Save the first coefficient of the critical value in the file critical_values_level_L where L=lvl is the current level 
            std::string filename = "critical_values_level" + std::to_string(lvl) + ".txt";
            std::ofstream outfile(filename, std::ios::app);
            if (outfile.is_open()) {
                outfile << criticalValue[1] << std::endl;
                outfile.close();
            } else {
                std::cerr << "Impossible open file " << filename << " to write." << std::endl;
            }
        }
    }

    return 0;
}

