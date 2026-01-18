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
//   THIS PROGRAM EVALUATE THE VARIANCE OF THE    //
//   FIRST COEFFICIENT OF THE CRITICAL VALUE OF   //
//   A FRESH CIPHERTEXT AFTER MODULO SWITCH       //
//   IN THE BGV SCHEME.                           //
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


///////// MAIN PROGRAM ////////////

int main() {
	

	int64_t n = 8192; // Degree of the polynomial X^n+1, i.e., Ring Dimension
	uint64_t samples = 50000; // Number of samples
	int aux_size = 12; // Length of vectors used for creating the plaintexts
	
    NTL::ZZ sum1; // Variable where we store the sum of squares of the first coefficient of the critical value for all sampled critical values
	
	//Set cryptographic context
    CCParams<CryptoContextBGVRNS> parameters;
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetStandardDeviation(3.2);
    parameters.SetRingDim(n);
    parameters.SetPlaintextModulus(65537);
    parameters.SetMultiplicativeDepth(1);
	
	parameters.SetScalingTechnique(FIXEDMANUAL); //Modulo switch disabled
	
	CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
	
    cryptoContext->Enable(PKE);
	cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);
		
		
	//Start loop
    for (size_t i = 0; i < samples; ++i) {

   
        //Key generation
        KeyPair<DCRTPoly> keyPair = cryptoContext->KeyGen();
		cryptoContext->EvalMultKeyGen(keyPair.secretKey);
        
        //Generation of a random plaintext
		auto rand_vec = generateRandomVector(aux_size);
        Plaintext pt = cryptoContext->MakePackedPlaintext(rand_vec);
		
        auto ct = cryptoContext->Encrypt(keyPair.publicKey, pt); //Encryption of the plaintext
		auto ciphertext = cryptoContext->ModReduce(ct); //Modulo switch of the ciphertext
		
		// Start computation critical value of the final ciphertext
		
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
		
	// Calcolo della varianza e suo logaritmo del primo coefficiente del criticalValue
	double sum1Double = to_double(sum1);
    double var1 = sum1Double / static_cast<double>(samples); //variance of the first coefficient of the critical value
    double logvar1 = log2(var1); //log of the variance of the first coefficient of the critical value
	
	//write the result in the file expMS.txt
	std::string filename = "expMS.txt";
            std::ofstream outfile(filename, std::ios::app);
            if (outfile.is_open()) {
                outfile << logvar1 << std::endl;
                outfile.close();
            } else {
                std::cerr << "Impossible open file " << filename << " to write." << std::endl;
			}

    return 0;
}
