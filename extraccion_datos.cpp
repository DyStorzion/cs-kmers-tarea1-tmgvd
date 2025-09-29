#include "DNAsequence.cpp"
#include <unordered_map>

DNASequence canonicalKMer(DNASequence &&sequence){
    bool original = true;
    size_t length = sequence.getLength();

    DNASequence result;
    result.getChain().reserve(length);

    for (long long i = length - 1 ; i >= 0 ; i--){
        uint8_t base = sequence.getCode(length-1 - i);
        uint8_t reverseBase = DNASequence::reverseBase(sequence.getCode(i));

        if (original){
            if (base < reverseBase){
                return std::move(sequence);
            }
            if (base > reverseBase){
                original = false;
            }
        }
        
        result.push_back(reverseBase);
    }

    return std::move(result);
}

int main(){
    DNASequence adn;

    // Crear mecanismo de inserción

    //Suponiendo que está creada la cadena
    std::unordered_map<std::vector<uint8_t>, std::pair<DNASequence, int>, DNASequence::VectorHash> k21mers;
    std::unordered_map<std::vector<uint8_t>, std::pair<DNASequence, int>, DNASequence::VectorHash> k31mers;
    
    // Genera el Ground truth
    for (int i=0 ; i<adn.getLength()-20 ; i++){
        DNASequence kmer = canonicalKMer( std::move( adn.subSequence(i, i + 20) ) );

        std::vector<uint8_t> sequence = kmer.getChain();
        auto pair = k21mers.insert({sequence, {std::move(kmer), 1}});
        // si la clave ya existe se incrementa contador
        if (!pair.second) {
            ++(pair.first->second.second);
        }
    }
    for (int i=0 ; i<adn.getLength()-30 ; i++){
        DNASequence kmer = canonicalKMer( std::move( adn.subSequence(i, i + 30) ) );

        std::vector<uint8_t> sequence = kmer.getChain();
        auto pair = k31mers.insert({sequence, {std::move(kmer), 1}});
        if (!pair.second){
            ++(pair.first->second.second);
        }
    }

    // Generar lista de heavy hitters exactos
    std::vector<std::pair<DNASequence, int>> heavyHitters21;
    std::vector<std::pair<DNASequence, int>> heavyHitters31;

    // Umbral y cantidad de k-mers
    int k21mersAmount = adn.getLength()-20;
    int k31mersAmount = adn.getLength()-30;

    float k21mersThreshold = 1e-3f;
    float k31mersThreshold = 1e-3f;

    int k21mersBoundary = k21mersAmount * k21mersThreshold;
    int k31mersBoundary = k31mersAmount * k21mersThreshold;

    // Agrega el 21-mer a la lista de heavy hitters si supera la frecuencia definida por el umbral
    for (auto &kv : k21mers) {
        auto &seq = kv.second.first;
        int count = kv.second.second;
        if (count >= k21mersBoundary) {
            heavyHitters21.emplace_back(seq, count);
        }
    }

    // Agrega el 31-mer a la lista de heavy hitters si supera la frecuencia definida por el umbral
    for (auto &kv : k31mers) {
        auto &seq = kv.second.first;
        int count = kv.second.second;
        if (count >= k31mersBoundary) {
            heavyHitters31.emplace_back(seq, count);
        }
    }
}