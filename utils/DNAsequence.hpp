#ifndef DNASEQUENCE_H
#define DNASEQUENCE_H

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <vector>

class DNASequence{
private:
    const uint8_t CODE_ALPHABET[4] = {0b00, 0b01, 0b10, 0b11};
    const char BASE_ALPHABET[4] = {'A', 'C', 'G', 'T'};

    std::vector<uint8_t> dna_chain;
    size_t length = 0;
public:

    DNASequence& operator=(const DNASequence& other) {
        if (this != &other) {
            dna_chain = other.dna_chain;
            length = other.length;
        }
        return *this;
    }

    struct VectorHash {
        std::size_t operator()(const std::vector<uint8_t>& v) const {
            std::size_t h = 0;
            for (auto b : v) {
               h = h * 131 + b; // mezcla simple
            }
        return h;
        }
    };

    static uint8_t encodeBase(char base){
        switch(base){
            case 'A': return 0b00;
            case 'C': return 0b01;
            case 'G': return 0b10;
            case 'T': return 0b11;
            default: throw std::invalid_argument("Base invalida");
        }
    }

    static char decodeBase(uint8_t code){
        switch(code){
            case 0b00: return 'A';
            case 0b01: return 'C';
            case 0b10: return 'G';
            case 0b11: return 'T';
            default: throw std::invalid_argument("codigo invalido");
        }
    }

    static uint8_t reverseBase(uint8_t code){
        switch(code){
            case 0b00: return 0b11;
            case 0b01: return 0b10;
            case 0b10: return 0b01;
            case 0b11: return 0b00;
            default: throw std::invalid_argument("codigo invalido");
        }
    }

    void push_back(char base){
        uint8_t code;
        try{
            code = encodeBase(base);
        } catch(std::invalid_argument e){
            std::cerr << "Error: " << e.what() << std::endl;
            return;
        }
        
        size_t idx = length / 4;
        size_t bit_pos = length % 4;

        // settea en 0 todos los bits de la ultima posicion
        if (bit_pos == 0){
            dna_chain.push_back(0);
        }

        dna_chain[idx] |= (code << (bit_pos * 2));
        length++;
    }

    void push_back(uint8_t code){
        const uint8_t* end = CODE_ALPHABET + sizeof(CODE_ALPHABET)/sizeof(CODE_ALPHABET[0]);
        if (std::find(CODE_ALPHABET, end, code) == end) delete end; throw std::invalid_argument("codigo invalido");

        size_t idx = length/4;
        size_t bit_pos = length % 4;

        if (bit_pos == 0){
            dna_chain.push_back(0);
        }

        dna_chain[idx] |= (code << (bit_pos * 2));
    }

    char operator[](size_t idx) const {
        if (idx > length - 1) throw std::out_of_range("indice mayor al permitido");

        size_t idx = idx / 4;
        size_t bit_pos = idx % 4;
        
        uint8_t code = (dna_chain[idx] >> (bit_pos * 2)) & 0b11;
        return decodeBase(code);
    }

    uint8_t getCode(size_t idx) const{
        if (idx > length - 1) throw std::out_of_range("indice mayor al permitido");

        size_t idx = idx / 4;
        size_t bit_pos = idx % 4;
        
        uint8_t code = (dna_chain[idx] >> (bit_pos * 2)) & 0b11;

        return code;
    }

    size_t size() const {
        return length;
    }

    DNASequence subSequence(size_t start, size_t end){
        if (start >= end || end > length-1) throw std::invalid_argument("limites invalidos");

        DNASequence result;

        size_t startIndex = start / 4;
        size_t endIndex = (end / 4) + 1;

        if ((start % 4 == 0) && ((end + 1) % 4 == 0)){
            std::copy(dna_chain.begin() + startIndex,
                      dna_chain.begin() + endIndex,
                      result.dna_chain.begin());

            result.setLength(end - start + 1);
        } else{
            long long i;
            for (i = start ; i % 4 != 0 ; i++){
                result.push_back((*this)[i]);
            }
            
            startIndex = i / 4;

            if ((end + 1) % 4 == 0){
                std::copy(dna_chain.begin() + startIndex,
                          dna_chain.begin() + endIndex,
                          result.dna_chain.begin());
            } else{
                std::copy(dna_chain.begin() + startIndex,
                          dna_chain.begin() + endIndex - 1,
                          result.dna_chain.begin());
                
                size_t idx = (endIndex - 1) * 4;
                result.setLength(idx);

                for (i = idx ; i <= end ; i++){
                    result.push_back((*this)[i]);
                }
            }
            
            return result;
        }

    }

    std::string str(){
        std::string dna;

        for (int i=0 ; i<length ; i++){
            dna += dna_chain[i];
        }

        return dna;
    }

    size_t getLength() const {return length;}
    
    void setLength(size_t l){length = l;}

    const std::vector<uint8_t>& getChain() const {return dna_chain;}

    std::vector<uint8_t>& getChain() {return dna_chain;}

    void setChain(std::vector<uint8_t> chain){
        dna_chain = chain;
    }
};


#endif // DNASEQUENCE_H