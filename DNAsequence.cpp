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

    void push_back(char base){
        uint8_t code = encodeBase(base);
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

    char operator[](size_t i) const {
        if (i > length - 1) throw std::out_of_range("indice mayor al permitido");

        size_t idx = i / 4;
        size_t bit_pos = i % 4;
        
        uint8_t code = (dna_chain[idx] >> (bit_pos * 2)) & 0b11;
        return decodeBase(code);
    }

    size_t size() const {
        return length;
    }

    void setLength(size_t l){length = l;}

    std::vector<uint8_t> getChain() {return dna_chain;}

    void setChain(std::vector<uint8_t> chain){
        dna_chain = chain;
    }
};