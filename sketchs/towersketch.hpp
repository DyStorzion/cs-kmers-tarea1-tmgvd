#ifndef TOWER_SKETCH_H
#define TOWER_SKETCH_H


#include "murmurhash32.hpp"
#include <climits>

/**
 * CountMin sketch con conservative update y template para el tamaño de los contadores
 * 
 */
template<typename T>
class CountMinCU
{
private:
    int d,w;
    std::vector<std::vector<T>> tabla;
    // retorna el reverso complementario de la cadena de ADN
std::string revComp(const std::string &s) {
    std::string rc = s;
    for (char &c : rc) {
        switch(c) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
        }
    }
    std::reverse(rc.begin(), rc.end());
    return rc;
}
// retorna el k-mer canónico (lexicográficamente mínimo entre la cadena y su reverso complementario)
std::string canonical(const std::string &s) {
    std::string rc = revComp(s);
    return std::min(s, rc); // lexicográficamente mínimo
}
public:
    CountMinCU(int d, int w): d(d), w(w), tabla(d, std::vector<T>(w,0)) {};

    void insert(const std::string &kmer) {
        std::string canon = canonical(kmer);
        int frec_estimada = estimate(kmer);
        for (int j = 0; j < d; j++) {
            if(tabla[j][murmurhash(canon, j) % w] == frec_estimada)
                tabla[j][murmurhash(canon, j) % w]++;
        }
    }

    T estimate(const std::string &kmer) {
        std::string canon = canonical(kmer);
        T frec_est = INT_MAX;
        for (int j = 0; j < d; j++) {
            frec_est = std::min(frec_est, tabla[j][murmurhash(canon, j) % w]);
        }
        return frec_est;
    }

    size_t getSize() {
        return d * w * sizeof(T);
    }

};

/**
 * Tower Sketch: combinación de CountMinCU con contadores de 8, 16 y 32 bits
 * 
 */
class TowerSketch
{
private:
    CountMinCU<uint8_t> countMin8;
    CountMinCU<uint16_t> countMin16;
    CountMinCU<uint32_t> countMin32;

public:
    TowerSketch(int d8, int w8, int d16, int w16, int d32, int w32):
      countMin8(CountMinCU<uint8_t>(d8, w8)),
      countMin16(CountMinCU<uint16_t>(d16, w16)),
      countMin32(CountMinCU<uint32_t>(d32, w32)) {}

    TowerSketch(int d, int w):
      countMin8(CountMinCU<uint8_t>(d, w)),
      countMin16(CountMinCU<uint16_t>(d, w)),
      countMin32(CountMinCU<uint32_t>(d, w)) {}

    void insert(const std::string &kmer) {
        
        uint8_t est8 = countMin8.estimate(kmer);
        if (est8 < UINT8_MAX) {
            countMin8.insert(kmer);
            return;
        }
        uint16_t est16 = countMin16.estimate(kmer);
        if (est16 < UINT16_MAX) {
            countMin16.insert(kmer);
            return;
        }
        countMin32.insert(kmer);
        return;
    }

    int estimate(const std::string &kmer) {
        uint32_t est32 = countMin32.estimate(kmer);
        if (est32 > 0) return est32;

        uint16_t est16 = countMin16.estimate(kmer);
        if (est16 > 0) return est16;

        uint8_t est8 = countMin8.estimate(kmer);
        return est8;
    }

    size_t getSize() {
        return countMin8.getSize() + countMin16.getSize() + countMin32.getSize();
    }
};

#endif // TOWER_SKETCH_H