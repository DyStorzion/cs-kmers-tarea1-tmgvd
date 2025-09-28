#include "murmurhash32.hpp"
#include <climits>

class TowerSketch
{
private:
    int d,w,m; 
    std::vector<std::vector<std::vector<int>>> tabla;

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
    TowerSketch(int d, int w, int m) : d(d), w(w), m(m), tabla(d, std::vector<std::vector<int>>(m, std::vector<int>(w, 0))) {}

    // Estimación mínima de un k-mer en toda la torre
    int estimate(const std::string &kmer) {
        std::string canon = canonical(kmer); // usar k-mer canónico
        int min_val = INT_MAX;
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < m; j++) {
                uint32_t h = murmurhash(canon, i * 100 + j) % w; // hash para la columna
                min_val = std::min(min_val, tabla[i][j][h]); // actualizar el mínimo
            }
        }
        return min_val;
    }

    // Insertar k-mer con Conservative Update
    void insert(const std::string &kmer) {
        std::string canon = canonical(kmer);
        int est = estimate(canon); // estimación actual
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < m; j++) {
                uint32_t h = murmurhash(canon, i * 100 + j) % w; 
                
                tabla[i][j][h] = std::max(tabla[i][j][h], est + 1); // actualizar solo si es mayor
            }
        }
    }
    
};