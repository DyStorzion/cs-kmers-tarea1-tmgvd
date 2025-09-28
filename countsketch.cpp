#include "murmurhash32.hpp"

class CountSketch
{
private:
    int d,w; 
    std::vector<std::vector<int>> tabla;

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
    // Crea el countsketch con d filas y w columnas
    CountSketch(int d, int w): d(d), w(w), tabla(d, std::vector<int>(w,0)) {};

    // Inserta un k-mer en el countsketch
    void insert(const std::string &kmer) {
        std::string canon = canonical(kmer); // usar k-mer canónico
        for (int j = 0; j < d; j++) {
            uint32_t h_j = murmurhash(canon, j) % w; // hash para la columna
            uint32_t s_j = murmurhash(canon, j + 1000); // hash para el signo
            int sign = (s_j & 1) ? 1 : -1; // signo basado en el hash

            tabla[j][h_j] += sign; // actualizar la tabla
        }
    }

    // Estima la frecuencia de un k-mer en el countsketch
    int estimate(const std::string &kmer) {
        std::string canon = canonical(kmer); // usar k-mer canónico
        std::vector<int> estimates; 
        estimates.reserve(d); 

        for (int j = 0; j < d; j++) {
            uint32_t h_j = murmurhash(canon, j) % w; // hash para la columna
            uint32_t s_j = murmurhash(canon, j + 1000); // hash para el signo
            int sign = (s_j & 1) ? 1 : -1; // signo basado en el hash

            estimates.push_back(sign * tabla[j][h_j]);
        }

        // retornar la mediana de las estimaciones
        std::nth_element(estimates.begin(), estimates.begin() + estimates.size()/2, estimates.end());
        return estimates[estimates.size()/2];
    }
};