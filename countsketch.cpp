#include "murmurhash32.hpp"
#include "LectorGenomas.cpp"
#include <unordered_map>
#include <iostream>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <vector>
#include <string>
#include <stdexcept>

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

int main() {
    try {
        // Parámetros del CountSketch
        int d = 8;       // Número de funciones hash
        int w = 50000;   // Tamaño de cada tabla hash
        int k = 21;      // Longitud de k-mers
        
        CountSketch sketch(d, w);
    
        LectorGenomas reader("Genomas");
        std::cout << "Archivos FASTA encontrados: " << reader.getTotalFiles() << std::endl;
        
        // Estadísticas
        long long totalKmers = 0;
        int processedFiles = 0;
        
        
        // Colección de k-mers candidatos para verificar
        std::unordered_map<std::string, bool> candidateKmers;
        
        do {
            std::cout << "Procesando archivo " << (processedFiles + 1) << ": " 
                     << reader.getCurrentFilename() << std::endl;
            
            long long fileKmers = 0;
            int candidateCounter = 0;
            
            while (reader.hasMoreKmers(k) && fileKmers < 2000000) { // Límite por archivo
                std::string kmer = reader.getNextKmer(k);
                if (!kmer.empty()) {
                    // Verificar bases válidas
                    bool validKmer = true;
                    for (char c : kmer) {
                        if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
                            validKmer = false;
                            break;
                        }
                    }
                    
                    if (validKmer) {
                        sketch.insert(kmer);
                        fileKmers++;
                        totalKmers++;
                        
                        // Recolectar k-mers como candidatos a heavy hitters
                        if (candidateCounter % 1000 == 0) { // Cada 1000 k-mers
                            candidateKmers[kmer] = true;
                        }
                        candidateCounter++;
                        
                        // Progreso
                        if (totalKmers % 500000 == 0) {
                            std::cout << "  Procesados: " << totalKmers << " k-mers..." << std::endl;
                        }
                    }
                }
            }
            
            std::cout << "  K-mers procesados: " << fileKmers << std::endl;
            processedFiles++;
            
        } while (reader.nextFile());
        
        
        std::cout << "\n=== Estadísticas del procesamiento ===" << std::endl;
        std::cout << "Total de k-mers procesados: " << totalKmers << std::endl;
        std::cout << "Archivos procesados: " << processedFiles << std::endl;
        std::cout << "K-mers candidatos a evaluar: " << candidateKmers.size() << std::endl;
        
        // EXTRACCIÓN DE HEAVY HITTERS
        std::cout << "\n=== Extrayendo Heavy Hitters ===" << std::endl;
        
        std::vector<std::pair<std::string, int>> heavyHitters;
        int k21mersAmount = candidateKmers.size();
        float k21mersThreshold = 1e-3f;

        int k21mersBoundary = k21mersAmount * k21mersThreshold;

        // Evaluar cada k-mer candidato
        for (const auto& pair : candidateKmers) {
            const std::string& kmer = pair.first;
            int estimatedFreq = sketch.estimate(kmer);
            
            if (estimatedFreq >= k21mersBoundary) {
                heavyHitters.emplace_back(kmer, estimatedFreq);
            }
        }
        
        // Ordenar por frecuencia descendente
        std::sort(heavyHitters.begin(), heavyHitters.end(), 
                 [](const auto& a, const auto& b) { return a.second > b.second; });
        
        // Mostrar resultados
        std::cout << "\n=== RESULTADOS: Heavy Hitters encontrados ===" << std::endl;
        std::cout << "Total de heavy hitters (freq >= " << k21mersBoundary << "): " << heavyHitters.size() << std::endl;
        
        if (!heavyHitters.empty()) {
            std::cout << "\nTop Heavy Hitters:" << std::endl;
            std::cout << "Rank\tK-mer\t\t\tFrecuencia Estimada" << std::endl;
            std::cout << "----------------------------------------------------" << std::endl;
            
            int showCount = std::min(20, (int)heavyHitters.size()); // Mostrar top 20
            for (int i = 0; i < showCount; i++) {
                std::cout << (i+1) << "\t" << heavyHitters[i].first 
                         << "\t" << heavyHitters[i].second << std::endl;
            }
            
            if (heavyHitters.size() > 20) {
                std::cout << "... y " << (heavyHitters.size() - 20) << " más." << std::endl;
            }
        } else {
            std::cout << "\nNo se encontraron heavy hitters con umbral >= " << k21mersBoundary << std::endl;
        }
        std::cout << "\n=== Extracción completada ===" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}