#include "sketchs/countsketch.hpp"
#include "utils/LectorGenomas.hpp"

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