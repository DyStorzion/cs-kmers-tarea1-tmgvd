#include "sketchs/countsketch.hpp"
#include "utils/LectorGenomas.hpp"
#include <unordered_set>
#include <algorithm>

// Función para obtener k-mer canónico
std::string getCanonicalKmer(const std::string& kmer) {
    std::string revComp = kmer;
    
    // Crear complemento reverso
    for (char& c : revComp) {
        switch(c) {
            case 'A': c = 'T'; break;
            case 'T': c = 'A'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
        }
    }
    std::reverse(revComp.begin(), revComp.end());
    
    // Retornar el lexicográficamente menor
    return std::min(kmer, revComp);
}

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
        
        
        // Colección de TODOS los k-mers únicos para evaluación completa
        std::unordered_set<std::string> uniqueKmers;
        
        do {
            std::cout << "Procesando archivo " << (processedFiles + 1) << ": " 
                     << reader.getCurrentFilename() << std::endl;
            
            long long fileKmers = 0;
            
            while (reader.hasMoreKmers(k) && fileKmers < 2000000) { // Límite por archivo para pruebas
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
                        // Obtener k-mer canónico
                        std::string canonical = getCanonicalKmer(kmer);
                        sketch.insert(canonical);
                        uniqueKmers.insert(canonical);  // Guardar todos los k-mers únicos
                        fileKmers++;
                        totalKmers++;
                        
                        // Progreso
                        if (totalKmers % 500000 == 0) {
                            std::cout << "  Procesados: " << totalKmers << " k-mers, Únicos: " << uniqueKmers.size() << std::endl;
                        }
                    }
                }
            }
            
            std::cout << "  K-mers procesados: " << fileKmers << std::endl;
            processedFiles++;
            
        } while (reader.nextFile());
        
        
        std::cout << "\n=== Estadísticas del procesamiento ===" << std::endl;
        std::cout << "Total de k-mers procesados: " << totalKmers << std::endl;
        std::cout << "K-mers únicos encontrados: " << uniqueKmers.size() << std::endl;
        std::cout << "Archivos procesados: " << processedFiles << std::endl;
        
        // EXTRACCIÓN DE HEAVY HITTERS
        std::cout << "\n=== Extrayendo Heavy Hitters ===" << std::endl;
        
        // Calcular umbral correcto basado en el total de k-mers procesados
        double phi = 2e-6;  // Umbral que sabemos funciona bien para 21-mers
        int k21mersBoundary = static_cast<int>(phi * totalKmers);
        
        std::cout << "Umbral φ = " << phi << std::endl;
        std::cout << "Umbral frecuencia = " << k21mersBoundary << std::endl;
        
        std::vector<std::pair<std::string, int>> heavyHitters;

        // Evaluar TODOS los k-mers únicos
        std::cout << "\nEvaluando " << uniqueKmers.size() << " k-mers únicos..." << std::endl;
        int evaluatedCount = 0;
        
        for (const std::string& kmer : uniqueKmers) {
            int estimatedFreq = sketch.estimate(kmer);
            
            if (estimatedFreq >= k21mersBoundary) {
                heavyHitters.emplace_back(kmer, estimatedFreq);
            }
            
            evaluatedCount++;
            if (evaluatedCount % 50000 == 0) {
                std::cout << "Evaluados: " << evaluatedCount << "/" << uniqueKmers.size() 
                         << ", HH encontrados: " << heavyHitters.size() << std::endl;
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