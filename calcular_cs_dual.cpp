#include "sketchs/countsketch.hpp"
#include "utils/LectorGenomas.hpp"
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <chrono>

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

// Función para procesar k-mers de una longitud específica
std::vector<std::pair<std::string, int>> procesarCountSketch(int k, double phi, const std::string& titulo) {
    std::cout << "\n=== " << titulo << " ===" << std::endl;
    
    // Parámetros del CountSketch
    int d = 8;       // Número de funciones hash
    int w = 50000;   // Tamaño de cada tabla hash
    
    CountSketch sketch(d, w);
    LectorGenomas reader("Genomas");
    
    // Estadísticas
    long long totalKmers = 0;
    int processedFiles = 0;
    std::unordered_set<std::string> uniqueKmers;
    
    std::cout << "Procesando k-mers de longitud " << k << "..." << std::endl;
    
    do {
        if (processedFiles % 5 == 0) {
            std::cout << "Procesando archivo " << (processedFiles + 1) << ": " 
                     << reader.getCurrentFilename() << std::endl;
        }
        
        long long fileKmers = 0;
        reader.reset(); // Reiniciar posición para este archivo
        
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
                    std::string canonical = getCanonicalKmer(kmer);
                    sketch.insert(canonical);
                    uniqueKmers.insert(canonical);
                    fileKmers++;
                    totalKmers++;
                    
                    // Progreso cada 1M k-mers
                    if (totalKmers % 1000000 == 0) {
                        std::cout << "  Procesados: " << totalKmers << " k-mers, Únicos: " << uniqueKmers.size() << std::endl;
                    }
                }
            }
        }
        
        processedFiles++;
        
    } while (reader.nextFile());
    
    std::cout << "\nEstadísticas " << k << "-mers:" << std::endl;
    std::cout << "Total procesados: " << totalKmers << std::endl;
    std::cout << "Únicos encontrados: " << uniqueKmers.size() << std::endl;
    std::cout << "Archivos procesados: " << processedFiles << std::endl;
    
    // Calcular umbral
    int umbralFrecuencia = static_cast<int>(phi * totalKmers);
    std::cout << "Umbral φ = " << phi << " → frecuencia >= " << umbralFrecuencia << std::endl;
    
    // Extraer heavy hitters
    std::vector<std::pair<std::string, int>> heavyHitters;
    std::cout << "Evaluando " << uniqueKmers.size() << " k-mers únicos..." << std::endl;
    
    int evaluatedCount = 0;
    for (const std::string& kmer : uniqueKmers) {
        int estimatedFreq = sketch.estimate(kmer);
        
        if (estimatedFreq >= umbralFrecuencia) {
            heavyHitters.emplace_back(kmer, estimatedFreq);
        }
        
        evaluatedCount++;
        if (evaluatedCount % 100000 == 0) {
            std::cout << "  Evaluados: " << evaluatedCount << "/" << uniqueKmers.size() 
                     << ", HH encontrados: " << heavyHitters.size() << std::endl;
        }
    }
    
    // Ordenar por frecuencia descendente
    std::sort(heavyHitters.begin(), heavyHitters.end(), 
             [](const auto& a, const auto& b) { return a.second > b.second; });
    
    std::cout << "Heavy hitters " << k << "-mers encontrados: " << heavyHitters.size() << std::endl;
    
    // Mostrar top 10
    if (!heavyHitters.empty()) {
        std::cout << "\nTop 10 Heavy Hitters " << k << "-mers:" << std::endl;
        std::cout << "Rank\tFrecuencia\tK-mer" << std::endl;
        std::cout << "----------------------------------------" << std::endl;
        
        int showCount = std::min(10, (int)heavyHitters.size());
        for (int i = 0; i < showCount; i++) {
            std::cout << (i+1) << "\t" << heavyHitters[i].second 
                     << "\t\t" << heavyHitters[i].first << std::endl;
        }
    }
    
    // Guardar en CSV
    std::string csvFilename = "countsketch_heavy_hitters_" + std::to_string(k) + "mers.csv";
    std::ofstream csvFile(csvFilename);
    
    if (csvFile.is_open()) {
        csvFile << "rank,kmer,estimated_frequency,threshold_used,total_kmers,phi_value,kmer_length\\n";
        
        for (size_t i = 0; i < heavyHitters.size(); i++) {
            csvFile << (i + 1) << "," 
                   << heavyHitters[i].first << "," 
                   << heavyHitters[i].second << "," 
                   << umbralFrecuencia << "," 
                   << totalKmers << "," 
                   << std::scientific << std::setprecision(1) << phi << ","
                   << k << "\\n";
        }
        
        csvFile.close();
        std::cout << "💾 Guardado en: " << csvFilename << std::endl;
        std::cout << "📊 Registros: " << heavyHitters.size() << std::endl;
    } else {
        std::cerr << "❌ Error creando CSV: " << csvFilename << std::endl;
    }
    
    return heavyHitters;
}

int main() {
    try {
        std::cout << "╔══════════════════════════════════════════════════════════════╗" << std::endl;
        std::cout << "║           COUNTSKETCH PARA 21-MERS Y 31-MERS               ║" << std::endl;
        std::cout << "╚══════════════════════════════════════════════════════════════╝" << std::endl;
        
        LectorGenomas testReader("Genomas");
        std::cout << "Archivos FASTA encontrados: " << testReader.getTotalFiles() << std::endl;
        
        // Procesar 21-mers
        auto heavyHitters21 = procesarCountSketch(21, 2e-6, "PROCESANDO 21-MERS CON COUNTSKETCH");
        
        // Procesar 31-mers  
        auto heavyHitters31 = procesarCountSketch(31, 6e-6, "PROCESANDO 31-MERS CON COUNTSKETCH");
        
        // Resumen final
        std::cout << "\n╔══════════════════════════════════════════════════════════════╗" << std::endl;
        std::cout << "║                    RESUMEN FINAL                           ║" << std::endl;
        std::cout << "╚══════════════════════════════════════════════════════════════╝" << std::endl;
        
        std::cout << "🎯 21-mers Heavy Hitters: " << heavyHitters21.size() << std::endl;
        std::cout << "🎯 31-mers Heavy Hitters: " << heavyHitters31.size() << std::endl;
        
        std::cout << "\n📁 Archivos generados:" << std::endl;
        std::cout << "   • countsketch_heavy_hitters_21mers.csv" << std::endl;
        std::cout << "   • countsketch_heavy_hitters_31mers.csv" << std::endl;
        
        std::cout << "\n✅ Procesamiento completado exitosamente" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}