#include "utils/LectorGenomas.hpp"
#include <unordered_map>
#include <iostream>
#include <chrono>
#include <algorithm>

// Función mejorada para obtener k-mer canónico usando strings (más simple)
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

int main(){
    try {
        std::cout << "=== Extracción de Ground Truth para Heavy Hitters ===" << std::endl;
        
        // Crear lector de genomas
        LectorGenomas reader("Genomas");
        
        // Mapas para contar frecuencias exactas
        std::unordered_map<std::string, int> k21mers;
        std::unordered_map<std::string, int> k31mers;
        
        // Estadísticas
        long long totalKmers21 = 0;
        long long totalKmers31 = 0;
        
        do {
            // Extraer k-mers de longitud 21
            reader.reset(); // Reiniciar posición para este archivo
            while (reader.hasMoreKmers(21)) {
                std::string kmer = reader.getNextKmer(21);
                if (!kmer.empty()) {
                    bool validKmer = true;
                    for (char c : kmer) {
                        if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
                            validKmer = false;
                            break;
                        }
                    }
                    
                    if (validKmer) {
                        std::string canonical = getCanonicalKmer(kmer);
                        k21mers[canonical]++;
                        totalKmers21++;
                    }
                }
            }
            
            // Extraer k-mers de longitud 31
            reader.reset(); // Reiniciar posición para este archivo
            while (reader.hasMoreKmers(31)) {
                std::string kmer = reader.getNextKmer(31);
                if (!kmer.empty()) {
                    bool validKmer = true;
                    for (char c : kmer) {
                        if (c != 'A' && c != 'T' && c != 'C' && c != 'G') {
                            validKmer = false;
                            break;
                        }
                    }
                    
                    if (validKmer) {
                        std::string canonical = getCanonicalKmer(kmer);
                        k31mers[canonical]++;
                        totalKmers31++;
                    }
                }
            }
        } while (reader.nextFile()); // Procesar TODOS los archivos
        
        std::cout << "\n=== Estadísticas del procesamiento ===" << std::endl;
        std::cout << "Total 21-mers únicos: " << k21mers.size() << std::endl;
        std::cout << "Total 31-mers únicos: " << k31mers.size() << std::endl;
        std::cout << "Total 21-mers procesados: " << totalKmers21 << std::endl;
        std::cout << "Total 31-mers procesados: " << totalKmers31 << std::endl;
        
        //HHϕ = { x ∈ U : f(x) ≥ ϕ N } donde N es el total de k-mers canónicos (incluyendo repetidos)
        
        // Definir ϕ (ej: 10^-3 a 10^-5) 
        double phi_21 = 2e-6;
        double phi_31 = 6e-6;  
        
        // Calcular umbrales: f(e) ≥ ϕN
        int k21mersBoundary = (int)(phi_21 * totalKmers21);
        int k31mersBoundary = (int)(phi_31 * totalKmers31);
        
        // Analizar distribución de frecuencias para entender mejor los datos
        std::vector<int> frequencies21, frequencies31;
        
        for (const auto& kv : k21mers) {
            frequencies21.push_back(kv.second);
        }
        for (const auto& kv : k31mers) {
            frequencies31.push_back(kv.second);
        }
        
        // Ordenar para obtener estadísticas
        std::sort(frequencies21.begin(), frequencies21.end(), std::greater<int>());
        std::sort(frequencies31.begin(), frequencies31.end(), std::greater<int>());
                
        std::cout << "\n=== Configuración de umbrales (Definición formal) ===" << std::endl;
        std::cout << "Fórmula: HHϕ = { x ∈ U : f(x) ≥ ϕ N }" << std::endl;
        std::cout << "21-mers: ϕ = " << phi_21 << ", N = " << totalKmers21 << std::endl;
        std::cout << "         Umbral = " << k21mersBoundary << " (= " << phi_21 << " × " << totalKmers21 << ")" << std::endl;
        std::cout << "31-mers: ϕ = " << phi_31 << ", N = " << totalKmers31 << std::endl;
        std::cout << "         Umbral = " << k31mersBoundary << " (= " << phi_31 << " × " << totalKmers31 << ")" << std::endl;
        
        // Generar listas de heavy hitters exactos
        std::vector<std::pair<std::string, int>> heavyHitters21;
        std::vector<std::pair<std::string, int>> heavyHitters31;
        
        // Identificar 21-mers heavy hitters
        for (const auto& kv : k21mers) {
            if (kv.second >= k21mersBoundary) {
                heavyHitters21.emplace_back(kv.first, kv.second);
            }
        }
        
        // Identificar 31-mers heavy hitters
        for (const auto& kv : k31mers) {
            if (kv.second >= k31mersBoundary) {
                heavyHitters31.emplace_back(kv.first, kv.second);
            }
        }
        
        // Ordenar por frecuencia descendente
        std::sort(heavyHitters21.begin(), heavyHitters21.end(),
                 [](const auto& a, const auto& b) { return a.second > b.second; });
        std::sort(heavyHitters31.begin(), heavyHitters31.end(),
                 [](const auto& a, const auto& b) { return a.second > b.second; });
        
        int threshold31 = (int)(phi_21 * totalKmers31);
        int threshold21 = (int)(phi_31 * totalKmers21);
        int count21 = 0, count31 = 0;

        for (const auto& kv : k21mers) {
            if (kv.second >= threshold21) count21++;
        }
        for (const auto& kv : k31mers) {
            if (kv.second >= threshold31) count31++;
        }
        
        
        // Mostrar resultados para ϕ seleccionado
        std::cout << "\n=== GROUND TRUTH: Heavy Hitters ===" << std::endl;
        std::cout << "21-mers heavy hitters: " << heavyHitters21.size() << std::endl;
        std::cout << "31-mers heavy hitters: " << heavyHitters31.size() << std::endl;
        
        if (!heavyHitters21.empty()) {
            std::cout << "\nTop 10 Heavy Hitters 21-mers:" << std::endl;
            std::cout << "Rank\tFrecuencia\tK-mer" << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            for (size_t i = 0; i < std::min(size_t(10), heavyHitters21.size()); i++) {
                std::cout << (i+1) << "\t" << heavyHitters21[i].second 
                         << "\t\t" << heavyHitters21[i].first << std::endl;
            }
        }
        
        if (!heavyHitters31.empty()) {
            std::cout << "\nTop 10 Heavy Hitters 31-mers:" << std::endl;
            std::cout << "Rank\tFrecuencia\tK-mer" << std::endl;
            std::cout << "----------------------------------------" << std::endl;
            for (size_t i = 0; i < std::min(size_t(10), heavyHitters31.size()); i++) {
                std::cout << (i+1) << "\t" << heavyHitters31[i].second 
                         << "\t\t" << heavyHitters31[i].first << std::endl;
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}