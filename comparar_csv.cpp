#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <iomanip>
#include <algorithm>

struct HeavyHitter {
    std::string kmer;
    int frequency;
    int rank;
};

class CompararHeavyHitters {
private:
    std::vector<HeavyHitter> countsketch_hh;
    std::vector<HeavyHitter> groundtruth_hh;
    
public:
    bool cargarCountSketch(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cout << "❌ No se pudo abrir: " << filename << std::endl;
            return false;
        }
        
        std::string linea;
        std::getline(file, linea); // Saltar header
        
        while (std::getline(file, linea)) {
            std::stringstream ss(linea);
            std::string rank_str, kmer, freq_str;
            
            if (std::getline(ss, rank_str, ',') && 
                std::getline(ss, kmer, ',') && 
                std::getline(ss, freq_str, ',')) {
                
                HeavyHitter hh;
                hh.rank = std::stoi(rank_str);
                hh.kmer = kmer;
                hh.frequency = std::stoi(freq_str);
                countsketch_hh.push_back(hh);
            }
        }
        
        file.close();
        std::cout << "✅ CountSketch cargado: " << countsketch_hh.size() << " heavy hitters" << std::endl;
        return true;
    }
    
    bool cargarGroundTruth(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cout << "❌ No se pudo abrir: " << filename << std::endl;
            return false;
        }
        
        std::string linea;
        std::getline(file, linea); // Saltar header
        
        while (std::getline(file, linea)) {
            std::stringstream ss(linea);
            std::string rank_str, kmer, freq_str;
            
            if (std::getline(ss, rank_str, ',') && 
                std::getline(ss, kmer, ',') && 
                std::getline(ss, freq_str, ',')) {
                
                HeavyHitter hh;
                hh.rank = std::stoi(rank_str);
                hh.kmer = kmer;
                hh.frequency = std::stoi(freq_str);
                groundtruth_hh.push_back(hh);
            }
        }
        
        file.close();
        std::cout << "✅ Ground Truth cargado: " << groundtruth_hh.size() << " heavy hitters" << std::endl;
        return true;
    }
    
    void compararResultados() {
        std::cout << "\n╔══════════════════════════════════════════════════════════════╗" << std::endl;
        std::cout << "║                COMPARACIÓN DE RESULTADOS                   ║" << std::endl;
        std::cout << "╚══════════════════════════════════════════════════════════════╝" << std::endl;
        
        // Crear mapas para búsqueda rápida
        std::unordered_map<std::string, int> cs_map, gt_map;
        
        for (const auto& hh : countsketch_hh) {
            cs_map[hh.kmer] = hh.frequency;
        }
        
        for (const auto& hh : groundtruth_hh) {
            gt_map[hh.kmer] = hh.frequency;
        }
        
        // Calcular métricas
        int verdaderos_positivos = 0;
        int falsos_positivos = 0;
        int falsos_negativos = 0;
        
        double sum_error_absoluto = 0.0;
        double sum_error_relativo = 0.0;
        int coincidencias = 0;
        
        std::cout << "\n🎯 COINCIDENCIAS EXACTAS (k-mers presentes en ambos):" << std::endl;
        std::cout << std::setw(25) << "K-mer" << std::setw(12) << "CountSketch" 
                  << std::setw(12) << "Ground Truth" << std::setw(12) << "Error Abs" 
                  << std::setw(12) << "Error Rel%" << std::endl;
        std::cout << std::string(73, '-') << std::endl;
        
        for (const auto& gt_hh : groundtruth_hh) {
            if (cs_map.count(gt_hh.kmer)) {
                verdaderos_positivos++;
                coincidencias++;
                
                int cs_freq = cs_map[gt_hh.kmer];
                int gt_freq = gt_hh.frequency;
                int error_abs = std::abs(cs_freq - gt_freq);
                double error_rel = (gt_freq > 0) ? (double(error_abs) / gt_freq * 100.0) : 0.0;
                
                sum_error_absoluto += error_abs;
                sum_error_relativo += error_rel;
                
                std::cout << std::setw(25) << gt_hh.kmer.substr(0, 15) + "..." 
                         << std::setw(12) << cs_freq
                         << std::setw(12) << gt_freq
                         << std::setw(12) << error_abs
                         << std::setw(11) << std::fixed << std::setprecision(1) << error_rel << "%" 
                         << std::endl;
            } else {
                falsos_negativos++;
            }
        }
        
        // Contar falsos positivos
        for (const auto& cs_hh : countsketch_hh) {
            if (!gt_map.count(cs_hh.kmer)) {
                falsos_positivos++;
            }
        }
        
        if (coincidencias == 0) {
            std::cout << "❌ No hay coincidencias exactas entre CountSketch y Ground Truth" << std::endl;
        }
        
        // Calcular métricas finales
        double precision = (verdaderos_positivos + falsos_positivos > 0) ? 
                          double(verdaderos_positivos) / (verdaderos_positivos + falsos_positivos) : 0.0;
        double recall = (verdaderos_positivos + falsos_negativos > 0) ? 
                       double(verdaderos_positivos) / (verdaderos_positivos + falsos_negativos) : 0.0;
        double f1_score = (precision + recall > 0) ? 2.0 * precision * recall / (precision + recall) : 0.0;
        
        double mae = (coincidencias > 0) ? sum_error_absoluto / coincidencias : 0.0;
        double mre = (coincidencias > 0) ? sum_error_relativo / coincidencias : 0.0;
        
        std::cout << "\n📊 MÉTRICAS DE EVALUACIÓN:" << std::endl;
        std::cout << "   • Verdaderos Positivos: " << std::setw(6) << verdaderos_positivos << std::endl;
        std::cout << "   • Falsos Positivos:     " << std::setw(6) << falsos_positivos << std::endl;
        std::cout << "   • Falsos Negativos:     " << std::setw(6) << falsos_negativos << std::endl;
        
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "   • Precisión:            " << std::setw(8) << precision << std::endl;
        std::cout << "   • Recall:               " << std::setw(8) << recall << std::endl;
        std::cout << "   • F1-Score:             " << std::setw(8) << f1_score << std::endl;
        
        std::cout << "   • Error Abs. Medio:     " << std::setw(8) << mae << std::endl;
        std::cout << "   • Error Rel. Medio:     " << std::setw(7) << mre << "%" << std::endl;
        
        // TOP comparativo
        std::cout << "\n🏆 TOP 5 COMPARATIVO:" << std::endl;
        
        std::cout << "\n🔵 CountSketch (Top 5):" << std::endl;
        for (int i = 0; i < std::min(5, (int)countsketch_hh.size()); i++) {
            std::cout << "   " << std::setw(2) << (i+1) << ". " 
                      << std::setw(4) << countsketch_hh[i].frequency << " - " 
                      << countsketch_hh[i].kmer.substr(0, 15) << "..." << std::endl;
        }
        
        std::cout << "\n🟢 Ground Truth (Top 5):" << std::endl;
        for (int i = 0; i < std::min(5, (int)groundtruth_hh.size()); i++) {
            std::cout << "   " << std::setw(2) << (i+1) << ". " 
                      << std::setw(4) << groundtruth_hh[i].frequency << " - " 
                      << groundtruth_hh[i].kmer.substr(0, 15) << "..." << std::endl;
        }
        
        // Interpretación
        std::cout << "\n📝 INTERPRETACIÓN:" << std::endl;
        if (precision > 0.8) {
            std::cout << "   ✅ Precisión EXCELENTE (>80%)" << std::endl;
        } else if (precision > 0.6) {
            std::cout << "   ⚠️  Precisión BUENA (60-80%)" << std::endl;
        } else {
            std::cout << "   ❌ Precisión BAJA (<60%)" << std::endl;
        }
        
        if (f1_score > 0.8) {
            std::cout << "   ✅ F1-Score EXCELENTE (>80%)" << std::endl;
        } else if (f1_score > 0.6) {
            std::cout << "   ⚠️  F1-Score BUENO (60-80%)" << std::endl;
        } else {
            std::cout << "   ❌ F1-Score BAJO (<60%)" << std::endl;
        }
    }
};

int main() {
    std::cout << "╔══════════════════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║           COMPARADOR RÁPIDO DE HEAVY HITTERS               ║" << std::endl;
    std::cout << "║         CountSketch vs Ground Truth (desde CSV)             ║" << std::endl;
    std::cout << "╚══════════════════════════════════════════════════════════════╝" << std::endl;
    
    CompararHeavyHitters comparador;
    
    // Cargar datos desde CSV
    if (!comparador.cargarCountSketch("countsketch_heavy_hitters.csv")) {
        std::cout << "⚠️  No se pudo cargar CountSketch CSV. Ejecuta calcular_cs.cpp primero." << std::endl;
        return 1;
    }
    
    if (!comparador.cargarGroundTruth("ground_truth_21mers.csv")) {
        std::cout << "⚠️  No se pudo cargar Ground Truth CSV. Ejecuta extraccion_datos.cpp primero." << std::endl;
        return 1;
    }
    
    // Realizar comparación
    comparador.compararResultados();
    
    std::cout << "\n═══════════════════ COMPARACIÓN COMPLETADA ════════════════════" << std::endl;
    
    return 0;
}