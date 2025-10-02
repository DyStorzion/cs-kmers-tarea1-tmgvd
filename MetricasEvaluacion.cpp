#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <cmath>

class MetricasEvaluacion {
private:
    struct Resultado {
        double precision;
        double recall;
        double f1_score;
        int verdaderos_positivos;
        int falsos_positivos;
        int falsos_negativos;
        int total_estimados;
        int total_reales;
    };

public:
    MetricasEvaluacion() {}

    Resultado evaluarHeavyHitters(const std::vector<std::pair<std::string, int>>& estimaciones,
        const std::vector<std::pair<std::string, int>>& valores_reales,
        int umbral_estimado = -1,
        int umbral_real = -1) {

        // Si no se especifican umbrales, usar los valores mínimos de cada conjunto
        if (umbral_estimado == -1 && !estimaciones.empty()) {
            umbral_estimado = (*std::min_element(estimaciones.begin(), estimaciones.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; })).second;
        }
        
        if (umbral_real == -1 && !valores_reales.empty()) {
            umbral_real = (*std::min_element(valores_reales.begin(), valores_reales.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; })).second;
        }

        // Crear conjuntos de heavy hitters estimados y reales
        std::unordered_set<std::string> hh_estimados = crearConjuntoHH(estimaciones, umbral_estimado);
        std::unordered_set<std::string> hh_reales = crearConjuntoHH(valores_reales, umbral_real);

        return calcularMetricas(hh_estimados, hh_reales);
    }

    // Método alternativo que recibe directamente los conjuntos de heavy hitters
    Resultado evaluarConjuntos(const std::unordered_set<std::string>& hh_estimados,
        const std::unordered_set<std::string>& hh_reales) {
        return calcularMetricas(hh_estimados, hh_reales);
    }

    // Método para evaluar con diferentes umbrales y encontrar el óptimo
    void evaluarMultiplesUmbrales(const std::vector<std::pair<std::string, int>>& estimaciones,
        const std::vector<std::pair<std::string, int>>& valores_reales,
        const std::vector<int>& umbrales_estimados,
        const std::vector<int>& umbrales_reales) {
        std::cout << "\n=== EVALUACIÓN CON MÚLTIPLES UMBRALES ===" << std::endl;
        std::cout << std::setw(12) << "Umbral Est" 
                  << std::setw(12) << "Umbral Real" 
                  << std::setw(10) << "Precision" 
                  << std::setw(10) << "Recall" 
                  << std::setw(10) << "F1-Score" 
                  << std::setw(8) << "TP" 
                  << std::setw(8) << "FP" 
                  << std::setw(8) << "FN" << std::endl;
        std::cout << std::string(88, '-') << std::endl;

        double mejor_f1 = 0.0;
        int mejor_umbral_est = 0, mejor_umbral_real = 0;

        for (int umbral_est : umbrales_estimados) {
            for (int umbral_real : umbrales_reales) {
                Resultado resultado = evaluarHeavyHitters(estimaciones, valores_reales, 
                                                        umbral_est, umbral_real);
                
                std::cout << std::setw(12) << umbral_est
                          << std::setw(12) << umbral_real
                          << std::setw(10) << std::fixed << std::setprecision(3) << resultado.precision
                          << std::setw(10) << std::fixed << std::setprecision(3) << resultado.recall
                          << std::setw(10) << std::fixed << std::setprecision(3) << resultado.f1_score
                          << std::setw(8) << resultado.verdaderos_positivos
                          << std::setw(8) << resultado.falsos_positivos
                          << std::setw(8) << resultado.falsos_negativos << std::endl;

                if (resultado.f1_score > mejor_f1) {
                    mejor_f1 = resultado.f1_score;
                    mejor_umbral_est = umbral_est;
                    mejor_umbral_real = umbral_real;
                }
            }
        }

        std::cout << "\nMEJOR CONFIGURACIÓN:" << std::endl;
        std::cout << "Umbral Estimado: " << mejor_umbral_est << std::endl;
        std::cout << "Umbral Real: " << mejor_umbral_real << std::endl;
        std::cout << "F1-Score: " << std::fixed << std::setprecision(3) << mejor_f1 << std::endl;
    }

    // Método para mostrar análisis detallado
    void mostrarAnalisisDetallado(const std::vector<std::pair<std::string, int>>& estimaciones,
        const std::vector<std::pair<std::string, int>>& valores_reales,
        int umbral_estimado,
        int umbral_real) {
        std::unordered_set<std::string> hh_estimados = crearConjuntoHH(estimaciones, umbral_estimado);
        std::unordered_set<std::string> hh_reales = crearConjuntoHH(valores_reales, umbral_real);
        
        Resultado resultado = calcularMetricas(hh_estimados, hh_reales);

        std::cout << "\n=== ANÁLISIS DETALLADO DE HEAVY HITTERS ===" << std::endl;
        std::cout << "Umbral para estimaciones: " << umbral_estimado << std::endl;
        std::cout << "Umbral para valores reales: " << umbral_real << std::endl;
        std::cout << "\nRESULTADOS:" << std::endl;
        std::cout << "Precision: " << std::fixed << std::setprecision(3) << resultado.precision << std::endl;
        std::cout << "Recall: " << std::fixed << std::setprecision(3) << resultado.recall << std::endl;
        std::cout << "F1-Score: " << std::fixed << std::setprecision(3) << resultado.f1_score << std::endl;
        
        std::cout << "\nMÉTRICAS BÁSICAS:" << std::endl;
        std::cout << "Verdaderos Positivos (TP): " << resultado.verdaderos_positivos << std::endl;
        std::cout << "Falsos Positivos (FP): " << resultado.falsos_positivos << std::endl;
        std::cout << "Falsos Negativos (FN): " << resultado.falsos_negativos << std::endl;
        std::cout << "Total HH Estimados: " << resultado.total_estimados << std::endl;
        std::cout << "Total HH Reales: " << resultado.total_reales << std::endl;

        mostrarEjemplos(hh_estimados, hh_reales, estimaciones, valores_reales);
    }

    // Método para comparar precisión de estimaciones individuales
    void evaluarPrecisionEstimaciones(const std::vector<std::pair<std::string, int>>& estimaciones,
        const std::vector<std::pair<std::string, int>>& valores_reales) {
        std::unordered_map<std::string, int> mapa_estimaciones;
        std::unordered_map<std::string, int> mapa_reales;

        for (const auto& par : estimaciones) {
            mapa_estimaciones[par.first] = par.second;
        }
        for (const auto& par : valores_reales) {
            mapa_reales[par.first] = par.second;
        }

        std::cout << "\n=== PRECISIÓN DE ESTIMACIONES INDIVIDUALES ===" << std::endl;
        std::cout << std::setw(25) << "K-mer" 
                  << std::setw(12) << "Estimado" 
                  << std::setw(12) << "Real" 
                  << std::setw(12) << "Error Rel %" << std::endl;
        std::cout << std::string(61, '-') << std::endl;

        double error_promedio = 0.0;
        int conteo = 0;

        for (const auto& par : estimaciones) {
            const std::string& kmer = par.first;
            int estimado = par.second;
            
            if (mapa_reales.find(kmer) != mapa_reales.end()) {
                int real = mapa_reales[kmer];
                double error_relativo = std::abs(estimado - real) * 100.0 / real;
                
                std::cout << std::setw(25) << kmer.substr(0, 20) + "..."
                          << std::setw(12) << estimado
                          << std::setw(12) << real
                          << std::setw(12) << std::fixed << std::setprecision(1) << error_relativo << std::endl;
                
                error_promedio += error_relativo;
                conteo++;
                
                if (conteo >= 10) break; // Mostrar solo los primeros 10
            }
        }

        if (conteo > 0) {
            std::cout << "\nError relativo promedio: " << std::fixed << std::setprecision(2) 
                      << (error_promedio / conteo) << "%" << std::endl;
        }
    }

private:
    // Crear conjunto de heavy hitters basado en umbral
    std::unordered_set<std::string> crearConjuntoHH(const std::vector<std::pair<std::string, int>>& datos,
        int umbral) {
        std::unordered_set<std::string> conjunto;
        for (const auto& par : datos) {
            if (par.second >= umbral) {
                conjunto.insert(par.first);
            }
        }
        return conjunto;
    }

    // Calcular métricas de precisión, recall y F1-score
    Resultado calcularMetricas(const std::unordered_set<std::string>& hh_estimados,
        const std::unordered_set<std::string>& hh_reales) {
        Resultado resultado;

        // Calcular verdaderos positivos
        resultado.verdaderos_positivos = 0;
        for (const std::string& kmer : hh_estimados) {
            if (hh_reales.find(kmer) != hh_reales.end()) {
                resultado.verdaderos_positivos++;
            }
        }

        // Calcular falsos positivos y negativos
        resultado.falsos_positivos = hh_estimados.size() - resultado.verdaderos_positivos;
        resultado.falsos_negativos = hh_reales.size() - resultado.verdaderos_positivos;
        
        resultado.total_estimados = hh_estimados.size();
        resultado.total_reales = hh_reales.size();

        // Calcular precision, recall y F1-score
        if (resultado.total_estimados > 0) {
            resultado.precision = static_cast<double>(resultado.verdaderos_positivos) / resultado.total_estimados;
        } else {
            resultado.precision = 0.0;
        }

        if (resultado.total_reales > 0) {
            resultado.recall = static_cast<double>(resultado.verdaderos_positivos) / resultado.total_reales;
        } else {
            resultado.recall = 0.0;
        }

        if (resultado.precision + resultado.recall > 0) {
            resultado.f1_score = 2.0 * (resultado.precision * resultado.recall) / (resultado.precision + resultado.recall);
        } else {
            resultado.f1_score = 0.0;
        }

        return resultado;
    }

    // Mostrar ejemplos de cada categoría
    void mostrarEjemplos(const std::unordered_set<std::string>& hh_estimados,
        const std::unordered_set<std::string>& hh_reales,
        const std::vector<std::pair<std::string, int>>& estimaciones,
        const std::vector<std::pair<std::string, int>>& valores_reales) {
        std::cout << "\nEJEMPLOS:" << std::endl;
        
        std::cout << "Verdaderos Positivos (estimados correctamente como HH):" << std::endl;
        int count = 0;
        for (const std::string& kmer : hh_estimados) {
            if (hh_reales.find(kmer) != hh_reales.end()) {
                std::cout << "  " << kmer.substr(0, 15) << "..." << std::endl;
                if (++count >= 3) break;
            }
        }

        std::cout << "Falsos Positivos (estimados como HH pero no son reales):" << std::endl;
        count = 0;
        for (const std::string& kmer : hh_estimados) {
            if (hh_reales.find(kmer) == hh_reales.end()) {
                std::cout << "  " << kmer.substr(0, 15) << "..." << std::endl;
                if (++count >= 3) break;
            }
        }

        std::cout << "Falsos Negativos (HH reales no detectados):" << std::endl;
        count = 0;
        for (const std::string& kmer : hh_reales) {
            if (hh_estimados.find(kmer) == hh_estimados.end()) {
                std::cout << "  " << kmer.substr(0, 15) << "..." << std::endl;
                if (++count >= 3) break;
            }
        }
    }
};