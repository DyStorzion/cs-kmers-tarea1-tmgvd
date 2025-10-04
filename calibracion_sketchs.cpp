#include <iostream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <cmath>
#include <algorithm>
#include "sketchs/towersketch.hpp"
#include "sketchs/countsketch.hpp"
#include "utils/LectorGenomas.hpp"

using namespace std;

struct ResultadosError {
    double mae;
    double mre;
};

template <typename SketchType>
ResultadosError calcularErrores(
    unordered_map<string,int>& groundTruth,
    SketchType& sketch
) {
    double mae = 0.0, mre = 0.0;
    int n = 0;

    for (auto& kv : groundTruth) {
        const string& kmer = kv.first;
        int real = kv.second;
        int estimado = sketch.estimate(kmer);

        mae += abs(real - estimado);
        if (real > 0) {
            mre += (double)abs(real - estimado) / real;
        }
        n++;
    }

    if (n > 0) {
        mae /= n;
        mre /= n;
    }

    return {mae, mre};
}

std::string revComp(const std::string &s) {
    std::string rc = s;
    for (char &c : rc) {
        switch(c) {
            case 'A': c='T'; break;
            case 'T': c='A'; break;
            case 'C': c='G'; break;
            case 'G': c='C'; break;
        }
    }
    std::reverse(rc.begin(), rc.end());
    return rc;
}

std::string canonical(const std::string &s) {
    std::string rc = revComp(s);
    return std::min(s, rc);
}

int main() {
    int k = 21;
    string dir = "Genomas";

    LectorGenomas lector(dir);
    unordered_map<string,int> groundTruth;

    cout << "Generando ground truth..." << endl;
    while (lector.hasMoreKmers(k)) {
        string kmer = lector.getNextKmer(k);
        if (kmer.empty()) break;
        groundTruth[canonical(kmer)]++;
    }
    cout << "Se cargaron " << groundTruth.size() << " k-mers únicos" << endl;

    std::string archivocsv = "results_calibracion/calibracion_countsketch_" + std::to_string(k) + "mer.csv";
    ofstream out(archivocsv);
    out << "sketch,d,w,tamano,mae,mre\n";

    vector<int> d_vals = {3, 5, 7};
    vector<int> w_vals = {25000, 35000, 40000};

    
    
    for (int d : d_vals) {
        for (int w : w_vals) {
            CountSketch cs(d, w);

            lector.reset();
            while (lector.hasMoreKmers(k)) {
                string kmer = lector.getNextKmer(k);
                if (kmer.empty()) break;
                cs.insert(kmer);
            }

            auto res = calcularErrores(groundTruth, cs);
            int totalSize = d * w * sizeof(int);
            out << "CS," << d << "," << w << "," << totalSize << ","
                << res.mae << "," << res.mre << "\n";

            cout << "[CS] d=" << d << ", w=" << w
                 << " -> MAE=" << res.mae << ", MRE=" << res.mre << endl;
        }
    }
    cout << "Resultados guardados en " << archivocsv << endl;
    
    archivocsv = "results_calibracion/calibracion_towersketch_" + std::to_string(k) + "mer.csv";
    ofstream tsout(archivocsv);
    tsout << "sketch,d,w8,w16,w32,tamano,mae,mre\n";

    std::vector<int> w8_vals = {90000, 100000};
    std::vector<int> w16_vals = {15000, 20000};
    std::vector<int> w32_vals = {2000, 5000};

    for (int d : d_vals) {
        for (int w8 : w8_vals) {
            for (int w16 : w16_vals) {
                for (int w32 : w32_vals) {
                    TowerSketch ts(d, w8, d, w16, d, w32);

                    lector.reset();
                    while (lector.hasMoreKmers(k)) {
                        string kmer = lector.getNextKmer(k);
                        if (kmer.empty()) break;
                        ts.insert(kmer);
                    }

                    auto res = calcularErrores(groundTruth, ts);
                    size_t totalSize = ts.getSize();
                    tsout << "TS," << d << "," << w8 << "," << w16 << "," << w32 << "," << totalSize << ","
                        << res.mae << "," << res.mre << "\n";

                    cout << "[TS] d=" << d << ", w8=" << w8 << ", w16=" << w16 << ", w32=" << w32
                         << " -> MAE=" << res.mae << ", MRE=" << res.mre << endl;
                }
            }
        }
    }
    
   
    out.close();
    cout << "Resultados guardados en " << archivocsv << endl;

    return 0;
}