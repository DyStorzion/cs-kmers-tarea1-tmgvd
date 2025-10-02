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

    ofstream out("calibracion_sketches.csv");
    out << "sketch,d,w,tamano,mae,mre\n";

    vector<int> d_vals = {3, 5, 7};
    vector<int> w_vals = {500, 1000, 2000, 5000, 10000};

    /*
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
    */
    
    for (int d : d_vals) {
        for (int w : w_vals) {
            TowerSketch ts(d, w);

            lector.reset();
            while (lector.hasMoreKmers(k)) {
                string kmer = lector.getNextKmer(k);
                if (kmer.empty()) break;
                ts.insert(kmer);
            }

            auto res = calcularErrores(groundTruth, ts);
            size_t totalSize = ts.getSize();
            out << "TS," << d << "," << w << "," << totalSize << ","
                << res.mae << "," << res.mre << "\n";

            cout << "[TS] d=" << d << ", w=" << w
                 << " -> MAE=" << res.mae << ", MRE=" << res.mre << endl;
        }
    }
    
   
    out.close();
    cout << "Resultados guardados en calibracion_sketches.csv" << endl;

    return 0;
}