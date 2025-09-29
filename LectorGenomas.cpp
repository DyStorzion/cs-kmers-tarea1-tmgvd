#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <filesystem>
#include <vector>

/**
 * Clase para leer archivos genómicos en formato FASTA
 * Permite extraer k-mers de forma secuencial avanzando posición por posición
 */
class LectorGenomas {
private:
    std::string genomicData;              // Almacena toda la secuencia genómica concatenada
    size_t currentPosition;               // Posición actual en la secuencia
    std::vector<std::string> fastaFiles;  // Lista de archivos FASTA
    size_t currentFileIndex;              // Índice del archivo actual
    std::string currentFilename;          // Nombre del archivo actual
    std::string genomasDirectory;         // Directorio de genomas

public:
    /**
     * Constructor que carga todos los archivos FASTA de un directorio
     * @param directory Ruta al directorio que contiene archivos FASTA
     */
    LectorGenomas(const std::string& directory = "Genomas") 
        : currentPosition(0), currentFileIndex(0), genomasDirectory(directory) {
        loadFastaDirectory(directory);
        if (!fastaFiles.empty()) {
            loadCurrentFile();
        }
    }

    /**
     * Carga la lista de archivos FASTA de un directorio
     * @param directory Ruta al directorio que contiene archivos FASTA
     */
    void loadFastaDirectory(const std::string& directory) {
        fastaFiles.clear();
        
        if (!std::filesystem::exists(directory)) {
            throw std::runtime_error("El directorio no existe: " + directory);
        }
        
        for (const auto& entry : std::filesystem::directory_iterator(directory)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().string();
                // Verificar que sea un archivo FASTA (extensiones .fna, .fa, .fasta)
                if (filename.size() >= 4 && 
                    (filename.substr(filename.size() - 4) == ".fna" ||
                     filename.substr(filename.size() - 3) == ".fa" ||
                     (filename.size() >= 6 && filename.substr(filename.size() - 6) == ".fasta"))) {
                    fastaFiles.push_back(filename);
                }
            }
        }
        
        if (fastaFiles.empty()) {
            throw std::runtime_error("No se encontraron archivos FASTA en el directorio: " + directory);
        }
        
        std::cout << "Se encontraron " << fastaFiles.size() << " archivos FASTA en " << directory << std::endl;
    }

    /**
     * Carga el archivo actual basado en currentFileIndex
     */
    void loadCurrentFile() {
        if (currentFileIndex >= fastaFiles.size()) {
            throw std::out_of_range("Índice de archivo fuera de rango");
        }
        
        currentFilename = fastaFiles[currentFileIndex];
        loadFastaFile(currentFilename);
        currentPosition = 0; // Reiniciar posición para el nuevo archivo
    }

    /**
     * Carga el contenido del archivo FASTA, omitiendo las líneas de cabecera
     * @param filename Ruta al archivo FASTA
     */
    void loadFastaFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("No se pudo abrir el archivo: " + filename);
        }

        std::string line;
        genomicData.clear();
        
        while (std::getline(file, line)) {
            // Saltamos las líneas que empiezan con '>' (cabeceras FASTA)
            if (!line.empty() && line[0] != '>') {
                // Concatenamos la línea eliminando saltos de línea
                genomicData += line;
            }
        }
        
        file.close();
        
        if (genomicData.empty()) {
            throw std::runtime_error("El archivo no contiene datos genómicos válidos");
        }
        
        std::cout << "Archivo cargado: " << filename << std::endl;
        std::cout << "Longitud total de la secuencia: " << genomicData.length() << " nucleótidos" << std::endl;
    }

    /**
     * Extrae un k-mer de longitud k desde la posición actual
     * Avanza la posición actual en 1 después de cada llamada
     * Si se acaba el archivo actual, automáticamente pasa al siguiente
     * @param k Longitud del k-mer a extraer
     * @return String con el k-mer extraído, o string vacío si no hay más k-mers en ningún archivo
     */
    std::string getNextKmer(int k) {
        if (k <= 0) {
            throw std::invalid_argument("El valor de k debe ser mayor que 0");
        }
        
        // Verificamos si hay suficientes nucleótidos restantes en el archivo actual
        while (currentPosition + k > genomicData.length()) {
            // Si no hay más nucleótidos en este archivo, intentar pasar al siguiente
            if (hasMoreFiles()) {
                std::cout << "Terminando archivo: " << currentFilename << std::endl;
                nextFile();
                std::cout << "Cambiando a archivo: " << currentFilename << std::endl;
            } else {
                return ""; // No hay más k-mers disponibles en ningún archivo
            }
        }
        
        // Extraemos el k-mer desde la posición actual
        std::string kmer = genomicData.substr(currentPosition, k);
        
        // Avanzamos la posición en 1 para la siguiente llamada
        currentPosition++;
        
        return kmer;
    }

    /**
     * Reinicia la posición actual al inicio de la secuencia
     */
    void reset() {
        currentPosition = 0;
    }

    /**
     * Obtiene la posición actual en la secuencia
     * @return Posición actual (base 0)
     */
    size_t getCurrentPosition() const {
        return currentPosition;
    }

    /**
     * Obtiene la longitud total de la secuencia genómica
     * @return Longitud en nucleótidos
     */
    size_t getSequenceLength() const {
        return genomicData.length();
    }

    /**
     * Verifica si hay más k-mers disponibles para extraer
     * Considera tanto el archivo actual como los archivos restantes
     * @param k Longitud del k-mer
     * @return true si hay más k-mers, false en caso contrario
     */
    bool hasMoreKmers(int k) const {
        // Hay k-mers en el archivo actual
        if (currentPosition + k <= genomicData.length()) {
            return true;
        }
        // O hay más archivos disponibles
        return hasMoreFiles();
    }

    /**
     * Obtiene un fragmento de la secuencia sin avanzar la posición
     * @param start Posición de inicio
     * @param length Longitud del fragmento
     * @return String con el fragmento solicitado
     */
    std::string getSequenceFragment(size_t start, size_t length) const {
        if (start + length > genomicData.length()) {
            throw std::out_of_range("El fragmento solicitado excede la longitud de la secuencia");
        }
        return genomicData.substr(start, length);
    }

    /**
     * Obtiene información del archivo cargado
     */
    void printInfo() const {
        std::cout << "=== Información del archivo FASTA ===" << std::endl;
        std::cout << "Archivo actual: " << currentFilename << std::endl;
        std::cout << "Archivo " << (currentFileIndex + 1) << " de " << fastaFiles.size() << std::endl;
        std::cout << "Longitud de la secuencia: " << genomicData.length() << " nucleótidos" << std::endl;
        std::cout << "Posición actual: " << currentPosition << std::endl;
        std::cout << "Primeros 50 nucleótidos: " << genomicData.substr(0, 50) << "..." << std::endl;
    }

    /**
     * Avanza al siguiente archivo FASTA
     * @return true si hay más archivos, false si ya no quedan más
     */
    bool nextFile() {
        if (currentFileIndex + 1 < fastaFiles.size()) {
            currentFileIndex++;
            loadCurrentFile();
            return true;
        }
        return false;
    }

    /**
     * Retrocede al archivo anterior
     * @return true si se pudo retroceder, false si ya estaba en el primero
     */
    bool previousFile() {
        if (currentFileIndex > 0) {
            currentFileIndex--;
            loadCurrentFile();
            return true;
        }
        return false;
    }

    /**
     * Obtiene el número total de archivos FASTA
     */
    size_t getTotalFiles() const {
        return fastaFiles.size();
    }

    /**
     * Obtiene el índice del archivo actual (base 0)
     */
    size_t getCurrentFileIndex() const {
        return currentFileIndex;
    }

    /**
     * Obtiene el nombre del archivo actual
     */
    std::string getCurrentFilename() const {
        return currentFilename;
    }

    /**
     * Verifica si hay más archivos después del actual
     */
    bool hasMoreFiles() const {
        return currentFileIndex + 1 < fastaFiles.size();
    }

    /**
     * Lista todos los archivos FASTA encontrados
     */
    void listFiles() const {
        std::cout << "=== Lista de archivos FASTA ===" << std::endl;
        for (size_t i = 0; i < fastaFiles.size(); i++) {
            std::string indicator = (i == currentFileIndex) ? " <- ACTUAL" : "";
            std::cout << (i + 1) << ". " << fastaFiles[i] << indicator << std::endl;
        }
    }
};

/**
 * Función de ejemplo que demuestra el uso de la clase LectorGenomas con múltiples archivos
 */



/**
 * Función principal de demostración
 */
/* int main() {
    //ejemploDeUso();
    return 0;
} */