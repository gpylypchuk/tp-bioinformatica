# BLAST

Este proyecto incluye un script que permite ejecutar BLAST tanto de forma remota (NCBI) como local, con una simple opción de línea de comandos.

## Archivos

- `blast.py` - Script principal unificado (remoto + local)
- `setup_local_blast.py` - Script de configuración automática para modo local
- `requirements.txt` - Dependencias de Python

## Requisitos

### 1. BLAST+ (NCBI BLAST Command Line Tools)

**Windows:**
1. Descargar desde: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
2. Extraer el archivo ZIP
3. Agregar la carpeta `bin` al PATH del sistema
4. Reiniciar terminal/IDE

**Linux:**
```bash
sudo apt-get update
sudo apt-get install ncbi-blast+
```

**macOS:**
```bash
brew install blast
```

### 2. Python Dependencies

```bash
pip install -r requirements.txt
```

## Configuración Automática

Ejecuta el script de configuración para descargar e instalar la base de datos Swiss-Prot:

```bash
python setup_local_blast.py
```

Este script:
- Verifica que BLAST+ esté instalado
- Descarga la base de datos Swiss-Prot desde NCBI
- Crea la base de datos BLAST local
- Limpia archivos temporales

## Uso

### BLAST Unificado - Modo Remoto

```bash
python blast.py --mode remote -i tu_archivo.fasta -o resultado.xml
```

### BLAST Unificado - Modo Local

```bash
python blast.py --mode local -i tu_archivo.fasta -o resultado.xml -d path/a/tu/db
```

### Opciones Disponibles

```bash
python blast.py --help
```

**Parámetros Principales:**
- `--mode {remote,local}`: Modo de ejecución (requerido)
- `-i, --input`: Archivo FASTA con secuencias (requerido)
- `-o, --output`: Archivo de salida (default: blast_results.xml)
- `-f, --outfmt`: Formato de salida (5=XML, 6=Tabular)
- `-n, --num_hits`: Número máximo de hits (default: 10)
- `-d, --database`: Base de datos BLAST local (solo modo local)
- `--info`: Mostrar información detallada del modo

### Ejemplos

```bash
# BLAST Remoto - Salida XML
python blast.py --mode remote -i secuencias.fasta -o resultados.xml

# BLAST Local - Salida tabular
python blast.py --mode local -i secuencias.fasta -f 6 -o resultados.tsv

# BLAST Local - Más hits y base de datos personalizada
python blast.py --mode local -i secuencias.fasta -n 50 -d mi_db -o resultados.xml

# Ver información de los modos
python blast.py --mode remote --info
python blast.py --mode local --info
```

## Ventajas de Cada Modo

### 🌐 Modo Remoto
✅ **Sin instalación**: No requiere BLAST+ local
✅ **Siempre actualizado**: Usa las últimas bases de datos de NCBI
✅ **Sin configuración**: Listo para usar inmediatamente
✅ **Acceso completo**: Todas las bases de datos de NCBI disponibles

### 💻 Modo Local
✅ **Velocidad**: Mucho más rápido que BLAST remoto
✅ **Sin límites**: No hay restricciones de NCBI
✅ **Privacidad**: Los datos no salen de tu computadora
✅ **Control**: Configuración completa de parámetros
✅ **Offline**: Funciona sin conexión a internet

## Solución de Problemas

### Error: "blastp no encontrado"
- Verifica que BLAST+ esté instalado
- Asegúrate de que esté en el PATH del sistema
- Reinicia tu terminal/IDE

### Error: "Base de datos no encontrada"
- Ejecuta `python setup_local_blast.py` para configurar la base de datos
- Verifica que los archivos `.phr`, `.pin`, `.psq`, `.pdb`  existan

### Error de permisos
- En Windows: Ejecuta como administrador
- En Linux/macOS: Usa `sudo` si es necesario

## Comparación: Modo Local vs Remoto

| Aspecto | Modo Local | Modo Remoto |
|---------|------------|-------------|
| Velocidad | ⚡ Muy rápido | 🐌 Lento |
| Límites | ❌ Ninguno | ⚠️ Restricciones NCBI |
| Privacidad | 🔒 Total | 🌐 Datos en NCBI |
| Configuración | 🔧 Completa | 📋 Limitada |
| Instalación | 📦 Requiere BLAST+ | 🚀 Listo para usar |
| Conexión | ❌ No necesaria | 🌐 Requerida |
| Bases de datos | 📁 Locales | 🌍 Todas las de NCBI |

## Guía de Inicio Rápido

### Para usar inmediatamente (Modo Remoto):
```bash
pip install biopython
python blast.py --mode remote -i secuencias.fasta -o resultados.xml
```

### Para máxima velocidad (Modo Local):
```bash
pip install biopython
python setup_local_blast.py  # Configurar base de datos
python blast_unified.py --mode local -i secuencias.fasta -o resultados.xml
```
