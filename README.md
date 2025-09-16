# 🧬 TP Bioinformática

Este proyecto contiene los programas y scripts desarrollados para el trabajo práctico de **Bioinformática**.  
Incluye ejercicios de análisis de secuencias genómicas y proteicas, búsquedas BLAST y alineamientos múltiples (MSA).

---

## ⚙️ Configuración del entorno

### Crear entorno virtual
```bash
python3 -m venv venv
```

### Activar entorno

- **Mac/Linux**
  ```bash
  source venv/bin/activate
  ```

- **Windows**
  ```bash
  venv\Scripts\activate.bat
  ```

### Instalar dependencias
```bash
pip install -r requirements.txt
```

---

## ▶️ Ejecución de los programas

### 📌 Ejercicio 1 – Traducción en 6 frames
Convierte un archivo GenBank con la secuencia del gen **HBB** a proteínas en los 6 marcos de lectura posibles.
```bash
python3 ex1.py -i HBB_NM000518_5.gbk -o HBB_AA_6frames.fasta
```
**Salida:** archivo FASTA con las secuencias proteicas.

---

### 📌 Ejercicio 2 – BLAST remoto
Ejecuta un **BLAST remoto contra la base de datos NCBI**, partiendo de una secuencia de entrada.
```bash
python3 ex2.py -i seq_1.fasta -o blast_remote.xml
```
**Salida:** archivo `blast_remote.xml` en formato XML (outfmt 5) con los resultados de BLAST.

---

### 📌 Ejercicio 3 – Top 10 BLAST + Alineamiento múltiple
Este ejercicio integra los resultados de BLAST con alineamientos múltiples de secuencias.

1. Asegurate de tener una consulta válida en formato **FAA** (por ejemplo, la ORF correcta de HBB ~147 aa):
   ```bash
   query.faa
   ```

2. Ejecutar BLAST remoto, recuperar los **10 mejores hits** y alinearlos con **MUSCLE**:
   ```bash
   export NCBI_EMAIL="tu_email@dominio.com" # requerido por Entrez
   python blast_top10_msa.py --query query.faa --remote --outdir out --msa-tool muscle
   ```

3. Si ya tenés un archivo XML de BLAST (outfmt 5), podés reutilizarlo:
   ```bash
   python blast_top10_msa.py --xml blast.xml --query query.faa --outdir out --msa-tool clustalo
   ```

**Salida:** alineamiento múltiple en FASTA/Clustal dentro de la carpeta `out/`.

---

## 👥 Contribuidores

| [![agmiguel-k](https://github.com/agmiguel-k.png?size=60)](https://github.com/agmiguel-k) | [![gpylypchuk](https://github.com/gpylypchuk.png?size=500)](https://github.com/gpylypchuk) | [![thiagosayegh](https://github.com/thiagosayegh.png?size=60)](https://github.com/thiagosayegh) | [![s-gruss](https://github.com/s-gruss.png?size=60)](https://github.com/s-gruss) |
|:---:|:---:|:---:|:---:|
| [agmiguel-k](https://github.com/agmiguel-k) | [gpylypchuk](https://github.com/gpylypchuk) | [thiagosayegh](https://github.com/thiagosayegh) | [s-gruss](https://github.com/s-gruss) |

---
