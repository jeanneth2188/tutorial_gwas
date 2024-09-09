import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Leer el archivo de datos
df = pd.read_csv('results_pyseer_RIF_continuo.tsv', sep='\t')

# Calcular -log10(p-valor)
df['-log10(P-value)'] = -np.log10(df['P-value'])

# Definir el umbral de significancia
umbral_significancia = -np.log10(0.00001)

# Filtrar SNPs significativos
significativos = df[df['-log10(P-value)'] > umbral_significancia]

# Crear una figura
plt.figure(figsize=(12, 6))

# Graficar todos los SNPs
plt.scatter(df['Position'], df['-log10(P-value)'], s=10, color='grey', label='Not Significant')

# Graficar SNPs significativos en rojo
plt.scatter(significativos['Position'], significativos['-log10(P-value)'], s=10, color='red', label='Significant')

# Añadir etiquetas a los SNPs significativos
for i, row in significativos.iterrows():
    plt.text(row['Position'], row['-log10(P-value)'], f'{row["Position"]}', fontsize=8, color='red', ha='right')

# Añadir líneas horizontales para los umbrales de significancia
plt.axhline(y=-np.log10(0.00001), color='blue', linestyle='--')

# Etiquetas y título
plt.xlabel('Genomic Position')
plt.ylabel('-log10(p-value)')
plt.title('Manhattan Plot')
plt.legend()
plt.savefig('Manhattan_Plot.png', dpi=300)

plt.show()
