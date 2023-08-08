# Módulo Bio.Affy.Celfile

Lee información de archivos CEL de Affymetrix (versiones 3 y 4)

### 1. `exception` `Bio.Affy.CelFile.ParserError(*args)`

Base: `ValueError`

Error del parser de Affymetrix.

#### `__init__(*args)`

Inicializa la clase.

### 2. `class` `Bio.Affy.Celfile.Record`

Base: `object`

Almacena la información en un archivo cel.

Ejemplo de uso:

```py
>>> from Bio.Affy import CelFile
>>> with open("Affy/affy_v3_example.CEL") as handle:
...    c = CelFile.read(handle)
...
>>> print(c.ncols, c.nrows)
5 5
>>> print(c.intensities)
[[   234.    170.  22177.    164.  22104.]
 [   188.    188.  21871.    168.  21883.]
 [   188.    193.  21455.    198.  21300.]
 [   188.    182.  21438.    188.  20945.]
 [   193.  20370.    174.  20605.    168.]]
>>> print(c.stdevs)
[[   24.     34.5  2669.     19.7  3661.2]
 [   29.8    29.8  2795.9    67.9  2792.4]
 [   29.8    88.7  2976.5    62.   2914.5]
 [   29.8    76.2  2759.5    49.2  2762. ]
 [   38.8  2611.8    26.6  2810.7    24.1]]
>>> print(c.npix)
[[25 25 25 25 25]
 [25 25 25 25 25]
 [25 25 25 25 25]
 [25 25 25 25 25]
 [25 25 25 25 25]]
```

#### `__init__(*args)`

Inicializa la clase.

### 3. `Bio.Affy.Celfile.read(handle, version=None)`

Lee un archivo CEL de Affymetrix y retorna un objeto Record.

Soporta las versiones 3 y 4 de los archivos de formato CEL. Por favor especifique el formato de archivo CEL como 3 o 4 (de saberlo) para el argumento "version". Si el número de versión no es especificado, el parser intentará detectar la versión a partir de los contenidos del archivo.

El objeto Record retornado por esta función almacena las intensidades del archivo CEL en record.intensities. Actualmente, record.mask y record.outliers no están establecidos al parsear archivos CEL versión 4.

Ejemplo de uso:

```py
>>> from Bio.Affy import CelFile
>>> with open("Affy/affy_v3_example.CEL") as handle:
...     record = CelFile.read(handle)
...
>>> record.version == 3
True
>>> print("%i by %i array" % record.intensities.shape)
5 by 5 array
```

```py
>>> with open("Affy/affy_v4_example.CEL", "rb") as handle:
...     record = CelFile.read(handle, version=4)
...
>>> record.version == 4
True
>>> print("%i by %i array" % record.intensities.shape)
5 by 5 array
```
