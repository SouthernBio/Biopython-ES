# Paquete Bio.Align

## 1. Subpaquetes


- Paquete Bio.Align.Applications
	- Contenidos del módulo
- Paquete Bio.Align.substitution_matrices
	- Contenidos del módulo

## 2. Submódulos

- [Bio.Align.AlignInfo]()
- [Bio.Align.a2m]()
- [Bio.Align.bed]()
- [Bio.Align.bigbed]()
- [Bio.Align bigmaf]()
- [Bio.Align.bigpsl]()
- [Bio.Align.clustal]()
- [Bio.Align.emboss]()
- [Bio.Align.exonerate]()
- [Bio.Align.fasta]()
- [Bio.Align.hhr]()
- [Bio.Align.interfaces]()
- [Bio.Align.maf]()
- [Bio.Align.mauve]()
- [Bio.Align.msf]()
- [Bio.Align.nexus]()
- [Bio.Align.phylip]()
- [Bio.Align.psl]()
- [Bio.Align.sam]()
- [Bio.Align.stockholm]()
- [Bio.Align.tabular]()

## 3. Contenidos del módulo

Código para tratar con alineamientos de secuencias.

Uno de los componentes más importantes en este módulo es la clase MultipleSeqAlignment, usada en el módulo Bio.AlignIO.

### `class Bio.Align.AlignmentCounts(gaps, identities, mismatches)`

Base: `tuple`

`__getnewargs__()`

Retorna `self` como una tupla simple. Usado por "copy and pickle".

`static __new__(cls, gaps, identities, mismatches)`

Crea una nueva instancia de AlignmentCounts(gaps, identities, mismatches)

`__repr__()`

Retorna un string de representación formateado.

`__slots__()`

`gaps`

Alias para número de campo 0.

`identities`

Alias para número de campo 1.

`mismatches`

Alias para número de campo 2.


### `class Bio.Align.MultipleSeqAlignment(records, alphabet=None, annotations=None, column_annotations=None)`

Base: `object`

Representa un alineamiento de secuencias múltiple (MSA, en inglés).

Nos referimos con esto a una colección de secuencias (usualmente vista como filas) equilongitudinales (usualmente poseen caracteres "gap" para las inserciones o padding). Los datos pueden ser luego tomados como una matriz de letras, con columnas bien definidas.

Uno crea usualmente un MSA mediante la carga de un archivo de alineamiento con el módulo AlignIO:

```py
>>> from Bio import AlignIO
>>> align = AlignIO.read("Clustalw/opuntia.aln", "clustal")
>>> print(align)
Alignment with 7 rows and 156 columns
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273285|gb|AF191659.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273284|gb|AF191658.1|AF191
TATACATTAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273287|gb|AF191661.1|AF191
TATACATAAAAGAAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273286|gb|AF191660.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273290|gb|AF191664.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273289|gb|AF191663.1|AF191
TATACATTAAAGGAGGGGGATGCGGATAAATGGAAAGGCGAAAG...AGA gi|6273291|gb|AF191665.1|AF191
```

En algunos sentidos se pueden tratar estos objetos como listas de objetos SeqRecord, donde cada uno representa una fila del alineamiento. El iterar sobre un alineamiento devuelve el objeto SeqRecord para cada fila:

```py
>>> len(align)
7
>>> for record in align:
...     print("%s %i" % (record.id, len(record)))
...
gi|6273285|gb|AF191659.1|AF191 156
gi|6273284|gb|AF191658.1|AF191 156
gi|6273287|gb|AF191661.1|AF191 156
gi|6273286|gb|AF191660.1|AF191 156
gi|6273290|gb|AF191664.1|AF191 156
gi|6273289|gb|AF191663.1|AF191 156
gi|6273291|gb|AF191665.1|AF191 156
```

También se puede acceder a filas individuales como objetos SeqRecord a través de sus índices:

```py
>>> print(align[0].id)
gi|6273285|gb|AF191659.1|AF191
>>> print(align[-1].id)
gi|6273291|gb|AF191665.1|AF191
```

Y extraer columnas como strings:

```py
>>> print(align[:, 1])
AAAAAAA
```

O tomar solamente las primeras diez columnas como sub-alineamientos:

```py
>>> print(align[:, :10])
Alignment with 7 rows and 10 columns
TATACATTAA gi|6273285|gb|AF191659.1|AF191
TATACATTAA gi|6273284|gb|AF191658.1|AF191
TATACATTAA gi|6273287|gb|AF191661.1|AF191
TATACATAAA gi|6273286|gb|AF191660.1|AF191
TATACATTAA gi|6273290|gb|AF191664.1|AF191
TATACATTAA gi|6273289|gb|AF191663.1|AF191
TATACATTAA gi|6273291|gb|AF191665.1|AF191
```

Combinando este slicing de alineamientos con adición de alineamientos se puede remover una sección del mismo. Por ejemplo, se puede tomar solamente la primeras y últimas diez columnas:

```py
>>> print(align[:, :10] + align[:, -10:])
Alignment with 7 rows and 20 columns
TATACATTAAGTGTACCAGA gi|6273285|gb|AF191659.1|AF191
TATACATTAAGTGTACCAGA gi|6273284|gb|AF191658.1|AF191
TATACATTAAGTGTACCAGA gi|6273287|gb|AF191661.1|AF191
TATACATAAAGTGTACCAGA gi|6273286|gb|AF191660.1|AF191
TATACATTAAGTGTACCAGA gi|6273290|gb|AF191664.1|AF191
TATACATTAAGTATACCAGA gi|6273289|gb|AF191663.1|AF191
TATACATTAAGTGTACCAGA gi|6273291|gb|AF191665.1|AF191
```

Nota - Este objeto NO intenta modelar el tipo de alineamientos usados en NGS con múltiples lecturas de secuencias que son más cortas que el alineamiento, y en donde hay usualmente una secuencia consenso o de referencia con un estatus especial.

#### `__init__(records, alphabet=None, column_annotations=None, annotations=None)`

Inicializar un nuevo objeto MultipleSeqAlignment.

**Argumentos:**

- **records - Una lista (o iterador) de objetos SeqRecord** cutas secuencias son todas equilongitudinales. 
- **alphabet** - Solo para retrocompatibilidad; su valor debe ser siempre None.
- **annotations** - Información sobre todo el alineamiento (diccionario).
- **column_annotations** - Annotation por columna (diccionario restringido). Esto contiene secuencias de Python (listas, strings, tuplas) cuya longitud es igual al número de columnas. Un uso típico sería usualmente un string consenso de una estructura secundaria.

Normalmente uno carga un MSA desde un archivo usando Bio.AlignIO, pero también se puede hacer desde una lista de objetos SeqRecord:

```py
>>> from Bio.Seq import Seq
>>> from Bio.SeqRecord import SeqRecord
>>> from Bio.Align import MultipleSeqAlignment
>>> a = SeqRecord(Seq("AAAACGT"), id="Alpha")
>>> b = SeqRecord(Seq("AAA-CGT"), id="Beta")
>>> c = SeqRecord(Seq("AAAAGGT"), id="Gamma")
>>> align = MultipleSeqAlignment([a, b, c],
...                             annotations={"tool": "demo"},
...                             column_annotations={"stats": "CCCXCCC"})
>>> print(align)
Alignment with 3 rows and 7 columns
AAAACGT Alpha
AAA-CGT Beta
AAAAGGT Gamma
>>> align.annotations
{'tool': 'demo'}
>>> align.column_annotations
{'stats': 'CCCXCCC'}
```

#### `property column_annotations`

Diccionario de "annotations por letra" para la secuencia.

#### `__str__()`

Retorna un string resumen multilínea del alineamiento.

Este output es legible pero los alineamientos largos se muestran truncados. Se muestra un máximo de 20 filas (de secuencias) y 50 columnas, con los identificadores de registro, por lo que debería poder visualizarse en una sola pantalla. Por ejemplo:

```py
>>> from Bio.Seq import Seq
>>> from Bio.SeqRecord import SeqRecord
>>> from Bio.Align import MultipleSeqAlignment
>>> a = SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha")
>>> b = SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta")
>>> c = SeqRecord(Seq("ACTGCTAGATAG"), id="Gamma")
>>> align = MultipleSeqAlignment([a, b, c])
>>> print(align)
Alignment with 3 rows and 12 columns
ACTGCTAGCTAG Alpha
ACT-CTAGCTAG Beta
ACTGCTAGATAG Gamma
```

Ver también el método de formato de alineamiento.

#### `__repr__()`

Retorna una representación del objeto para depuración.

Esta representación no puede ser usada con `eval()` para recrear el objeto, lo cual es usualmente posible con objetos simples de Python. Por ejemplo:

`Bio.Align.MultipleSeqAlignment instance (2 records of length 14) at a3c184c`

El hex string es la dirección de memoria del objeto, ver `help(id)`. Esto provee una forma simple de distinguir visualmente alineamientos del mismo tamaño.

#### `__format__(format_spec)`

Retorna el alineamiento como un string en el formato especificado.

El formato debería ser un string en minúsculas soportado como un formato de output por Bio.AlignIO (como "fasta", "clustal", "phylip", "stockholm", etcétera), que es usado para convertir el alineamiento en un string. Por ejemplo:

```py
>>> from Bio.Seq import Seq
>>> from Bio.SeqRecord import SeqRecord
>>> from Bio.Align import MultipleSeqAlignment
>>> a = SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha", description="")
>>> b = SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta", description="")
>>> c = SeqRecord(Seq("ACTGCTAGATAG"), id="Gamma", description="")
>>> align = MultipleSeqAlignment([a, b, c])
>>> print(format(align, "fasta"))
>Alpha
ACTGCTAGCTAG
>Beta
ACT-CTAGCTAG
>Gamma
ACTGCTAGATAG

>>> print(format(align, "phylip"))
 3 12
Alpha      ACTGCTAGCT AG
Beta       ACT-CTAGCT AG
Gamma      ACTGCTAGAT AG
```

#### `__iter__()`

Iterar sobre las filas de un alineamiento como objetos SeqRecord. Por ejemplo:

```py
>>> from Bio.Seq import Seq
>>> from Bio.SeqRecord import SeqRecord
>>> from Bio.Align import MultipleSeqAlignment
>>> a = SeqRecord(Seq("ACTGCTAGCTAG"), id="Alpha")
>>> b = SeqRecord(Seq("ACT-CTAGCTAG"), id="Beta")
>>> c = SeqRecord(Seq("ACTGCTAGATAG"), id="Gamma")
>>> align = MultipleSeqAlignment([a, b, c])
>>> for record in align:
...    print(record.id)
...    print(record.seq)
...
Alpha
ACTGCTAGCTAG
Beta
ACT-CTAGCTAG
Gamma
ACTGCTAGATAG
```

#### `__len__()`

Retornar el número de secuencias en el alineamiento.

Use `len(alignment)` para obtener el número de secuencias (es decir, el número de filas), y `alignment.get_alignment_length()` para obtener el largo de la secuencia más larga (lo que equivale al número de columnas).

Esto es fácil de recordar si se visualiza al alineamiento como una lista de objetos SeqRecord.

#### `get_alignment_length()`

Retorna el largo máximo del alineamiento.

Todos los objetos en el alineamiento deberían (supuestamente) tener el mismo largo. Esta función encuentra este largo a través de encontrar el largo máximo de las secuencias en este alineamiento.

