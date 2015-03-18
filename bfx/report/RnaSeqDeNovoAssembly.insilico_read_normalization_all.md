### Normalization

Normalization has been performed in order to reduce memory requirement and decrease assembly runtime by reducing the number of reads, using the [Trinity] normalization utility inspired by the Diginorm algorithm [@diginorm]: first, a catalog of k-mers from all reads is created; then, each RNA-Seq fragment (single read or pair of reads) is probabilistically selected based on its k-mer coverage value and the targeted maximum coverage value. Studies detailed in the Nature Protocol article [@trinity_protocol] found that normalization results in full-length reconstruction to an extent approaching that based on the entire read set.

Table: Normalization Metrics

Trimming Surviving $read_type$ Reads #|Normalization Surviving $read_type$ Reads #|%
----:|----:|----:
$normalization_table$

* Trimming Surviving $read_type$ Reads #: number of remaining $read_type$ Reads after the trimming step
* Normalization Surviving $read_type$ Reads #: number of remaining $read_type$ Reads after the normalization step
* %: percentage of Surviving $read_type$ Reads after normalization / Surviving $read_type$ Reads after trimming
