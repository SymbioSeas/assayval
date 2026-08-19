# assay-design worked example

A tiny, self-contained example that runs the full assay-design workflow
(ortholog identification → sequence extraction) in a few seconds.

The data is a **small subset** of the FMS2026 *Vibrio mediterranei* Panaroo
pangenome used in Smith 2026 — 15 genomes across 3 clades and 79 orthogroups.
It is for demonstration and testing only, **not** the full manuscript dataset.

## Run it

From this directory (with the `assayval` conda environment active):

```bash
bash run_example.sh
```

Outputs are written to `output/` (gitignored; regenerated on each run).

## Input files

| File | Description |
|------|-------------|
| `example_gene_presence_absence.csv` | Panaroo presence/absence matrix: 14 metadata columns + 15 genome columns × 79 orthogroups. Cells hold the locus ID for that genome (empty = gene absent). |
| `isolate_groups/*.txt` | One file per inclusion group (`Vmed_pathogenic`, `Vmed_protective`, `Vmed_intermediate-avirulent`), one genome name per line. Genome names match matrix column headers. |
| `gene_data.csv` | Panaroo `gene_data` slice for the three representative assemblies only (`gff_file`, `annotation_id`, `prot_sequence`, `dna_sequence`, `gene_name`, `description`). |
| `representatives.tsv` | Maps each group to the representative assembly whose sequences are extracted. |

## Expected outputs (`output/`)

For each of the 3 groups:

| File | Contents |
|------|----------|
| `<group>_conserved_orthologs.tsv` | Orthogroups present in ≥90% of the group's genomes, with inclusion/exclusion coverage metrics. |
| `<group>_specific_orthologs.tsv` | The conserved subset also absent from ≥90% of the exclusion genomes (clade-specific candidates). |
| `<group>_conserved_proteins.faa` / `<group>_conserved_genes.fna` | Representative protein / nucleotide sequences for the conserved orthologs. |
| `<group>_specific_proteins.faa` / `<group>_specific_genes.fna` | Representative protein / nucleotide sequences for the group-specific orthologs (the clade-specific candidate targets). |

With the shipped data you should see roughly: `Vmed_pathogenic` 37 conserved /
12 specific, `Vmed_protective` 47 / 17, `Vmed_intermediate-avirulent` 30 / 4,
and non-empty conserved and specific `.faa`/`.fna` files (representative match
rate 92–100%).

To run the tool on your own Panaroo output, see the
[assay-design README](../README.md).
