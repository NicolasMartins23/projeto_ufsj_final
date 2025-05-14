import re
from collections import Counter

from .DataBiochemistry import (
    CODON_SIZE,
    STOP_CODON_DNA,
    STOP_CODON_RNA,
    TABLE_DNA_CODON_TO_AMINOACID,
    TABLE_RNA_CODON_TO_AMINOACID,
)
from .protein import Protein


class NucleicAcid:
    """
    Abstract base class for nucleic acid sequences (DNA or RNA).
    Handles common sequence-based operations like counting and translation.
    """

    def __init__(self, sequence: str, codon_table: dict, stop_codons: list[str]) -> None:
        """
        Initialize a nucleic acid sequence.

        Args:
            sequence: A string representing the nucleotide sequence.
            codon_table: Dictionary mapping codons to amino acids.
        """
        self.sequence: str = sequence
        self.sequence_size: int = len(self.sequence)
        self.codon_table: dict = codon_table
        self.stop_codons: list[str] = stop_codons

    @property
    def count_nucleotides(self) -> dict[str, int]:
        """
        Count each nucleotide and return total length.

        Returns:
            Dictionary of nucleotide counts with 'total' included.
        """
        counter = Counter(self.sequence)
        return {**counter}

    def _trim_at_stop_codon(self) -> str:
        """
        Truncate the sequence at the first stop codon encountered.

        Returns:
            Sequence up to and including the first stop codon.
        """
        trimmed_sequence: str = ""
        for i in range(0, self.sequence_size, CODON_SIZE):
            codon: str = self.sequence[i:i + CODON_SIZE]
            trimmed_sequence += codon
            if codon in self.stop_codons:
                return trimmed_sequence
        return trimmed_sequence


    def gc_skew(self, decimal_precision: int = 3) -> float:
        """
        Calculate GC skew = (G - C) / (G + C).

        Returns:
            GC skew value scaled and rounded.
        """
        g = self.count_nucleotides.get("G", 0)
        c = self.count_nucleotides.get("C", 0)
        gc_skew: float = (g - c) / (g + c) if (g + c) != 0 else 0
        return round(gc_skew, decimal_precision)

    def gc_content(
        self,
        decimal_precision: int = 3) -> float:
        """
        Calculate GC content as a percentage or ratio.

        Args:
            multiply_by: Scale the result (e.g. 100 for percentage).
            decimal_precision: Number of decimal places to round to.

        Returns:
            GC content as a float.
        """
        gc_count: int = self.count_nucleotides.get("C", 0) + self.count_nucleotides.get("G", 0)
        gc_content: float = (gc_count / self.sequence_size)
        return round(gc_content, decimal_precision)


    def peptide_sequence(
        self,
        show_stop_codon: bool = False,
        ignore_stop_codons: bool = False) -> str:
        """
        Translate nucleotide sequence into a peptide string using codons.

        Args:
            show_stop_codon: If False, the final stop codon is excluded.

        Returns:
            The resulting amino acid sequence as a string.
        """
        peptide_sequence: str = ""

        if ignore_stop_codons:
            for i in range(0, self.sequence_size, CODON_SIZE):
                codon: str = self.sequence[i:i + CODON_SIZE]
                peptide_sequence += self.codon_table.get(codon, "")

        else:
            trim_peptide_sequence: str = self._trim_at_stop_codon()
            for i in range(0, len(trim_peptide_sequence), CODON_SIZE):
                codon: str = trim_peptide_sequence[i:i + CODON_SIZE]
                peptide_sequence += self.codon_table.get(codon, "")

        return peptide_sequence if show_stop_codon else peptide_sequence[:-1]

    def to_protein(self) -> Protein:
        """
        Convert the nucleotide sequence into a Protein object.

        Returns:
            A Protein object with the translated amino acid sequence.
        """
        peptide_sequence: str = self.peptide_sequence()
        return Protein(peptide_sequence)


class DNA(NucleicAcid):
    """
    DNA class extending NucleicAcid with DNA-specific methods.
    """

    def __init__(self, sequence: str) -> None:
        """
        Initialize a DNA object with cleaned sequence and DNA codon table.
        """
        fasta_sequence: str = self._fasta_sequence(sequence)
        super().__init__(fasta_sequence, TABLE_DNA_CODON_TO_AMINOACID, STOP_CODON_DNA)

    def _fasta_sequence(self, dna_seq: str = "") -> str:
        """
        Clean the sequence to contain only valid DNA bases.

        Returns:
            A cleaned uppercase string of A, C, G, T.
        """
        return re.sub(r'[^ACGT]', '', dna_seq.upper())

    def template_strand(self, reverse_string: bool = True) -> str:
        """
        Generate the complementary DNA strand.

        Args:
            reverse_string: Whether to reverse the strand (5'→3').

        Returns:
            The complementary strand as a string.
        """
        comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        template_strand = ''.join(comp.get(nt, '') for nt in self.sequence)
        return template_strand[::-1] if reverse_string else template_strand

    def rna_sequence(self) -> str:
        """
        Convert DNA sequence to RNA (T → U).

        Returns:
            RNA sequence as a string.
        """
        return self.sequence.replace("T", "U")

    def to_rna(self) -> "RNA":
        """
        Convert the DNA object to an RNA object.

        Returns:
            An RNA object with the transcribed sequence.
        """
        return RNA(self.rna_sequence())


class RNA(NucleicAcid):
    """
    RNA class extending NucleicAcid with RNA-specific methods.
    """

    def __init__(self, sequence: str) -> None:
        """
        Initialize an RNA object with cleaned sequence and RNA codon table.
        """
        fasta_sequence: str = self._fasta_sequence(sequence)
        super().__init__(fasta_sequence, TABLE_RNA_CODON_TO_AMINOACID, STOP_CODON_RNA)

    def _fasta_sequence(self, sequence: str = "") -> str:
        """
        Clean the sequence to contain only valid RNA bases.

        Returns:
            A cleaned uppercase string of A, C, G, U.
        """
        return re.sub(r'[^ACGU]', '', sequence.upper())

    def dna_sequence(self) -> str:
        """
        Convert RNA sequence to DNA (U → T).

        Returns:
            DNA sequence as a string.
        """
        return self.sequence.replace("U", "T")

    def to_dna(self) -> DNA:
        """
        Convert the RNA object to a DNA object.

        Returns:
            A DNA object with the reverse-transcribed sequence.
        """
        return DNA(self.dna_sequence())

