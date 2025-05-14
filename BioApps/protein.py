from collections import Counter
from re import sub

from .DataBiochemistry import (
    AMINOACIDS_AROMATIC,
    AMINOACIDS_NEGATIVE,
    AMINOACIDS_NONPOLAR,
    AMINOACIDS_POLAR,
    AMINOACIDS_POSITIVE,
    AMINOACIDS_TABLE,
    C_TERM_PKA,
    N_TERM_PKA,
    WATER_MASS,
)


class Protein:
    """
    Represents a protein sequence and provides methods for biochemical analysis.
    """

    def __init__(self, sequence: str = "") -> None:
        """
        Initialize a Protein object with a validated amino acid sequence.

        Args:
            sequence: Amino acid sequence using 1-letter codes.
        """
        self.sequence: str = self._fasta_sequence(sequence)
        self.sequence_size: int = len(self.sequence)

    def _fasta_sequence(self, sequence: str) -> str:
        """
        Clean sequence to contain only valid amino acid characters.

        Args:
            sequence: Raw sequence string.

        Returns:
            Uppercase string with valid amino acid letters only.
        """
        return sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence.upper())

    @property
    def count_aminoacids(self) -> Counter:
        """
        Count occurrences of each amino acid.

        Returns:
            A Counter object with counts of each residue.
        """
        return Counter(self.sequence)

    @property
    def category_counts(self) -> dict:
        """
        Categorize amino acids by biochemical properties.

        Returns:
            A dictionary with counts of aromatic, polar, nonpolar, etc.
        """
        counts = self.count_aminoacids
        return {
            'aromatic':
            sum(counts[aa] for aa in AMINOACIDS_AROMATIC),
            'nonpolar':
            sum(counts[aa] for aa in AMINOACIDS_NONPOLAR),
            'polar':
            sum(counts[aa] for aa in AMINOACIDS_POLAR),
            'polar_negative':
            sum(counts[aa] for aa in AMINOACIDS_NEGATIVE),
            'polar_positive':
            sum(counts[aa] for aa in AMINOACIDS_POSITIVE),
            'polar_neutral':
            sum(counts[aa] for aa in AMINOACIDS_POLAR
                if aa not in AMINOACIDS_NEGATIVE + AMINOACIDS_POSITIVE)
        }

    def aromaticity(self,
                    multiply_by: float = 1.0,
                    decimal_places: int = 2) -> float:
        """
        Calculate fraction of aromatic residues in the protein.

        Args:
            multiply_by: Scale factor (e.g., 100 for percentage).
            decimal_places: Decimal places to round result.

        Returns:
            Scaled and rounded aromaticity value.
        """
        aromatic_count: int = self.category_counts['aromatic']
        return round((aromatic_count / self.sequence_size) * multiply_by,
                     decimal_places)

    def charge_at_pH(self, pH: float = 7.0) -> float:
        """
        Calculate net protein charge at a given pH.

        Args:
            pH: Desired pH value for charge estimation.

        Returns:
            Net charge (rounded to 2 decimal places).
        """
        net_charge: float = 0.0
        pH = round(pH, 3)
        counts = self.count_aminoacids

        # Positive residues
        for aa in AMINOACIDS_POSITIVE:
            if aa in counts:
                pKa = AMINOACIDS_TABLE[aa]['pKr']
                if pKa is not None:
                    net_charge += counts[aa] * (1.0 / (1.0 + 10**(pH - pKa)))

        # Negative residues
        for aa in AMINOACIDS_NEGATIVE:
            if aa in counts:
                pKa = AMINOACIDS_TABLE[aa]['pKr']
                if pKa is not None:
                    net_charge -= counts[aa] * (1.0 / (1.0 + 10**(pKa - pH)))

        # Terminal groups
        net_charge += 1.0 / (1.0 + 10**(pH - N_TERM_PKA))
        net_charge -= 1.0 / (1.0 + 10**(C_TERM_PKA - pH))

        result = round(net_charge, 2)
        return abs(result) if result == 0.00 else result

    def composition_ratio(self,
                          multiply_by: float = 1.0,
                          decimal_places: int = 2) -> dict:
        """
        Calculate ratio of each amino acid in the protein.

        Args:
            multiply_by: Scale factor for each ratio.
            decimal_places: Number of decimals for rounding.

        Returns:
            Dictionary of amino acid ratios.
        """
        return {
            aa: round((count / self.sequence_size) * multiply_by,
                      decimal_places)
            for aa, count in self.count_aminoacids.items()
        }

    def extinction_coefficient(self) -> dict:
        """
        Estimate protein extinction coefficient at 280 nm.

        Based on Trp, Tyr, and half of Cys (as disulfide bonds).

        Returns:
            Dictionary with values assuming reduced and oxidized cysteines.
        """
        counts = self.count_aminoacids
        c_count: int = counts.get("C", 0)
        w_count: int = counts.get("W", 0)
        y_count: int = counts.get("Y", 0)

        # Coefficients in M^-1 cm^-1
        c_coeff = 125
        w_coeff = 5500
        y_coeff = 1490

        disulfide_bonds: int = c_count // 2
        coeff_no_c: int = (w_count * w_coeff) + (y_count * y_coeff)
        coeff_with_c: int = coeff_no_c + (disulfide_bonds * c_coeff)

        return {
            "cys_cystines": round(coeff_with_c, 2),
            "cys_reduced": round(coeff_no_c, 2),
        }

    def hydrophobic_index(self) -> float:
        """
        Calculate average hydrophobicity of the protein.

        Returns:
            Hydrophobic index as a float.
        """
        total: float = sum(AMINOACIDS_TABLE[aa]['hydrophobicity']
                           for aa in self.sequence)
        return round(total / self.sequence_size, 2)

    def isoelectric_point(self) -> float:
        """
        Estimate isoelectric point (pI) using binary search.

        Returns:
            pH at which the net charge is approximately zero.
        """
        low: float = 0.0
        high: float = 14.0
        tolerance: float = 0.01

        while high - low > tolerance:
            mid: float = (low + high) / 2.0
            net_charge: float = self.charge_at_pH(mid)
            if abs(net_charge) < tolerance:
                return round(mid, 2)
            if net_charge > 0:
                low = mid
            else:
                high = mid

        return round((low + high) / 2.0, 2)

    def molecular_weight(self) -> float:
        """
        Compute molecular weight of the protein.

        Returns:
            Total weight with water mass removed for each peptide bond.
        """
        total: float = sum(AMINOACIDS_TABLE[aa]['weight']
                           for aa in self.sequence)
        total -= (self.sequence_size - 1) * WATER_MASS
        return round(total, 2)

    def secondary_structures(self) -> dict:
        """
        Estimate relative propensities for secondary structures.

        Returns:
            Dictionary of predicted % for alpha helix, beta sheet, and coil.
        """
        alpha: float = sum(AMINOACIDS_TABLE[aa]['alpha_helix']
                           for aa in self.sequence) / self.sequence_size
        beta: float = sum(AMINOACIDS_TABLE[aa]['beta_sheet']
                          for aa in self.sequence) / self.sequence_size
        coil: float = 1.0 - (alpha + beta)

        return {
            'alpha_helix': round(alpha * 100, 2),
            'beta_sheet': round(beta * 100, 2),
            'coil': round(coil * 100, 2),
        }
