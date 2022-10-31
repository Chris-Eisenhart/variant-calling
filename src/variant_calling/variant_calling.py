import csv
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Set, Tuple


# Frozen so we can use it as a dict key
@dataclass(order=True, frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str

    def to_str(self):
        return f"{self.chrom}-{self.pos}-{self.ref}-{self.alt}"


@dataclass(order=True, frozen=True)
class GenomicPosition:
    chrom: str
    pos: int

    def to_str(self):
        return f"{self.chrom}-{self.pos}"


def get_reverse_complement(dna):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join([complement[base] for base in dna[::-1]])


def parse_md_string(md_string: str) -> List[Tuple[int, str]]:
    parsed_md_list = []
    for matched_region in re.finditer(r"([0-9]+)([\^A,C,T,G]+)", md_string[5:]):
        base_count = int(matched_region.group(1))
        ref_base = matched_region.group(2)
        parsed_md_list.append((base_count, ref_base))
    final_base_count = ""
    char = md_string[-1]
    count = -1
    string_numbers = [str(x) for x in range(0, 10)]
    while char in string_numbers:
        final_base_count += char
        count -= 1
        char = md_string[count]
    parsed_md_list.append((int(final_base_count), "N"))
    return parsed_md_list


def parse_sam_flag(flag: int) -> bool:
    if flag < 16:
        return False
    else:
        # Parse the binary flag, the 5th bit indicates reverse complement
        return bool(int(bin(flag)[-5]))


def identify_and_validate_reference_bases(
    parsed_md_list: List[Tuple[int, str]], seq: str, qual: str, reverse_complement: bool
) -> List[Tuple[int, str, str]]:
    read_variant_list = []  # Parsed variant information from the reads perspective
    for md_element in parsed_md_list:
        base_count, ref_base = md_element
        if ref_base == "N":
            continue
        alt_base = seq[base_count]
        base_quality = ord(qual[base_count]) - 33
        if base_quality <= 20:
            continue
        if reverse_complement:
            alt_base = get_reverse_complement(alt_base)
        read_variant_list.append((base_count, ref_base, alt_base))
    return read_variant_list


def variant_call(
    chrom: str,
    alignment_start_pos: int,
    reverse_complement: bool,
    read_variant_list: List[Tuple[int, str, str]],
) -> Set[Variant]:
    variant_set: Set[Variant] = set()
    total_bases_seen = 0
    for read_variant in read_variant_list:
        base_count, ref_base, alt_base = read_variant
        total_bases_seen += base_count
        if ref_base == "N":
            continue
        if reverse_complement:
            variant_position = alignment_start_pos - total_bases_seen
            variant = Variant(chrom, variant_position, ref_base, alt_base)
        else:
            variant_position = alignment_start_pos + total_bases_seen
            variant = Variant(chrom, variant_position, ref_base, alt_base)
        variant_set.add(variant)
    return variant_set


def variant_calling_for_one_sam_record(
    sam_record: Dict[str, Any]
) -> Optional[Set[Variant]]:
    """
    Identifies all variants in a single BWA alignment record. Returns a list of variant objects.
    """
    # If there are no variants in the alignment MD then return asap
    no_variant = True
    for base in ["A", "T", "C", "G"]:
        if base in sam_record["MD"]:
            no_variant = False
    if no_variant:
        return None

    # Gather the necessary variables from the alignment
    chrom: str = sam_record["RNAME"]
    alignment_start_pos: int = int(sam_record["POS"])
    reverse_complement: bool = parse_sam_flag(int(sam_record["FLAG"]))
    seq: str = sam_record["SEQ"]
    qual: str = sam_record["QUAL"]

    # Call variants in this alignment
    parsed_md_list: List[Tuple[int, str]] = parse_md_string(sam_record["MD"])
    read_variant_list: List[
        Tuple[int, str, str]
    ] = identify_and_validate_reference_bases(
        parsed_md_list, seq, qual, reverse_complement
    )
    variant_set = variant_call(
        chrom, alignment_start_pos, reverse_complement, read_variant_list
    )
    return variant_set


def parse_cigar_string(cigar_string: str) -> List[Tuple[int, str]]:
    parsed_cigar_list: List[Tuple[int, str]] = []
    for matched_region in re.finditer("([0-9]+)([M,I,D,N,S,H,P,=,X])", cigar_string):
        base_count = int(matched_region.group(1))
        cigar_op = matched_region.group(2)
        parsed_cigar_list.append((base_count, cigar_op))
    return parsed_cigar_list


def get_coverage_data_for_one_sam_record(
    sam_record: Dict[str, Any]
) -> Set[GenomicPosition]:
    """
    Determine the valid reference bases that are covered by this alignment where a valid base
    has a quality greater than or equal to 20 phred.

    BWA alignments are done from the reads perspective, every base is considered and reported.
    """
    # Gather the necessary variables from the alignment
    alignment_start_pos: int = int(sam_record["POS"])
    reverse_complement: bool = parse_sam_flag(int(sam_record["FLAG"]))
    qual: str = sam_record["QUAL"]
    chrom: str = sam_record["RNAME"]

    # Process the string cigar into a machine readable list of tuples
    parsed_cigar_list: List[Tuple[int, str]] = parse_cigar_string(sam_record["CIGAR"])
    position_set: Set[GenomicPosition] = set()  # The reference positions which are covered
    current_position_in_ref = (
        alignment_start_pos  # The current position in the reference genome
    )
    current_position_in_read = 0  # The current position in the read

    # Go over the elements of the cigar string, update the current read position
    # as well as the current reference position each iteration so they are always accurate
    for parsed_cigar in parsed_cigar_list:
        base_count, cigar_op = parsed_cigar
        if cigar_op == "M":  # Matched sequence
            start = current_position_in_read
            stop = current_position_in_read + base_count
            for char in qual[start:stop]:
                if ord(char) >= 20:
                    position_set.add(GenomicPosition(chrom, current_position_in_ref))
                if reverse_complement:
                    current_position_in_ref -= 1
                else:
                    current_position_in_ref += 1
        elif cigar_op == "I":  # Insertions
            current_position_in_read += base_count
        elif (
            cigar_op == "D" or cigar_op == "S" or cigar_op == "H"
        ):  # Skipped or deleted sequence
            current_position_in_read += base_count
            if reverse_complement:
                current_position_in_ref -= base_count
            else:
                current_position_in_ref += base_count
    return position_set


def evaluate_sam_file(
    sam_file: str
) -> Tuple[Dict[Variant, int], Dict[GenomicPosition, int]]:
    """
    Evaluate a list of sam records, report two dicts, one that lists all valid variants and their read depth, another
    that lists all valid positions and their read depth.
    """
    variant_to_read_depth: Dict[Variant, int] = {}
    position_to_read_depth: Dict[GenomicPosition, int] = {}
    sam_header = [
        "QNAME",
        "FLAG",
        "RNAME",
        "POS",
        "MAPQ",
        "CIGAR",
        "RNEXT",
        "PNEXT",
        "TLEN",
        "SEQ",
        "QUAL",
        "NM",
        "MD",
        "AS",
        "XS",
        "XA",
    ]
    with open(sam_file, "r") as in_file:
        while in_file.readline().startswith("@"):
            continue # Skip headers
        for sam_record in csv.DictReader(in_file, delimiter="\t", fieldnames=sam_header):
            position_set = get_coverage_data_for_one_sam_record(sam_record)
            for position in position_set:
                if position_to_read_depth.get(position):
                    position_to_read_depth[position] += 1
                else:
                    position_to_read_depth[position] = 1
            variant_set = variant_calling_for_one_sam_record(sam_record)
            if variant_set:
                for variant in variant_set:
                    v_gen_pos = GenomicPosition(variant.chrom, variant.pos)
                    if v_gen_pos in position_set:  # The base has passing quality
                        if variant_to_read_depth.get(variant):
                            variant_to_read_depth[variant] += 1
                        else:
                            variant_to_read_depth[variant] = 1
    return variant_to_read_depth, position_to_read_depth


def write_variant_out_file(
    variant_to_read_depth: Dict[Variant, int],
    position_to_read_depth: Dict[GenomicPosition, int],
    variant_out_file: str,
):
    with open(variant_out_file, "w") as out_file:
        out_file.write("variant\tvar_read_depth\tfull_read_depth\n")
        for variant, variant_read_depth in variant_to_read_depth.items():
            read_depth = position_to_read_depth.get(
                GenomicPosition(variant.chrom, variant.pos), 0
            )
            out_file.write(f"{variant.to_str()}\t{variant_read_depth}\t{read_depth}\n")


def write_position_depth_out_file(
    position_to_read_depth: Dict[GenomicPosition, int], position_out_file: str
):
    with open(position_out_file, "w") as out_file:
        out_file.write("position\tread_depth\n")
        for position, read_depth in position_to_read_depth.items():
            out_file.write(f"{position.to_str()}\t{read_depth}\n")


def call_variants_on_sam_file(
    sam_file: str, variant_out_file: str, position_out_file: str
):
    """
    Stream through the sam alignment file to gather variant and coverage information.
    Write that information to file.
    """
    variant_to_read_depth, position_to_read_depth = evaluate_sam_file(sam_file)
    write_variant_out_file(
        variant_to_read_depth, position_to_read_depth, variant_out_file
    )
    write_position_depth_out_file(position_to_read_depth, position_out_file)
    return
