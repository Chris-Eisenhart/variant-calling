import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple


@dataclass(order=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str


def get_reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])


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


def identify_and_validate_reference_bases(parsed_md_list: List[Tuple[int, str]], seq: str, qual: str, reverse_complement: bool) -> List[Tuple[int, str, str]]:
    read_variant_list = [] # Parsed variant information from the reads perspective
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
):
    variant_list = []
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
        variant_list.append(variant)
    return variant_list


def variant_calling_for_one_sam_record(
    sam_record: Dict[str, Any]
) -> Optional[List[Variant]]:
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
    seq: str = sam_record['SEQ']
    qual: str = sam_record['QUAL']

    # Call variants in this alignment
    parsed_md_list: List[Tuple[int, str]] = parse_md_string(sam_record["MD"])
    read_variant_list: List[Tuple[int, str, str]] = identify_and_validate_reference_bases(parsed_md_list, seq, qual, reverse_complement)
    variant_list = variant_call(
        chrom, alignment_start_pos, reverse_complement, read_variant_list
    )
    return variant_list


def parse_cigar_string(cigar_string: str) -> List[Tuple[int, str]]:
    parsed_cigar_list: List[Tuple[int, str]] = []
    for matched_region in re.finditer("([0-9]+)([M,I,D,N,S,H,P,=,X])", cigar_string):
        base_count = int(matched_region.group(1))
        cigar_op = matched_region.group(2)
        parsed_cigar_list.append((base_count, cigar_op))
    return parsed_cigar_list


def get_coverage_data_for_one_sam_record(sam_record: Dict[str, Any]) -> List[int]:
    """
    Determine the valid reference bases that are covered by this alignment where a valid base
    has a quality greater than or equal to 20 phred. 

    BWA alignments are done from the reads perspective, every base is considered and reported. 
    """
    # Gather the necessary variables from the alignment
    chrom: str = sam_record["RNAME"]
    alignment_start_pos: int = int(sam_record["POS"])
    reverse_complement: bool = parse_sam_flag(int(sam_record["FLAG"]))
    seq: str = sam_record['SEQ']
    qual: str = sam_record['QUAL']

    # Process the string cigar into a machine readable list of tuples
    parsed_cigar_list: List[Tuple[int, str]] = parse_cigar_string(sam_record["CIGAR"])
    position_list: List[int] = [] # The reference positions which are covered
    current_position_in_ref = alignment_start_pos # The current position in the reference genome
    current_position_in_read = 0 # The current position in the read

    # Go over the elements of the cigar string, update the current read position
    # as well as the current reference position each iteration so they are always accurate
    for parsed_cigar in parsed_cigar_list:
        base_count, cigar_op = parsed_cigar
        if cigar_op == "M": # Matched sequence
            start = current_position_in_read
            stop = current_position_in_read + base_count
            if reverse_complement:
                start = current_position_in_read - base_count
                stop = current_position_in_read
            for char in qual[start:stop]:
                if ord(char) >= 20:
                    position_list.append(current_position_in_ref)
                if reverse_complement:
                    current_position_in_ref -= 1
                else:
                    current_position_in_ref += 1
        elif cigar_op == "I": # Insertions
            current_position_in_read += base_count
        elif cigar_op == "D" or cigar_op == "S" or cigar_op == "H": # Skipped or deleted sequence
            current_position_in_read += base_count
            if reverse_complement:
                current_position_in_ref -= base_count
            else:
                current_position_in_ref += base_count
    return position_list


def evaluate_all_sam_records(sam_record_list: List[Dict[str, Any]]):
    """
    """
    position_to_read_depth: Dict[int, int] = {}
    variant_to_read_depth: Dict[Variant, int] = {}
    for sam_record in sam_record_list:
        variant_list = variant_calling_for_one_sam_record(sam_record)
        position_list = get_coverage_data_for_one_sam_record(sam_record)
 
    return
