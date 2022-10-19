import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple


@dataclass(frozen=True, order=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str


string_numbers = [str(x) for x in range(0, 10)]


def variant_call(
    chrom: str,
    alignment_start_pos: int,
    reverse_complement: bool,
    parsed_md_list: List[Tuple[int, str]],
):
    variant_list = []
    for md_element in parsed_md_list:
        base_count, alt_base = md_element
        if alt_base == "N":
            continue
        if reverse_complement:
            variant_position = alignment_start_pos - base_count
            variant = Variant(chrom, variant_position, "N", alt_base)
        else:
            variant_position = alignment_start_pos + base_count
            variant = Variant(chrom, variant_position, "N", alt_base)
        variant_list.append(variant)
    return variant_list


def parse_md_string(md_string: str) -> List[Tuple[int, str]]:
    parsed_md_list = []
    for matched_region in re.finditer(r"([0-9]+)([\^A,C,T,G]+)", md_string[5:]):
        base_count = int(matched_region.group(1))
        alt_base = matched_region.group(2)
        parsed_md_list.append((base_count, alt_base))
    final_base_count = ""
    char = md_string[-1]
    count = -1
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


def variant_calling_for_one_alignment(
    bwa_alignment: Dict[str, Any]
) -> Optional[List[Variant]]:
    # If there are no variants in the alignment MD then return
    no_variant = True
    for base in ["A", "T", "C", "G"]:
        if base in bwa_alignment["MD"]:
            no_variant = False
    if no_variant:
        return None

    # Gather the necessary variables from the alignment
    chrom: str = bwa_alignment["RNAME"]
    alignment_start_pos: int = int(bwa_alignment["POS"])
    reverse_complement: bool = parse_sam_flag(int(bwa_alignment["FLAG"]))
    parsed_md_list: List[Tuple[int, str]] = parse_md_string(bwa_alignment["MD"])

    # Call variants in this alignment
    variant_list = variant_call(
        chrom, alignment_start_pos, reverse_complement, parsed_md_list
    )

    return variant_list
