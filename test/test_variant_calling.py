from variant_calling.variant_calling import (
    Variant,
    get_coverage_data_for_one_sam_record,
    variant_calling_for_one_sam_record,
)

test_sam_file = "test_data/chr19_small_test.sam"


def test_variant_calling_null():
    """
    Test the null case where there are no variants in the alignment
    """
    sam_record = {
        "QNAME": "SRR1518133.27",
        "FLAG": "16",
        "RNAME": "chr15",
        "POS": "102500878",
        "MAPQ": "0",
        "CIGAR": "76M",
        "RNEXT": "*",
        "PNEXT": "0",
        "TLEN": "0",
        "SEQ": "TGCTGGACTTTGGACTGATGATGCTCTTTTTATATATTTAACTTTTTAAAAAAGCCTCTTTTCTTTCTTTTTACCA",
        "QUAL": "DFEFFFHHHHGIHGIJJJJJIJJJJJJJJJJJIJJJIGJJIIIJJJIJIJJJIIIGHGIIIHHDGHHHFFFFDCB@",
        "NM": "NM:i:0",
        "MD": "MD:Z:76",
        "AS": "AS:i:76",
        "XS": "XS:i:76",
        "XA": "XA:Z:chr1,+30212,76M,0;chr19,+71819,76M,0;chr2,-114340752,76M,0;chr9,+29990,76M,0;chr12,-73401,10M2I64M,2;",
    }
    variant_list = variant_calling_for_one_sam_record(sam_record)
    assert not variant_list


def test_variant_calling_for_one_sam_record():
    """
    Test a basic case, there is a single variant with high quality supporting base
    """
    sam_record = {
        "QNAME": "test_read",
        "FLAG": "0",
        "RNAME": "chr15",
        "POS": "102500878",
        "MAPQ": "0",
        "CIGAR": "76M",
        "RNEXT": "*",
        "PNEXT": "0",
        "TLEN": "0",
        "SEQ": "TGCTGGACTTTGGACTGATGATGCTCTTTTTAAAAAGAAAAACTTTTTAAAAAAGCCTCTTTTCTTTCTTTTTACCA",
        "QUAL": "DFEFFFHHHHGIHGIJJJJJIJJJJJJJJJJJIJJJIGJJIIIJJJIJIJJJIIIGHGIIIHHDGHHHFFFFDCB@",
        "NM": "NM:i:0",
        "MD": "MD:Z:36T36",
        "AS": "AS:i:76",
        "XS": "XS:i:76",
        "XA": "XA:Z:chr1,+30212,76M,0;chr19,+71819,76M,0;chr2,-114340752,76M,0;chr9,+29990,76M,0;chr12,-73401,10M2I64M,2;",
    }
    variant_list = variant_calling_for_one_sam_record(sam_record)
    variant_1 = Variant(chrom="chr15", pos=102500914, ref="T", alt="G")
    assert variant_list == [variant_1]


def test_variant_calling_for_one_sam_record_multi_variants():
    """
    Test an alignment that has multiple variants, some of which do not have high enough supporting quality.

    The first variant at read position 14 should not be called since the quality is '4' which converts to 19
    in the phred scaling.
    """
    sam_record = {
        "QNAME": "test_read",
        "FLAG": "0",
        "RNAME": "chr15",
        "POS": "102500878",
        "MAPQ": "0",
        "CIGAR": "76M",
        "RNEXT": "*",
        "PNEXT": "0",
        "TLEN": "0",
        "SEQ": "TGCTGGACTTTGGACTGATGATGCTCTTTTTAAAAAGAAAAACTTTTTAAAAAAGCCTCTTTTCTTTCTTTTTACCA",
        "QUAL": "DFEFFFHHHHGIHA4AJJJJIJJJJJJJJJJJIJJJIGJJIIIJJJIJIJJJIIIGHGIIIHHDGHHHFFFFDCB@",
        "NM": "NM:i:0",
        "MD": "MD:Z:14G22T3G3C30",
        "AS": "AS:i:76",
        "XS": "XS:i:76",
        "XA": "XA:Z:chr1,+30212,76M,0;chr19,+71819,76M,0;chr2,-114340752,76M,0;chr9,+29990,76M,0;chr12,-73401,10M2I64M,2;",
    }
    variant_list = variant_calling_for_one_sam_record(sam_record)
    variant_1 = Variant(chrom="chr15", pos=102500900, ref="T", alt="G")
    variant_2 = Variant(chrom="chr15", pos=102500903, ref="G", alt="T")
    variant_3 = Variant(chrom="chr15", pos=102500906, ref="C", alt="T")
    assert variant_list == [variant_1, variant_2, variant_3]


def test_get_coverage_data_for_one_sam_record():
    """
    A basic alignment with a single match segment
    """
    sam_record = {
        "QNAME": "test_read",
        "FLAG": "0",
        "RNAME": "chr15",
        "POS": "102500878",
        "MAPQ": "0",
        "CIGAR": "76M",
        "RNEXT": "*",
        "PNEXT": "0",
        "TLEN": "0",
        "SEQ": "TGCTGGACTTTGGACTGATGATGCTCTTTTTAAAAAGAAAAACTTTTTAAAAAAGCCTCTTTTCTTTCTTTTTACCA",
        "QUAL": "DFEFFFHHHHGIHA4AJJJJIJJJJJJJJJJJIJJJIGJJIIIJJJIJIJJJIIIGHGIIIHHDGHHHFFFFDCB@",
        "NM": "NM:i:0",
        "MD": "MD:Z:14G22T3G3C30",
        "AS": "AS:i:76",
        "XS": "XS:i:76",
        "XA": "XA:Z:chr1,+30212,76M,0;chr19,+71819,76M,0;chr2,-114340752,76M,0;chr9,+29990,76M,0;chr12,-73401,10M2I64M,2;",
    }
    position_list = get_coverage_data_for_one_sam_record(sam_record)
    assert position_list == [
        102500878,
        102500879,
        102500880,
        102500881,
        102500882,
        102500883,
        102500884,
        102500885,
        102500886,
        102500887,
        102500888,
        102500889,
        102500890,
        102500891,
        102500892,
        102500893,
        102500894,
        102500895,
        102500896,
        102500897,
        102500898,
        102500899,
        102500900,
        102500901,
        102500902,
        102500903,
        102500904,
        102500905,
        102500906,
        102500907,
        102500908,
        102500909,
        102500910,
        102500911,
        102500912,
        102500913,
        102500914,
        102500915,
        102500916,
        102500917,
        102500918,
        102500919,
        102500920,
        102500921,
        102500922,
        102500923,
        102500924,
        102500925,
        102500926,
        102500927,
        102500928,
        102500929,
        102500930,
        102500931,
        102500932,
        102500933,
        102500934,
        102500935,
        102500936,
        102500937,
        102500938,
        102500939,
        102500940,
        102500941,
        102500942,
        102500943,
        102500944,
        102500945,
        102500946,
        102500947,
        102500948,
        102500949,
        102500950,
        102500951,
        102500952,
        102500953,
    ]


def test_get_coverage_data_for_one_sam_record_complex_cigar():
    """
    A complex alignment with two match segments, soft clipping and a deleted and inserted segment
    """
    sam_record = {
        "QNAME": "Test_read",
        "FLAG": "0",
        "RNAME": "chr15",
        "POS": "102500878",
        "MAPQ": "0",
        "CIGAR": "22S16M3I16M2D16M",
        "RNEXT": "*",
        "PNEXT": "0",
        "TLEN": "0",
        "SEQ": "TGCTGGACTTTGGACTGATGATGCTCTTTTTAAAAAGAAAAACTTTTTAAAAAAGCCTCTTTTCTTTCTTTTTACCA",
        "QUAL": "DFEFFFHHHHGIHA4AJJJJIJJJJJJJJJJJIJJJIGJJIIIJJJIJIJJJIIIGHGIIIHHDGHHHFFFFDCB@",
        "NM": "NM:i:0",
        "MD": "MD:Z:14G22T3G3C30",
        "AS": "AS:i:76",
        "XS": "XS:i:76",
        "XA": "XA:Z:chr1,+30212,76M,0;chr19,+71819,76M,0;chr2,-114340752,76M,0;chr9,+29990,76M,0;chr12,-73401,10M2I64M,2;",
    }
    position_list = get_coverage_data_for_one_sam_record(sam_record)
    assert position_list == [
        102500900,
        102500901,
        102500902,
        102500903,
        102500904,
        102500905,
        102500906,
        102500907,
        102500908,
        102500909,
        102500910,
        102500911,
        102500912,
        102500913,
        102500914,
        102500915,
        102500916,
        102500917,
        102500918,
        102500919,
        102500920,
        102500921,
        102500922,
        102500923,
        102500924,
        102500925,
        102500926,
        102500927,
        102500928,
        102500929,
        102500930,
        102500931,
        102500934,
        102500935,
        102500936,
        102500937,
        102500938,
        102500939,
        102500940,
        102500941,
        102500942,
        102500943,
        102500944,
        102500945,
        102500946,
        102500947,
        102500948,
        102500949,
    ]


def test_evaluate_all_sam_records():
    return 0
