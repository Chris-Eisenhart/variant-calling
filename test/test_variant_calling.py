from variant_calling.variant_calling import (
    Variant,
    get_coverage_data_for_one_sam_record,
    variant_calling_for_one_sam_record,
    evaluate_sam_record_list
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
    """
    Test with three of the same alignments, everything should be tripled.
    """
    sam_record_1 = {
        "QNAME": "Test_read_1",
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
    sam_record_2 = {
        "QNAME": "Test_read_2",
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
    sam_record_3 = {
        "QNAME": "Test_read_3",
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
    sam_record_list = [sam_record_1, sam_record_2, sam_record_3]
    variant_to_read_depth, position_to_read_depth = evaluate_sam_record_list(sam_record_list)
    assert variant_to_read_depth == {Variant(chrom='chr15', pos=102500900, ref='T', alt='G'): 3, Variant(chrom='chr15', pos=102500903, ref='G', alt='T'): 3, Variant(chrom='chr15', pos=102500906, ref='C', alt='T'): 3}
    assert position_to_read_depth == {102500900: 3, 102500901: 3, 102500902: 3, 102500903: 3, 102500904: 3, 102500905: 3, 102500906: 3, 102500907: 3, 102500908: 3, 102500909: 3, 102500910: 3, 102500911: 3, 102500912: 3, 102500913: 3, 102500914: 3, 102500915: 3, 102500916: 3, 102500917: 3, 102500918: 3, 102500919: 3, 102500920: 3, 102500921: 3, 102500922: 3, 102500923: 3, 102500924: 3, 102500925: 3, 102500926: 3, 102500927: 3, 102500928: 3, 102500929: 3, 102500930: 3, 102500931: 3, 102500934: 3, 102500935: 3, 102500936: 3, 102500937: 3, 102500938: 3, 102500939: 3, 102500940: 3, 102500941: 3, 102500942: 3, 102500943: 3, 102500944: 3, 102500945: 3, 102500946: 3, 102500947: 3, 102500948: 3, 102500949: 3}


def test_evaluate_all_sam_records_variable():
    """
    Pass in three sam records that cover different positions and list different variants.  One variant
    is shared in all three. 
    """
    sam_record_1 = {
        "QNAME": "Test_read_1",
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
    sam_record_2 = {
        "QNAME": "Test_read_2",
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
        "MD": "MD:Z:14A22T3C3G30",
        "AS": "AS:i:76",
        "XS": "XS:i:76",
        "XA": "XA:Z:chr1,+30212,76M,0;chr19,+71819,76M,0;chr2,-114340752,76M,0;chr9,+29990,76M,0;chr12,-73401,10M2I64M,2;",
    }
    sam_record_3 = {
        "QNAME": "Test_read_3",
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
        "MD": "MD:Z:36T36",
        "AS": "AS:i:76",
        "XS": "XS:i:76",
        "XA": "XA:Z:chr1,+30212,76M,0;chr19,+71819,76M,0;chr2,-114340752,76M,0;chr9,+29990,76M,0;chr12,-73401,10M2I64M,2;",
    }
    sam_record_list = [sam_record_1, sam_record_2, sam_record_3]
    variant_to_read_depth, position_to_read_depth = evaluate_sam_record_list(sam_record_list)
    assert variant_to_read_depth == {Variant(chrom='chr15', pos=102500900, ref='T', alt='G'): 2, Variant(chrom='chr15', pos=102500903, ref='G', alt='T'): 1, Variant(chrom='chr15', pos=102500906, ref='C', alt='T'): 1, Variant(chrom='chr15', pos=102500903, ref='C', alt='T'): 1, Variant(chrom='chr15', pos=102500906, ref='G', alt='T'): 1, Variant(chrom='chr15', pos=102500914, ref='T', alt='G'): 1}
    assert position_to_read_depth == {102500900: 3, 102500901: 3, 102500902: 3, 102500903: 3, 102500904: 3, 102500905: 3, 102500906: 3, 102500907: 3, 102500908: 3, 102500909: 3, 102500910: 3, 102500911: 3, 102500912: 3, 102500913: 3, 102500914: 3, 102500915: 3, 102500916: 3, 102500917: 3, 102500918: 3, 102500919: 3, 102500920: 3, 102500921: 3, 102500922: 3, 102500923: 3, 102500924: 3, 102500925: 3, 102500926: 3, 102500927: 3, 102500928: 3, 102500929: 3, 102500930: 3, 102500931: 3, 102500934: 3, 102500935: 3, 102500936: 3, 102500937: 3, 102500938: 3, 102500939: 3, 102500940: 3, 102500941: 3, 102500942: 3, 102500943: 3, 102500944: 3, 102500945: 3, 102500946: 3, 102500947: 3, 102500948: 3, 102500949: 3, 102500878: 2, 102500879: 2, 102500880: 2, 102500881: 2, 102500882: 2, 102500883: 2, 102500884: 2, 102500885: 2, 102500886: 2, 102500887: 2, 102500888: 2, 102500889: 2, 102500890: 2, 102500891: 2, 102500892: 2, 102500893: 2, 102500894: 2, 102500895: 2, 102500896: 2, 102500897: 2, 102500898: 2, 102500899: 2, 102500932: 2, 102500933: 2, 102500950: 2, 102500951: 2, 102500952: 2, 102500953: 2}
