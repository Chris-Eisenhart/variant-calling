from variant_calling.variant_calling import Variant, variant_calling_for_one_alignment

test_sam_file = "test_data/chr19_small_test.sam"


def test_variant_calling_for_one_alignment_1():
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
    variant_list = variant_calling_for_one_alignment(sam_record)
    assert not variant_list


def test_variant_calling_for_one_alignment_2():
    sam_record = {
        "QNAME": "Test_read_2",
        "FLAG": "0",
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
        "MD": "MD:Z:36T36",
        "AS": "AS:i:76",
        "XS": "XS:i:76",
        "XA": "XA:Z:chr1,+30212,76M,0;chr19,+71819,76M,0;chr2,-114340752,76M,0;chr9,+29990,76M,0;chr12,-73401,10M2I64M,2;",
    }
    variant_list = variant_calling_for_one_alignment(sam_record)
    variant_1 = Variant(chrom="chr15", pos=102500914, ref="N", alt="T")
    assert variant_list == [variant_1]
