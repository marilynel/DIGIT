import json

# >Ds-GFP_TAG21 (2165 bp)
dsgfp = "TAGGGATGAAAGTAGGATGGGAAAATCCCGTACCGACCGTTATCGTATAACCGATTTTGTTAGTT" + \
        "TTATCCCGATCGATTTCGAACCCGAGGTAAAAAACGAAAACGGAACGGAAACGGGATATACAAAACGGT" + \
        "AAACGGAAACGGAAACGGTAGAGCTAGTTTCCCGACCGTTTCACCGGGATCCCGTTTTTAATCGGAATG" + \
        "ATCCCGTTTCGTTACCGTATTTTCTAATTCGGGATGACTGCAATATGGCCAGCTCCAACTCCCATCCAT" + \
        "AACCACTGAGGCCCAGCCCATGTAAGAAATACCTAGCGAACGCTGCTCTGCCTCTCTCCCAGGCGGCCA" + \
        "GGCACCACACGAGTAACAGCATCACACATTCACACGCCGCCACGCGCCCACGCCGGAGTCCGGACGCCG" + \
        "CCAGCCGCACGCCGACGCCGGCGACGCGTCTCGCTCTCGCCTGCTCTCTCCGACTCTCCCTGTCTCCCA" + \
        "GCCGGCCGGCCGCTGGGCTGCACCAGGCACCACACGCGGTGACGGCCGTGACGCGGCACGCCGGACGCA" + \
        "GACGCCGCCATCCACGGTCCGCCCTCCACTCCACTGCTCGGCTCTAGCGAAGGGGTTCGAGCTTGGagc" + \
        "ttACATGTGTAAAGGTGAAGAGATCTTGCATGTCATTCCACGTAGATAAAAAGAATGCCTATATAAAAA" + \
        "TGGCACATTTTCTTGTAGGTAGTGGAAAGTATCTTTCCAGCAAAGACCATATAATCCGATAAAGGTGAT" + \
        "AACTAAATGTTGAAATCGAGTAGATGCCATATCATCTATACCTTATCTGTTGTTTGGAAAAAAAAGACA" + \
        "AAATCCAAAAAATTTATATGAGATCTCACCTATATAAATAGCGCCCAAATCAGTAGTTAATCCATCACC" + \
        "CATCTAGAGGATCCCCGGGTACCGGTCGCCACCATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTG" + \
        "GTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAG" + \
        "GGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGG" + \
        "CCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAG" + \
        "CACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGAC" + \
        "GGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAG" + \
        "GGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAAC" + \
        "GTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAG" + \
        "GACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTG" + \
        "CCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATG" + \
        "GTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAAAGCGGC" + \
        "CGCGACTCTAGAGTCGACCTGCAGGCATGCAAGCTCGAGTTTCTCCATAATAATGTGTGAGTAGTTCCC" + \
        "AGATAAGGGAATTAGGGTTCCTATAGGGTTTCGCTCATGTGTTGAGCATATAAGAAACCCTTAGTATGT" + \
        "ATTTGTATTTGTAAAATACTTCTATCAATAAAATTTCTAATTCCTAAAACCAAAATCCAGTACTAAAAT" + \
        "CCAGATCCCCCGAATTCCAAGCTCGAACCCCTTCGCTAGAGCTAAGACTTGTGTTTACAATTTTTTATA" + \
        "TTTGTTTTTAAGTTTTGAATATATGTTTTCATGTGTGATTTTACCGAACAAAAATACCGGTTCCCGTCC" + \
        "GATTTCGACTTTAACCCGACCGGATCGTATCGGTTTTCGATTACCGTATTTATCCCGTTCGTTTTCGTT" + \
        "ACCGGTATATCCCGTTTTCGTTTCCGTCCCGCAAGTTAAATATGAAAATGAAAACGGTAGAGGTATTTT" + \
        "ACCGACCGTTCCCGACCGTTTTCATCCCTA"

# >Ds-GFP_TAG21_ReverseComplement  (2165 bp)
dsgfp_revcomp = "TAGGGATGAAAACGGTCGGGAACGGTCGGTAAAATACCTCTACCGTTTTCATTTTCA" + \
                "TATTTAACTTGCGGGACGGAAACGAAAACGGGATATACCGGTAACGAAAACGAACGGGATAAATACGGT" + \
                "AATCGAAAACCGATACGATCCGGTCGGGTTAAAGTCGAAATCGGACGGGAACCGGTATTTTTGTTCGGT" + \
                "AAAATCACACATGAAAACATATATTCAAAACTTAAAAACAAATATAAAAAATTGTAAACACAAGTCTTA" + \
                "GCTCTAGCGAAGGGGTTCGAGCTTGGAATTCGGGGGATCTGGATTTTAGTACTGGATTTTGGTTTTAGG" + \
                "AATTAGAAATTTTATTGATAGAAGTATTTTACAAATACAAATACATACTAAGGGTTTCTTATATGCTCA" + \
                "ACACATGAGCGAAACCCTATAGGAACCCTAATTCCCTTATCTGGGAACTACTCACACATTATTATGGAG" + \
                "AAACTCGAGCTTGCATGCCTGCAGGTCGACTCTAGAGTCGCGGCCGCTTTACTTGTACAGCTCGTCCAT" + \
                "GCCGAGAGTGATCCCGGCGGCGGTCACGAACTCCAGCAGGACCATGTGATCGCGCTTCTCGTTGGGGTC" + \
                "TTTGCTCAGGGCGGACTGGGTGCTCAGGTAGTGGTTGTCGGGCAGCAGCACGGGGCCGTCGCCGATGGG" + \
                "GGTGTTCTGCTGGTAGTGGTCGGCGAGCTGCACGCTGCCGTCCTCGATGTTGTGGCGGATCTTGAAGTT" + \
                "CACCTTGATGCCGTTCTTCTGCTTGTCGGCCATGATATAGACGTTGTGGCTGTTGTAGTTGTACTCCAG" + \
                "CTTGTGCCCCAGGATGTTGCCGTCCTCCTTGAAGTCGATGCCCTTCAGCTCGATGCGGTTCACCAGGGT" + \
                "GTCGCCCTCGAACTTCACCTCGGCGCGGGTCTTGTAGTTGCCGTCGTCCTTGAAGAAGATGGTGCGCTC" + \
                "CTGGACGTAGCCTTCGGGCATGGCGGACTTGAAGAAGTCGTGCTGCTTCATGTGGTCGGGGTAGCGGCT" + \
                "GAAGCACTGCACGCCGTAGGTCAGGGTGGTCACGAGGGTGGGCCAGGGCACGGGCAGCTTGCCGGTGGT" + \
                "GCAGATGAACTTCAGGGTCAGCTTGCCGTAGGTGGCATCGCCCTCGCCCTCGCCGGACACGCTGAACTT" + \
                "GTGGCCGTTTACGTCGCCGTCCAGCTCGACCAGGATGGGCACCACCCCGGTGAACAGCTCCTCGCCCTT" + \
                "GCTCACCATGGTGGCGACCGGTACCCGGGGATCCTCTAGATGGGTGATGGATTAACTACTGATTTGGGC" + \
                "GCTATTTATATAGGTGAGATCTCATATAAATTTTTTGGATTTTGTCTTTTTTTTCCAAACAACAGATAA" + \
                "GGTATAGATGATATGGCATCTACTCGATTTCAACATTTAGTTATCACCTTTATCGGATTATATGGTCTT" + \
                "TGCTGGAAAGATACTTTCCACTACCTACAAGAAAATGTGCCATTTTTATATAGGCATTCTTTTTATCTA" + \
                "CGTGGAATGACATGCAAGATCTCTTCACCTTTACACATGTaagctCCAAGCTCGAACCCCTTCGCTAGA" + \
                "GCCGAGCAGTGGAGTGGAGGGCGGACCGTGGATGGCGGCGTCTGCGTCCGGCGTGCCGCGTCACGGCCG" + \
                "TCACCGCGTGTGGTGCCTGGTGCAGCCCAGCGGCCGGCCGGCTGGGAGACAGGGAGAGTCGGAGAGAGC" + \
                "AGGCGAGAGCGAGACGCGTCGCCGGCGTCGGCGTGCGGCTGGCGGCGTCCGGACTCCGGCGTGGGCGCG" + \
                "TGGCGGCGTGTGAATGTGTGATGCTGTTACTCGTGTGGTGCCTGGCCGCCTGGGAGAGAGGCAGAGCAG" + \
                "CGTTCGCTAGGTATTTCTTACATGGGCTGGGCCTCAGTGGTTATGGATGGGAGTTGGAGCTGGCCATAT" + \
                "TGCAGTCATCCCGAATTAGAAAATACGGTAACGAAACGGGATCATTCCGATTAAAAACGGGATCCCGGT" + \
                "GAAACGGTCGGGAAACTAGCTCTACCGTTTCCGTTTCCGTTTACCGTTTTGTATATCCCGTTTCCGTTC" + \
                "CGTTTTCGTTTTTTACCTCGGGTTCGAAATCGATCGGGATAAAACTAACAAAATCGGTTATACGATAAC" + \
                "GGTCGGTACGGGATTTTCCCATCCTACTTTCATCCCTA"


class Query:
    def __init__(self, line, genome, num_hits):
        if line:
            items = line.split('\t')
        else:
            items = [i == 0 for i in range(12)]
        self.query = items[0]  # R id number
        self.chromosome = items[1]  # blast output
        self.per_identity = items[2]  # blast output
        self.alignment_length = items[3]  # blast output
        self.mismatches = items[4]  # blast output
        self.gap_opens = items[5]  # blast output
        self.q_start = int(items[6])  # blast output
        self.q_end = int(items[7])  # blast output
        self.s_start = int(items[8])  # blast output
        self.s_end = int(items[9])  # blast output
        self.evalue = items[10]  # blast output
        self.bit_score = float(items[11])  # blast output
        self.genome = genome  # blast output
        self.num_hits = num_hits  # blast output
        self.diff = 0  # calculated with parser
        self.strand = 0  # self.strand_direction()                           # calculated with class method
        self.q_start_status = None  # self.set_q_start_status()                 # calculated with class method
        self.bit_score_status = None  # self.set_bit_score_status()             # calculated with class method
        self.primer_name_left = None  # created for primer3 input
        self.primer_name_right = None  # created for primer3 input
        self.side_gfp3utr = None  # blast output
        self.side_3dsgg = None  # blast output
        self.upper_sequence = None  # from filterfasta
        self.lower_sequence = None  # from filter fasta
        self.insertion_sequence = None  # from filter fasta
        self.wildtype_sequence = None  # from filterfasta
        self.primer_sequence_left = "FAIL"  # from Primer3 output
        self.primer_sequence_right = "FAIL"  # from Primer3 output
        self.left_primer_position = -1  # from Primer3 output
        self.right_primer_position = -1  # from Primer3 output
        self.product_size_gfp3utr = -1  # from Primer3 output
        self.product_size_3dsgg = -1  # from Primer3 output
        self.productWT_verify_size = -1  # from Primer3 verification output
        self.upper_coordinates = [-1, -1]  # self.get_upper_coordinates()           # blast output
        self.lower_coordinates = [-1, -1]  # self.get_lower_coordinates()           # blast output
        self.wildtype_coordinates = [-1, -1]  # self.get_wt_coordinates()           # blast output
        self.tm = -1  # TODO: where does this come from?
        self.primer_pair_penalty = -1  # from Primer3 output
        self.primer_left_penalty = -1  # from Primer3 output
        self.primer_right_penalty = -1  # from Primer3 output
        self.best_for_genome = False
        self.best_query = False  # TODO: set this to determine working set
        # self.set_primer_side()

    def __set_values__(self):
        self.strand = self.__strand_direction__()  # calculated with class method
        self.q_start_status = self.__set_q_start_status__()  # calculated with class method
        self.bit_score_status = self.__set_bit_score_status__()
        self.upper_coordinates = self.__get_upper_coordinates__()  # blast output
        self.lower_coordinates = self.__get_lower_coordinates__()  # blast output
        self.wildtype_coordinates = self.__get_wt_coordinates__()  # blast output
        self.__set_primer_side__()

    def __make_best_genome_string__(self):
        return f"{self.bit_score},{self.num_hits},{self.q_start_status},"

    def __print_query__(self):
        print(f"Query: {self.query}\tDatabase: {self.genome}\tBit Score: {self.bit_score}")
        print(f"Upper Sequence: {self.upper_sequence}")
        print(f"Insertion Sequence: {self.insertion_sequence}")

    def __get_upper_coordinates__(self):
        if self.strand == 1:
            return [self.s_start, self.s_start + 2000]
        elif self.strand == -1:
            return [self.s_start - 7, self.s_start + 2000]
        else:
            return [-1, -1]

    def __get_lower_coordinates__(self):
        if self.strand == 1:
            return [self.s_start - 2000, self.s_start + 7]
        elif self.strand == -1:
            return [self.s_start - 2000, self.s_start]
        else:
            return [-1, -1]

    def __get_wt_coordinates__(self):
        return [self.s_start - 2000, self.s_start + 2000]

    def __set_bit_score_status__(self):
        if self.bit_score >= 80:
            return True
        return False

    def __strand_direction__(self):
        # Determine if sequence is on the plus or minus strand
        if int(self.s_start < self.s_end):
            return 1
        else:
            return -1

    def __set_q_start_status__(self):
        # Deterime if sequence match has a q. start of 1
        if self.q_start == 1:
            return True
        else:
            return False

    def __build_insertion_sequence__(self):
        if self.strand == 1:
            self.insertion_sequence = self.lower_sequence + dsgfp + self.upper_sequence
        elif self.strand == -1:
            self.insertion_sequence = self.lower_sequence + dsgfp_revcomp + self.upper_sequence
        else:
            self.insertion_sequence = "FAILED TO MAKE INSERTION SEQUENCE"

    def __set_primer_side__(self):
        if self.strand == 1:
            self.side_gfp3utr = "right"
            self.side_3dsgg = "left"
        elif self.strand == -1:
            self.side_gfp3utr = "left"
            self.side_3dsgg = "right"

    def __toJSON__(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    def __iter__(self):
        # first start by grabbing the Class items
        iters = dict((x, y) for x, y in Query.__dict__.items() if x[:2] != '__')

        # then update the class items with the instance items
        iters.update(self.__dict__)

        # now 'yield' through the items
        for x, y in iters.items():
            yield x, y

    # def __print_Query__(self):
    #    q = self.__iter__()
    #    print(q)

    def __make_filterfasta_input(self):
        return f"{self.query} {self.chromosome} " + \
               f"{self.wildtype_coordinates[0]} {self.wildtype_coordinates[1]} " + \
               f"{self.upper_coordinates[0]} {self.upper_coordinates[1]} " + \
               f"{self.lower_coordinates[0]} {self.lower_coordinates[1]} " + \
               f"{self.genome[:-2]}"

    def __Query_from_JSON__(self, jsonObject):
        self.query = jsonObject["query"]
        self.chromosome = jsonObject["chromosome"]
        self.per_identity = float(jsonObject["per_identity"])
        self.alignment_length = int(jsonObject["alignment_length"])
        self.mismatches = int(jsonObject["mismatches"])
        self.gap_opens = int(jsonObject["gap_opens"])
        self.q_start = int(jsonObject["q_start"])
        self.q_end = int(jsonObject["q_end"])
        self.s_start = int(jsonObject["s_start"])
        self.s_end = int(jsonObject["s_end"])
        self.evalue = jsonObject["evalue"]
        self.bit_score = float(jsonObject["bit_score"])
        self.genome = jsonObject["genome"]
        self.num_hits = int(jsonObject["num_hits"])
        self.diff = float(jsonObject["diff"])
        self.strand = int(jsonObject["strand"])
        # the following few items may have issues with type casting
        self.q_start_status = jsonObject["q_start_status"]
        self.bit_score_status = jsonObject["bit_score_status"]
        self.primer_name_left = jsonObject["primer_name_left"]
        self.primer_name_right = jsonObject["primer_name_right"]
        self.side_gfp3utr = jsonObject["side_gfp3utr"]
        self.side_3dsgg = jsonObject["side_3dsgg"]
        self.upper_sequence = jsonObject["upper_sequence"]
        self.lower_sequence = jsonObject["lower_sequence"]
        self.insertion_sequence = jsonObject["insertion_sequence"]
        self.wildtype_sequence = jsonObject["wildtype_sequence"]
        self.primer_sequence_left = jsonObject["primer_sequence_left"]
        self.primer_sequence_right = jsonObject["primer_sequence_right"]
        self.left_primer_position = int(jsonObject["left_primer_position"])
        self.right_primer_position = int(jsonObject["right_primer_position"])
        self.product_size_gfp3utr = int(jsonObject["product_size_gfp3utr"])
        self.product_size_3dsgg = int(jsonObject["product_size_3dsgg"])
        self.productWT_verify_size = int(jsonObject["productWT_verify_size"])
        # may have issues with typecasting
        self.upper_coordinates = jsonObject["upper_coordinates"]
        self.lower_coordinates = jsonObject["lower_coordinates"]
        self.wildtype_coordinates = jsonObject["wildtype_coordinates"]
        self.tm = jsonObject["tm"]
        self.primer_pair_penalty = jsonObject["primer_pair_penalty"]
        self.primer_left_penalty = jsonObject["primer_left_penalty"]
        self.primer_right_penalty = jsonObject["primer_right_penalty"]
        self.best_for_genome = jsonObject["best_for_genome"]
        self.best_query = jsonObject["best_query"]

    def __create_primer3_input__(self):

        # primer_option, primer_sequence, product_size, hairpin, task = "", "", "", "", ""

        input_str = f"SEQUENCE_ID={self.query}\n"
        input_str += f"SEQUENCE_TEMPLATE={self.insertion_sequence}\n"

        # TODO: do this right!!!!!!
        input_str += f"SEQUENCE_PRIMER={left_side}\nSEQUENCE_PRIMER_REVCOMP={right_side}\n"

        # with open(outfile_name, "a+") as outfile:
        # outfile.write("SEQUENCE_ID=" + query + "\n")
        # outfile.write("SEQUENCE_TEMPLATE=" + sequence + "\n")
        # verify will be boolean, T/F value
        ##if verify:
        #  left_side = side[0]["primer_seq"]
        #  right_side = side[1]["primer_seq"]
        #  outfile.write(f"SEQUENCE_PRIMER={left_side}\nSEQUENCE_PRIMER_REVCOMP={right_side}\n")
        #  product_size = "300-1100"
        #  task = "check_primers"

        if not verify:
            task = "generic"
            if side == "left":
                primer_option = "SEQUENCE_PRIMER_REVCOMP="
            if side == "right":
                primer_option = "SEQUENCE_PRIMER="
            if primer == "a":
                primer_sequence = "TGCAAGCTCGAGTTTCTCCA"  # gfp3utr
                product_size = "750-1000"
                hairpin = "24"
            if primer == "b":
                primer_sequence = "TTGGAGCTGGCCATATTGCAG"  # dsgg3
                product_size = "425-775"
                hairpin = "38"
            outfile.write(primer_option + primer_sequence + "\n")

        outfile.write(
            f"PRIMER_MASK_KMERLIST_PATH=resources/part2/kmer_lists/zea_mays\n")
        outfile.write(
            f"PRIMER_TASK={task}\nPRIMER_PICK_LEFT_PRIMER=1\nPRIMER_PICK_INTER"
            "NAL_OLIGO=0\nPRIMER_PICK_RIGHT_PRIMER=1\nPRIMER_MASK_FAILURE_RATE"
            "=0.1\n")
        outfile.write(f"PRIMER_PRODUCT_SIZE_RANGE={product_size}\n")

        outfile.write("PRIMER_OPT_SIZE=22\nPRIMER_MIN_SIZE=20\nPRIMER_MAX_SIZE"
                      "=24\nPRIMER_OPT_TM=62.0\nPRIMER_MIN_TM=59.0\nPRIMER_MAX_TM=65.0\nPRIM"
                      "ER_TM_FORMULA=1\nPRIMER_PAIR_MAX_DIFF_TM=5.0\n")

        if hairpin:
            outfile.write(f"PRIMER_MAX_HAIRPIN_TH={hairpin}\n")

        outfile.write("PRIMER_EXPLAIN_FLAG=1\nPRIMER_MIN_GC=30.0\nPRIMER_SECON"
                      "DARY_STRUCTURE_ALIGNMENT=1\nPRIMER_MAX_END_STABILITY=9.0\nPRIMER_MIN_"
                      "LEFT_THREE_PRIME_DISTANCE=3\nPRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=3"
                      "\nPRIMER_LIBERAL_BASE=1\nPRIMER_FIRST_BASE_INDEX=1\nPRIMER_MAX_TEMPLA"
                      "TE_MISPRIMING=12.00\nPRIMER_MAX_TEMPLATE_MISPRIMING_TH=47.00\nPRIMER_"
                      "PAIR_MAX_TEMPLATE_MISPRIMING=24.00\nPRIMER_PAIR_MAX_TEMPLATE_MISPRIMI"
                      "NG_TH=47.00\n")

        outfile.write("=\n")

    return