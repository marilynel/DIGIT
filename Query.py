import json
from Sequences import Sequences


class Query:
    def __init__(self, line, genome, num_hits):
        if line:
            items = line.split('\t')
        else:
            items = [i == 0 for i in range(12)]
        self.query = items[0]  # allele / R number / query
        self.chromosome = items[1]  # chr location of allele instance
        self.per_identity = items[2]  # blast output
        self.alignment_length = items[3]  # blast output
        self.mismatches = items[4]  # blast output
        self.gap_opens = items[5]  # blast output
        self.q_start = int(items[6])  # blast output
        self.q_end = int(items[7])  # blast output
        self.s_start = int(items[8])  # beginning coordinate for query
        self.s_end = int(items[9])  # ending coordinate for query
        self.evalue = items[10]  # blast output
        self.bit_score = float(items[11])  # blast output
        self.genome = genome  # genome database where allele instance was found
        self.num_hits = num_hits  # total number of instances for query in genome
        self.diff = 0  # TODO
        self.strand = 0  # plus or minus strand, calculated with s_start and s_end
        self.q_start_status = False  # true if s_start = 1, else false
        self.bit_score_status = False  # true if bit_score >= 80, else false
        self.primer_name_left = None  # left primer name, determined by strand direction
        self.primer_name_right = None  # right primer name, determined by strand direction
        self.side_match_gfp3utr = None  # left = matches left primer, is right; right = matches right primer, is left
        self.side_match_3dsgg = None  # left = matches left primer, is right; right = matches right primer, is left
        self.upper_sequence = None  # from filterfasta
        self.lower_sequence = None  # from filterfasta
        self.insertion_sequence = None  # calculated with filterfasta results
        self.wildtype_sequence = None  # from filterfasta
        self.primer_sequence_left = "FAIL"  # from primer3 output
        self.primer_sequence_right = "FAIL"  # from primer3 output
        self.left_primer_position = -1  # from primer3 output
        self.right_primer_position = -1  # from primer3 output
        self.product_size_left = -1  # from primer3 output
        self.product_size_right = -1  # from primer3 output
        self.productWT_verify_size = -1  # TODO
        self.upper_coordinates = [-1, -1]  # calculated with strand and s_start
        self.lower_coordinates = [-1, -1]  # calculated with strand and s_start
        self.wildtype_coordinates = [-1, -1]  # calculated with strand and s_start
        self.tm = -1  # TODO: where does this come from?
        self.primer_pair_penalty = -1  # TODO
        self.primer_left_penalty = -1  # from primer3 output
        self.primer_right_penalty = -1  # from primer3 output
        self.best_for_genome = False  # true if this allele has best bit score across all instances in genome
        self.best_query = False  # true if this instance in this genome is the best for all versions of the allele

    def __set_values__(self):
        self.__strand_direction__()
        self.__set_q_start_status__()
        self.__set_bit_score_status__()
        self.__get_upper_coordinates__()
        self.__get_lower_coordinates__()
        self.__get_wt_coordinates__()

    def __make_best_genome_string__(self):
        # For printing to best_queries_by_genome.csv
        return f"{self.bit_score},{self.num_hits},{self.q_start_status},"

    def __print_query__(self):
        # For testing
        print(f"Query: {self.query}\tDatabase: {self.genome}\tBit Score: {self.bit_score}")

    # TODO: write out explanation of coordinates
    def __get_upper_coordinates__(self):
        if self.strand == 1:
            self.upper_coordinates = [self.s_start, self.s_start + 2000]
        elif self.strand == -1:
            self.upper_coordinates = [self.s_start - 7, self.s_start + 2000]
        else:
            self.upper_coordinates = [-1, -1]

    def __get_lower_coordinates__(self):
        if self.strand == 1:
            self.lower_coordinates = [self.s_start - 2000, self.s_start + 7]
        elif self.strand == -1:
            self.lower_coordinates = [self.s_start - 2000, self.s_start]
        else:
            self.lower_coordinates = [-1, -1]

    def __get_wt_coordinates__(self):
        self.wildtype_coordinates = [self.s_start - 2000, self.s_start + 2000]

    def __set_bit_score_status__(self):
        if self.bit_score >= 80:
            self.bit_score_status = True

    def __strand_direction__(self):
        # Determine if sequence is on the plus or minus strand
        if int(self.s_start < self.s_end):  # plus strand
            self.side_match_gfp3utr = "right"  # gfp3utr is a left primer, match with right
            self.side_match_3dsgg = "left"  # 3dsgg is a right primer, match with left
            self.primer_name_left = self.query + "b"  # left primer matches iwth 3dsgg
            self.primer_name_right = self.query + "a"  # right primer matches with gfp3utr
            self.strand = 1
        else:  # minus strand
            self.side_match_gfp3utr = "left"  # gfp3utr is a right primer, match with left
            self.side_match_3dsgg = "right"  # 3dsgg is a left primer, match with right
            self.primer_name_left = self.query + "a"  # left primer matches with gfp3utr
            self.primer_name_right = self.query + "b"  # right primer matches witj 3dsgg
            self.strand = -1

    def __set_q_start_status__(self):
        # Deterime if sequence match has a q. start of 1
        if self.q_start == 1:
            self.q_start_status = True

    def __build_insertion_sequence__(self):
        seqObj = Sequences()
        if self.strand == 1:
            self.insertion_sequence = self.lower_sequence + seqObj.dsgfp + self.upper_sequence
        elif self.strand == -1:
            self.insertion_sequence = self.lower_sequence + seqObj.dsgfp_revcomp + self.upper_sequence
        else:
            self.insertion_sequence = "FAILED TO MAKE INSERTION SEQUENCE"

    def __make_filterfasta_input__(self):
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
        self.side_match_gfp3utr = jsonObject["side_match_gfp3utr"]
        self.side_match_3dsgg = jsonObject["side_match_3dsgg"]
        self.upper_sequence = jsonObject["upper_sequence"]
        self.lower_sequence = jsonObject["lower_sequence"]
        self.insertion_sequence = jsonObject["insertion_sequence"]
        self.wildtype_sequence = jsonObject["wildtype_sequence"]
        self.primer_sequence_left = jsonObject["primer_sequence_left"]
        self.primer_sequence_right = jsonObject["primer_sequence_right"]
        self.left_primer_position = int(jsonObject["left_primer_position"])
        self.right_primer_position = int(jsonObject["right_primer_position"])
        self.product_size_left = int(jsonObject["product_size_left"])
        self.product_size_right = int(jsonObject["product_size_right"])
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

    # TODO: do I want to make a "primer3 input" object? would help clean up this code
    def __create_primer3_input__(self):
        seqObj = Sequences()
        right_primer_input = f"SEQUENCE_ID={self.primer_name_right}\n" + \
                             f"SEQUENCE_TEMPLATE={self.insertion_sequence}\n"
        left_primer_input = f"SEQUENCE_ID={self.primer_name_left}\n" + \
                            f"SEQUENCE_TEMPLATE={self.insertion_sequence}\n"

        if self.strand == 1:
            # 3dsgg is on right / revcomp       ->  matches with left, b
            # gfp3utr is on left / normal       ->  matches with right, a
            right_primer_input += f"SEQUENCE_PRIMER={seqObj.gfp3utr}\n"
            left_primer_input += f"SEQUENCE_PRIMER_REVCOMP={seqObj.dsgg3}\n"
        elif self.strand == -1:
            # 3dsgg is left / normal        ->  matches with right, b
            # gfp3utr is right / revcomp    ->  matches with left, a
            right_primer_input += f"SEQUENCE_PRIMER={seqObj.dsgg3}\n"
            left_primer_input += f"SEQUENCE_PRIMER_REVCOMP={seqObj.gfp3utr}\n"
        else:
            print("Error in Query method __create_primer3_input__()")
            right_primer_input += f"__error__\n"
            left_primer_input += f"__error__\n"

        input_str = \
            f"PRIMER_MASK_KMERLIST_PATH=genomes/kmer_lists/zea_mays\n" + \
            f"PRIMER_TASK=generic\n" + \
            f"PRIMER_PICK_LEFT_PRIMER=1\n" + \
            f"PRIMER_PICK_INTERNAL_OLIGO=0\n" + \
            f"PRIMER_PICK_RIGHT_PRIMER=1\n" + \
            f"PRIMER_MASK_FAILURE_RATE=0.1\n"

        if self.strand == 1:
            right_primer_input += f"{input_str}PRIMER_PRODUCT_SIZE_RANGE=750-1000\n"  # pair wtih gfp3utr
            left_primer_input += f"{input_str}PRIMER_PRODUCT_SIZE_RANGE=425-775\n"  # pair with dsgg3
        elif self.strand == -1:
            right_primer_input += f"{input_str}PRIMER_PRODUCT_SIZE_RANGE=425-775\n"  # pair with dsgg3
            left_primer_input += f"{input_str}PRIMER_PRODUCT_SIZE_RANGE=750-1000\n"  # pair with gfp3utr
        else:
            print("Error in Query method __create_primer3_input__()")
            right_primer_input += f"__error__\n"
            left_primer_input += f"__error__\n"

        input_str = f"PRIMER_OPT_SIZE=22\n" + \
                    f"PRIMER_MIN_SIZE=20\n" + \
                    f"PRIMER_MAX_SIZE=24\n" + \
                    f"PRIMER_OPT_TM=62.0\n" + \
                    f"PRIMER_MIN_TM=59.0\n" + \
                    f"PRIMER_MAX_TM=65.0\n" + \
                    f"PRIMER_TM_FORMULA=1\n" + \
                    f"PRIMER_PAIR_MAX_DIFF_TM=5.0\n"

        # 24 for gfp3utr, 38 for 3dsgg
        if self.strand == 1:
            right_primer_input += f"{input_str}PRIMER_MAX_HAIRPIN_TH=24\n"  # pair wtih gfp3utr
            left_primer_input += f"{input_str}PRIMER_MAX_HAIRPIN_TH=38\n"  # pair with 3dsgg
        elif self.strand == -1:
            right_primer_input += f"{input_str}PRIMER_MAX_HAIRPIN_TH=38\n"  # pair wtih 3dsgg
            left_primer_input += f"{input_str}PRIMER_MAX_HAIRPIN_TH=24\n"  # pair with gfp3utr
        else:
            print("Error in Query method __create_primer3_input__()")
            right_primer_input += f"__error__\n"
            left_primer_input += f"__error__\n"

        input_str = f"PRIMER_EXPLAIN_FLAG=1\n" + \
                    f"PRIMER_MIN_GC=30.0\n" + \
                    f"PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=1\n" + \
                    f"PRIMER_MAX_END_STABILITY=9.0\n" + \
                    f"PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=3\n" + \
                    f"PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=3\n" + \
                    f"PRIMER_LIBERAL_BASE=1\n" + \
                    f"PRIMER_FIRST_BASE_INDEX=1\n" + \
                    f"PRIMER_MAX_TEMPLATE_MISPRIMING=12.00\n" + \
                    f"PRIMER_MAX_TEMPLATE_MISPRIMING_TH=47.00\n" + \
                    f"PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00\n" + \
                    f"PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=47.00\n" + \
                    f"=\n"

        right_primer_input += input_str
        left_primer_input += input_str

        return right_primer_input + left_primer_input

    # TODO: double check directions!!!
    def __update_Query_primer3__(self, output_dict):
        # print(output_dict.keys())
        # right primer
        if "SEQUENCE_PRIMER" in output_dict and output_dict["PRIMER_PAIR_NUM_RETURNED"] != "0":
            self.primer_right_penalty = output_dict["PRIMER_RIGHT_0_PENALTY"]
            self.primer_sequence_right = output_dict["PRIMER_RIGHT_0_SEQUENCE"]
            self.right_primer_position = output_dict["PRIMER_RIGHT_0"].split(",")[0]
            self.product_size_right = output_dict["PRIMER_PAIR_0_PRODUCT_SIZE"]
        else:
            print(f"{self.query} has no right primer")
        # left primer
        if "SEQUENCE_PRIMER_REVCOMP" in output_dict and output_dict["PRIMER_PAIR_NUM_RETURNED"] != "0":
            # print("left")
            # print(output_dict.keys())
            # print("\n")
            self.primer_left_penalty = output_dict["PRIMER_LEFT_0_PENALTY"]
            self.primer_sequence_left = output_dict["PRIMER_LEFT_0_SEQUENCE"]
            self.left_primer_position = output_dict["PRIMER_LEFT_0"].split(",")[0]
            self.product_size_left = output_dict["PRIMER_PAIR_0_PRODUCT_SIZE"]
        else:
            print(f"{self.query} has no left primer")

    def __toJSON__(self):
        json0bj = json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)
        return json0bj

    def __iter__(self):
        # first start by grabbing the Class items
        iters = dict((x, y) for x, y in Query.__dict__.items() if x[:2] != '__')

        # then update the class items with the instance items
        iters.update(self.__dict__)

        # now 'yield' through the items
        for x, y in iters.items():
            yield x, y
