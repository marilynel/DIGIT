import json

class Query:
    def __init__(self, line, genome, num_hits):
        items = line.split('\t')
        self.query = items[0]                                           # R id number
        self.chromosome = items[1]                                      # blast output
        self.per_identity = items[2]                                    # blast output
        self.alignment_length = items[3]                                # blast output
        self.mismatches = items[4]                                      # blast output
        self.gap_opens = items[5]                                       # blast output
        self.q_start = int(items[6])                                    # blast output
        self.q_end = int(items[7])                                      # blast output
        self.s_start = int(items[8])                                    # blast output
        self.s_end = int(items[9])                                      # blast output
        self.evalue = items[10]                                         # blast output
        self.bit_score = float(items[11])                               # blast output
        self.genome = genome                                             # blast output
        self.num_hits = num_hits                                        # blast output
        self.diff = 0                                                   # calculated with parser
        self.strand = 0 #self.strand_direction()                           # calculated with class method
        self.q_start_status = None #self.set_q_start_status()                 # calculated with class method
        self.bit_score_status = None #self.set_bit_score_status()             # calculated with class method
        self.primer_name_left = None                                    # created for primer3 input
        self.primer_name_right = None                                   # created for primer3 input
        self.side_gfp3utr = None                                        # blast output
        self.side_3dsgg = None                                          # blast output
        self.upper_sequence = None                                      # from filterfasta
        self.lower_sequence = None                                      # from filter fasta
        self.insertion_sequence = None                                  # from filter fasta
        self.wildtype_sequence = None                                   # from filterfasta
        self.primer_sequence_left = "FAIL"                              # from Primer3 output
        self.primer_sequence_right = "FAIL"                             # from Primer3 output
        self.left_primer_position = -1                                  # from Primer3 output
        self.right_primer_position = -1                                 # from Primer3 output
        self.product_size_gfp3utr = -1                                  # from Primer3 output
        self.product_size_3dsgg = -1                                    # from Primer3 output
        self.productWT_verify_size = -1                                 # from Primer3 verification output
        self.upper_coordinates = [-1,-1] #self.get_upper_coordinates()           # blast output
        self.lower_coordinates = [-1,-1] #self.get_lower_coordinates()           # blast output
        self.wildtype_coordinates = [-1,-1] #self.get_wt_coordinates()           # blast output
        self.tm = -1                                                    # TODO: where does this come from?
        self.primer_pair_penalty = -1                                   # from Primer3 output
        self.primer_left_penalty = -1                                   # from Primer3 output
        self.primer_right_penalty = -1                                  # from Primer3 output
        self.best_for_genome = False
        #self.set_primer_side()

    def __set_values__(self):
        self.strand = self.__strand_direction__()                           # calculated with class method
        self.q_start_status = self.__set_q_start_status__()                 # calculated with class method
        self.bit_score_status = self.__set_bit_score_status__() 
        self.upper_coordinates = self.__get_upper_coordinates__()           # blast output
        self.lower_coordinates = self.__get_lower_coordinates__()           # blast output
        self.wildtype_coordinates = self.__get_wt_coordinates__()           # blast output
        self.__set_primer_side__()

    def __make_best_genome_string__(self):
        return f"{self.bit_score},{self.num_hits},{self.q_start_status},"

    def __print_query__(self):
        return f"Query: {self.query}\tDatabase: {self.genome}\tBit Score: {self.bit_score}"

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

    def __build_insertion_sequence__(self, dsgfp_seq):
        return self.lower_sequence + dsgfp_seq + self.upper_sequence

    def __set_primer_side__(self):
        if self.strand == 1:
            self.side_gfp3utr = "right"
            self.side_3dsgg = "left"
        elif self.strand == -1:
            self.side_gfp3utr = "left"
            self.side_3dsgg = "right"

    # TODO:
    # write as json object
    def __toJSON__(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    def __iter__(self):
        # first start by grabbing the Class items
        iters = dict((x,y) for x,y in Query.__dict__.items() if x[:2] != '__')

        # then update the class items with the instance items
        iters.update(self.__dict__)

        # now 'yield' through the items
        for x,y in iters.items():
            yield x,y