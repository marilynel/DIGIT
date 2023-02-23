class Query:
    def __init__(self, line, db_id, num_hits):
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
        self.genome = db_id                                             # blast output
        self.num_hits = num_hits                                        # blast output
        self.diff = 0                                                   # calculated with parser
        self.strand = self.strand_direction()                           # calculated with class method
        self.q_start_status = self.set_q_start_status()                 # calculated with class method
        self.bit_score_status = self.set_bit_score_status()             # calculated with class method
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
        self.upper_coordinates = self.get_upper_coordinates()           # blast output
        self.lower_coordinates = self.get_lower_coordinates()           # blast output
        self.wildtype_coordinates = self.get_wt_coordinates()           # blast output
        self.tm = -1                                                    # TODO: where does this come from?
        self.primer_pair_penalty = -1                                   # from Primer3 output
        self.primer_left_penalty = -1                                   # from Primer3 output
        self.primer_right_penalty = -1                                  # from Primer3 output
        self.best_for_genome = False
        self.set_primer_side()


    def make_best_genome_string(self):
        return f"{self.bit_score},{self.num_hits},{self.q_start_status},"

    def print_query(self):
        return f"Query: {self.query}\tDatabase: {self.genome}\tBit Score: {self.bit_score}"

    def get_upper_coordinates(self):
        if self.strand == 1:
            return [self.s_start, self.s_start + 2000]
        elif self.strand == -1:
            return [self.s_start - 7, self.s_start + 2000]
        else:
            return [-1, -1]

    def get_lower_coordinates(self):
        if self.strand == 1:
            return [self.s_start - 2000, self.s_start + 7]
        elif self.strand == -1:
            return [self.s_start - 2000, self.s_start]
        else:
            return [-1, -1]

    def get_wt_coordinates(self):
        return [self.s_start - 2000, self.s_start + 2000]

    def set_bit_score_status(self):
        if self.bit_score >= 80:
            return True
        return False

    def strand_direction(self):
        # Determine if sequence is on the plus or minus strand
        if int(self.s_start < self.s_end):
            return 1
        else:
            return -1

    def set_q_start_status(self):
        # Deterime if sequence match has a q. start of 1
        if self.q_start == 1:
            return True
        else:
            return False

    def build_insertion_sequence(self, dsgfp_seq):
        return self.lower_sequence + dsgfp_seq + self.upper_sequence

    def set_primer_side(self):
        if self.strand == 1:
            self.side_gfp3utr = "right"
            self.side_3dsgg = "left"
        elif self.strand == -1:
            self.side_gfp3utr = "left"
            self.side_3dsgg = "right"