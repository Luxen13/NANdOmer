IUPAC_DICT = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "W": ["A", "T"],
    "K": ["G", "T"],
    "M": ["A", "C"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
}

CANONICAL_BASES = ["A", "T", "C", "G"]


ALIGN_MATCH_SCORE = 10
ALIGN_MISMATCH_SCORE = -20
ALIGN_N_MATCH_SCORE = 6
ALIGN_DEGEN_MATCH_SCORE = 5
ALIGN_GAP_OPEN_COST = 20
ALIGN_GAP_EXT_COST = 10
ALIGN_SCORE_THRESHOLD = 50

CALLER_ANCHOR_LEN = 5
CALLER_MAX_ERR_RATIO = 2
CALLER_MIN_GAP_OVERLAP = (1, 1)
CALLER_MIN_PCT_IDENT =0 # relic, is not called anymore
CALLER_MIN_ANCHOR_PCT_IDENT = 0.6 
CALLER_MIN_ALIGN_PCT_IDENT_END_MIN = 0 # used to be 0.7


