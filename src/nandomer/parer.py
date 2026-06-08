#alignment with parasail happens here
import parasail
import re

from .defaults import ( ALIGN_DEGEN_MATCH_SCORE,
                        ALIGN_MATCH_SCORE,
                        ALIGN_MISMATCH_SCORE,
                        ALIGN_N_MATCH_SCORE,
                        ALIGN_GAP_EXT_COST,
                        ALIGN_GAP_OPEN_COST,
                        ALIGN_SCORE_THRESHOLD,
                        CANONICAL_BASES,
                        IUPAC_DICT)

#TODO: decide if i want to even if only local makes sense
_ALIGN_FN = {
    "local":       parasail.sw_trace_striped_16,
    "semi_global": parasail.sg_trace_striped_16,
    "global":      parasail.nw_trace_striped_16,
}


def align_config(match=ALIGN_MATCH_SCORE, mismatch=ALIGN_MISMATCH_SCORE, n_match=ALIGN_N_MATCH_SCORE,
                 degen_match=ALIGN_DEGEN_MATCH_SCORE, GAP_open=ALIGN_GAP_OPEN_COST,
                 GAP_ext=ALIGN_GAP_EXT_COST, score_threshold=ALIGN_SCORE_THRESHOLD):
    return {
        "match": match,
        "mismatch": mismatch,
        "n_match": n_match,
        "degen_match": degen_match,
        "GAP_open": GAP_open,
        "GAP_ext": GAP_ext,
        "score_threshold": score_threshold
    }   

full_alphabet = "".join(IUPAC_DICT.keys())

def _score_pair(b1, b2, config):
    """Score a single base pair using IUPAC expansion sets."""
    if b1 == b2 and b1 in CANONICAL_BASES:
        return config["match"]
    if "N" in (b1, b2):
        return config["n_match"]
    if set(IUPAC_DICT.get(b1, [])) & set(IUPAC_DICT.get(b2, [])):
        # At least one must be degenerate for degen_match
        if b1 in CANONICAL_BASES and b2 in CANONICAL_BASES:
            return config["mismatch"]
        return config["degen_match"]
    return config["mismatch"]


def build_matrix(config, alphabet):
    matrix = parasail.matrix_create(alphabet,
                                    config["match"],
                                    config["mismatch"])

    keys = list(alphabet)
    for x, b1 in enumerate(keys):
        for y, b2 in enumerate(keys):
            matrix[x, y] = _score_pair(b1, b2, config)

    return matrix




"""
Old alignment logic, produced a lot of trailing I ops in cigar leading to incorrect ref/query begin positions. Keeping for now for reference, but will need to be removed eventually.

#TODO: add modes for alignment (local, global, semi-global)
def align(query, reference, config, matrix, mode="local"):
    align_fn = _ALIGN_FN[mode]
    result = align_fn(query,
                      reference,
                      config["GAP_open"],
                      config["GAP_ext"],
                      matrix)
    
    cigar = result.cigar
    return {
        "score": result.score,
        "cigar_string": cigar.decode.decode() if isinstance(cigar.decode, bytes) else str(cigar.decode),
        "ref_begin": cigar.beg_ref,
        "query_begin": cigar.beg_query,
        #the full alignment is not needed, i need to extract the gap stats later.
        #"query_aligned": result.traceback.query,
        #"ref_aligned": result.traceback.ref,
        #"match_line": result.traceback.comp,  # '|' for match, ' ' for mismatch
    }
"""

def align(query, reference, config, matrix, mode="local"):
    align_fn = _ALIGN_FN[mode]
    result = align_fn(query, reference,
                      config["GAP_open"], config["GAP_ext"], matrix)

    cigar = result.cigar
    cigar_string = cigar.decode.decode() if isinstance(cigar.decode, bytes) else str(cigar.decode)
    query_begin = cigar.beg_query
    ref_begin = cigar.beg_ref

    # Parasail's traceback can extend the CIGAR beyond the local alignment, encoding unaligned query prefix/suffix as I ops.
    # This should strip those ops providing a non slippy right side of the gap or anchor. 
    leading = re.match(r'^(\d+)I', cigar_string)
    if leading:
        ins_len = int(leading.group(1))
        cigar_string = cigar_string[leading.end():]
        query_begin += ins_len

    trailing = re.search(r'(\d+)I$', cigar_string)
    if trailing:
        cigar_string = cigar_string[:trailing.start()]

    return {
        "score": result.score,
        "cigar_string": cigar_string,
        "ref_begin": ref_begin,
        "query_begin": query_begin,
    }