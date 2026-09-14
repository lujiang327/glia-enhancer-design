# Pure functions extracted from ChromBPNet eaa0fe58 and modisco 0.5.16.0.
# See adjacent licenses and provenance.json. No TensorFlow dependency.
import numpy as np
import itertools


def compute_per_position_ic(ppm, background, pseudocount):
    """Compute information content at each position of ppm.

    Arguments:
        ppm: should have dimensions of length x alphabet. Entries along the
            alphabet axis should sum to 1.
        background: the background base frequencies
        pseudocount: pseudocount to be added to the probabilities of the ppm
            to prevent overflow/underflow.

    Returns:
        total information content at each positon of the ppm.
    """
    assert len(ppm.shape)==2
    assert ppm.shape[1]==len(background),\
            "Make sure the letter axis is the second axis"
    if (not np.allclose(np.sum(ppm, axis=1), 1.0, atol=1.0e-5)):
        print("WARNING: Probabilities don't sum to 1 in all the rows; this can"
              +" be caused by zero-padding. Will renormalize. PPM:\n"
              +str(ppm)
              +"\nProbability sums:\n"
              +str(np.sum(ppm, axis=1)))
        ppm = ppm/np.sum(ppm, axis=1)[:,None]

    alphabet_len = len(background)
    ic = ((np.log((ppm+pseudocount)/(1 + pseudocount*alphabet_len))/np.log(2))
          *ppm - (np.log(background)*background/np.log(2))[None,:])
    return np.sum(ic,axis=1)

def convolve(to_scan, longer_seq):
    # Convolve to_scan matrix against longer_seq matrix
    vals = []
    for i in range(len(longer_seq) - len(to_scan) + 1):
        vals.append(np.sum(to_scan*longer_seq[i:i+len(to_scan)]))
    return vals

def get_ref_pwms(ref_motifs_path):
    """
    Expected format of file: 
    Contains motifs one position per line, 4 columns per base tab-separated. First 
    line of each motif starts with ">" followed by name. "_" is used to name motifs 
    and name should end with "_plus" or "_minus".

    For ATAC motifs the reference motifs were constructed using the get_pwms 
    function on +4/-5 shifted tagalign files and then taking central 20 bases.
    """
    pwms = {'+':{}, '-':{}}
    cur_orient = None
    cur_motif = None
    with open(ref_motifs_path) as f:
        for x in f:
            x = x.strip()
            if x.startswith(">"):
                # format current as numpy array before starting new
                if cur_motif is not None:
                    pwms[cur_orient][cur_motif] = np.array(pwms[cur_orient][cur_motif])

                if x.endswith("_plus"):
                    cur_orient = "+"
                elif x.endswith("_minus"):
                    cur_orient = "-"
                else:
                    raise ValueError("Invalid reference motif file")
                cur_motif = x[1:]
                pwms[cur_orient][cur_motif] = []
            else:
                pwms[cur_orient][cur_motif].append([float(y) for y in x.split('\t')])
    pwms[cur_orient][cur_motif] = np.array(pwms[cur_orient][cur_motif])

    return pwms['+'], pwms['-']

def compute_shift_ATAC(ref_plus_pwms, ref_minus_pwms, plus_pwm, minus_pwm):
    plus_shifts = set()
    minus_shifts = set()

    for x in ref_plus_pwms:
        # 14 is the value when comparing unshifted BAM pwm
        shift = 14 - np.argmax(convolve(ic_scale(ref_plus_pwms[x]), ic_scale(plus_pwm)))
        plus_shifts.add(shift)
    for x in ref_minus_pwms:
        shift = 5 - np.argmax(convolve(ic_scale(ref_minus_pwms[x]), ic_scale(minus_pwm)))
        minus_shifts.add(shift)

    if len(plus_shifts) != 1 or len(minus_shifts) != 1:
        raise ValueError("Input file shifts inconsistent. Please post an issue")
    
    plus_shift = list(plus_shifts)[0]
    minus_shift = list(minus_shifts)[0]

    if (plus_shift,minus_shift) not in [(0,0)]+ list(itertools.product([3,4,5],[-4,-5,-6])):
        raise ValueError("Input shift is non-standard ({:+}/{:+}). Please post an Issue.".format(plus_shift, minus_shift))

    return plus_shift, minus_shift

def ic_scale(pwm):
    pwm = pwm / np.sum(pwm, axis=-1, keepdims=True)
    return pwm * compute_per_position_ic(pwm, [.25]*4, .001)[:, None]

