import evaluationtools as et


def loadDAFS(path, ref="xanes", corrected=False, absorption=False, 
             pol="mixed"):
    """
        Load convoluted DAFS intensities from FDMNES output file.
        
        Inputs:
            ref : str
                miller index of reflection
            corrected : int
                0 - no absorption correction
                1 - normal absorption correction
                2 - double absorption correction
            absorption : bool
                if True, return the self absorption spectra as third output
            pol : str
                Defines which columns to use with respect to polarization of
                the diffracted beam. Can be:
                    "sigma" - polarized perpendicular to scattering plane
                    "pi"    - polarized parallel to scattering plane
                    "mixed  - sum of sigma and pi
        
    """
    if hasattr(ref, "__iter__"):
        ref = "".join(map(str, ref))
    ref = ref.strip("()")
    data = et.loaddat(path, todict=True, comment="")
    if "xanes" in ref.lower():
        return data["Energy"], data["<xanes>"]
    elif all([(s.isdigit() or s=="-") for s in ref]):
        col = "Ic" if corrected else "I"
        col += "(%s)"%ref
        colss = col+"ss_0"
        colsp = col+"sp_0"
        if not (colss in data and colsp in data):
            raise ValueError("Reflection not found in data.")
        
        I = 0
        if pol in ["sigma", "mixed"]:
            I =+ data[colss]
        if pol in ["pi", "mixed"]:
            I =+ data[colss]
        
        if absorption:
            colA = "A(%s)"%ref
            A = data[colA+"in_0"] + data[colA+"ou_0"] # in case of symmetric reflection
    else:
        raise ValueError("Invalid input for `ref'.")
    if absorption:
        return data["Energy"], I, A
    else:
        return data["Energy"], I