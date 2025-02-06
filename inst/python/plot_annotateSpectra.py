import functools
import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus
import spectrum_utils.fragment_annotation as sufa
import spectrum_utils.utils as suu


# function for custom annotation
#1. If precursor charge is 2, annotate with only fragment type with no (+)
#2. Option to or no include H2O and NH3 neutral losses
#3. Option to or no include diagnostic ions
def annotate_ion_type(annotation, prec_z, ion_types, annot_diag = False, neutraloss=False):

    if annotation.ion_type[0] in ion_types and annotation.ion_type[0] != 'I' :
        if abs(annotation.isotope) == 1:
            iso = "+i" if annotation.isotope > 0 else "-i"
        elif annotation.isotope != 0:
            iso = f"{annotation.isotope:+}i"
        else:
            iso = ""
        if prec_z == 2:
            if neutraloss:
                nl = {"-NH3": "*", "-H2O": "o"}.get(annotation.neutral_loss, "")
                return f"{annotation.ion_type}{iso}{nl}"
            else:
                return f"{annotation.ion_type}{iso}"
        if prec_z > 2:
            if neutraloss:
                nl = {"-NH3": "*", "-H2O": "o"}.get(annotation.neutral_loss, "")
                return f"{annotation.ion_type}{iso}{'+' * annotation.charge}{nl}"
            else:
                return f"{annotation.ion_type}{iso}{'+' * annotation.charge}"
            
    elif annotation.ion_type[0] == 'I':
        if hasattr(annotation, "label") and annot_diag:
            return annotation.label
        else:
            return sup.annotate_ion_type(annotation, "I")

    return ""

#plotting function
#
def plot_annotateSpectra(identifier: str, 
    peptide: str, 
    precursor_mz: float, 
    precursor_charge: int, 
    mz: list[float], 
    intensity: list[float],  
    output_plot_name: str,
    min_int: float,
    tol_mz: float,
    tol_mode: str,
    ion_types: str,
    annot_diag: bool, 
    nloss: bool):


    annotation_settings = {
    "fragment_tol_mass": float(tol_mz),
    "fragment_tol_mode": tol_mode,
    "ion_types": ion_types,
    "max_ion_charge": 2,
    "neutral_losses": {"NH3": -17.026549, "H2O": -18.010565},
    
    }

    
    my_spec = sus.MsmsSpectrum(identifier = identifier, 
                               precursor_mz=  precursor_mz , 
                               precursor_charge = precursor_charge , 
                               mz=mz, 
                               intensity = intensity)
    
    my_spec = my_spec.filter_intensity(min_intensity = min_int)
    my_spec.annotate_proforma(proforma_str = peptide, **annotation_settings)
    
    sup.colors["I"] = "#455a64"
    if annot_diag:   # Annotate each peak with diagnostic m/z values
        sup.colors["I"] = "#03DAC6"
        diag_mz_annotations = {
        98.0970: "me1K",
        112.0757: "foK",
        126.0917: "acK",
        140.1075: "prK",
        154.1226: "buK",
        152.1070: "crK",
        156.1019: "laK",
        170.1176: "hbK"}
        for i, peak_mz in enumerate(my_spec.mz):
            for diag_mz, diag_label in diag_mz_annotations.items():
                # Compute mass difference and check tolerance
                if float(abs(md := suu.mass_diff(diag_mz, peak_mz, tol_mode == "Da"))) < float(tol_mz):
                    my_spec.annotation[i]=[
                        sufa.FragmentAnnotation(
                            ion_type="I",   
                            charge=1,
                            mz_delta=(md, tol_mode),
                        )
                    ]
                    #add new attribute 'label' containing custom label
                    my_spec.annotation[i][-1].label = diag_label  
                    break  


    fig, ax = plt.subplots()

    sup.spectrum(my_spec, 
                 grid=False, 
                 ax=ax, 
                 annot_fmt= functools.partial(annotate_ion_type, 
                                              prec_z = my_spec.precursor_charge, 
                                              ion_types = ion_types,
                                              annot_diag = annot_diag, 
                                              neutraloss= nloss), 
                annot_kws={'size': 12})



# plot aesthetics
    ax.set_title(identifier, fontdict={"fontsize": "xx-large"}, pad = 30)
    ax.set_xlabel("m/z", fontsize=16, weight = 'bold')  # Set x-axis label font size
    ax.set_ylabel("Intensity", fontsize=16, weight = 'bold')  # Set y-axis label font size

    ## Hide the right and top spines
    ax.spines[['right', 'top']].set_visible(False)

    ## adjust thickness
    for axis in ['bottom','left']:
        ax.spines[axis].set_linewidth(2)

    ax.tick_params(axis="both", which="major", labelsize=14, width = 2, length = 6) 

    peaks = plt.gca()
    for line in peaks.lines:
        line.set_linewidth(2) 

    
    if output_plot_name:
        plt.savefig(output_plot_name, dpi=300, transparent=False, bbox_inches="tight")
        print(f"'{output_plot_name}' spectrum is saved sucessfully.")
    
    plt.show()
    plt.close()
   


