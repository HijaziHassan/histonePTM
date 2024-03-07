
import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus





annotation_settings = {
    "fragment_tol_mass": 20,
    "fragment_tol_mode": "ppm",
    "ion_types": "abyIm",
    "max_ion_charge": 2,
    "neutral_losses": {"NH3": -17.026549, "H2O": -18.010565},
}

def plot_annotateSpectra(identifier, peptide, precursor_mz, precursor_charge, mz, intensity):
    my_spec = sus.MsmsSpectrum(identifier = identifier, precursor_mz=  precursor_mz , precursor_charge = precursor_charge , mz=mz, intensity = intensity)
    my_spec = my_spec.filter_intensity(min_intensity = 0.05)
    my_spec.annotate_proforma(proforma_str = peptide, **annotation_settings)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    sup.spectrum(my_spec, grid=False, ax=ax)
    ax.set_title(identifier, fontdict={"fontsize": "xx-large"})
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.show()
  


