
import matplotlib.pyplot as plt
import spectrum_utils.plot as sup
import spectrum_utils.spectrum as sus





annotation_settings = {
    "fragment_tol_mass": 10,
    "fragment_tol_mode": "ppm",
    "ion_types": "by",
    "max_ion_charge": 2,
    "neutral_losses": {"NH3": -17.026549, "H2O": -18.010565},
    
}

annot_kws={"fontsize": "xx-large"}

def plot_annotateSpectra(identifier, peptide, precursor_mz, precursor_charge, mz, intensity, output_plot_name):
    my_spec = sus.MsmsSpectrum(identifier = identifier, precursor_mz=  precursor_mz , precursor_charge = precursor_charge , mz=mz, intensity = intensity)
    my_spec = my_spec.filter_intensity(min_intensity = 0.05)
    my_spec.annotate_proforma(proforma_str = peptide, **annotation_settings)
    
    fig, ax = plt.subplots(figsize=(7,6))
    sup.spectrum(my_spec, grid=False, ax=ax, annot_kws={'size': 12})
    ax.set_title(identifier, fontdict={"fontsize": "x-large"})
    ax.set_xlabel("m/z", fontsize=12)  # Set x-axis label font size
    ax.set_ylabel("Intensity", fontsize=12)  # Set y-axis label font size
    ax.tick_params(axis="both", which="major", labelsize=12) 
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    #ax.text(size = 15)
    #ax.annotate(fon)
    if output_plot_name:
        plt.savefig(output_plot_name, dpi=300)
        print(f"Plot saved as '{output_plot_name}'")
    plt.show()
    plt.close()
    


