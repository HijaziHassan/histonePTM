#' Quotes from the histone literature
#' @description
#' Print quote(s) cherry-picked from the histone literature.
#'
#'
#' @param quotes_num Number of quotes to print out (\code{integer}). To read all the quotes, use "All".
#'
#' @return A character string containing a quote about histones and histone PTMs.
#' @examples litRicher()
#' @export

litRicher <- function(quotes_num = 1){

    if (is.character(quotes_num)) {
      if (quotes_num != "All") {
        stop("Error: Argument must be either a number or 'All'.")
      }
    } else if (is.numeric(quotes_num)) {
      if (quotes_num %% 1 != 0) {
        warning("Warning: Argument was a double and has been converted to an integer.")
        quotes_num <- as.integer(quotes_num)
      }
    } else if(quotes_num <= 0){
      stop("Error: `quotes_num` must be either a non-negative integer.")
    }else{
      stop("Error: `quotes_num` only accepts a non-negative integer or 'All'.")

    }



  excerpts <- c(
    'The more PTMs the software is allowed to put on the peptide sequence, the great the error rate. The most modified peptides are the most likely to produce errors compared to the less modified peptide under the same global FDR.',

    'The PSMs that bear the most PTMs are the class with the highest false discovery rate.',

    'The more degrees of freedom given a PSM, the greater the proportion of false discoveries at a given score.',

    'A boring explanation is more often correct than a brilliant explanation.',

    'Biological databases contribute perspective to evaluate the reasonableness of PTMs.',

    'PTM-rich samples require much computational resources and long search times.',

    'Histone extracts consist of limited protein composition and only a small protein database suffices for feature annotation.',

    'Propionylation did allow us to more clearly distinguish in-source decay from enzyme specificity.',

    'PTM site(s) identified on peptides cannot be localized onto a protein when the peptide sequence itself is not unique to one protein.',

    'The components of such computational pipelines [database search algorithms] dramatically affect the identification of PTM sites in the analysis.',

    'The increase in the number of search parameters affects the scores of PDMs, which can itself complicate analysis.',

    'Search for all PTMs is still an area of research, and the most appropriate strategies that balance sensitivity and false discovery rates have not yet been fully explored, although they are being tried.',

    'Some PTMs are known to be altered during sample processing itself.',

    'Peptides with PTMs are generally of low abundance, and their detection is adversely affected by the vast excess of non-modified peptides.',

    'The DIA method should conceptually have the advantage to comprehensively recover most of peptide ions for sequencing.',

    'The use of multiple database search algorithms can help increase peptide identifications as well as PTM analysis, as no single search engine is perfect.',

    'A change in abundance of the modified peptide does not necessarily imply a change in the extent of the PTM site (change can come from change of the protein abundance itself, addition of adjacent PTM, etc.).',

    'One should understand that a standard MS-based PTM analysis does not provide occupancy of the PTM site (i.e., stoichiometry). This is because modified and unmodified peptides behave unequally in the proteomic pipeline, making it difficult to compare them directly based on their relative MS responses.',

    'MS-based histone analysis requires distinct sample preparations, acquisition, and data analysis workflows when compared to traditional MS-based approaches.',

    'The vast amount of highly similar peptides leads to several unconventional situations, including more coelution of isobaric peptidoforms during acquisition. Consequently, this results in a higher amount of chimeric MS/MS spectra and ambiguous spectrum annotations leading to problematic identifications.',

    'From a quantitative point of view, the total histone protein abundance can be overall considered constant in a cell; at least, their upper limit is proportional to the amount of DNA that is present in the cell. In fact, hPTMs induce changes in peptide abundance, thereby creating an intrinsic peptide-centric setting for histone analysis. Thus, while proteoforms are resolved using peptide data, peptidoforms need to be resolved by transition data.',

    'Standard FDR control is inapplicable for searching combinatorial PTMs.',

    'As far as false positives go, we believe that no accurate target-decoy strategy exists to assess FDR for combinatorial PTM searches in traditional DDA approaches.',

    'Extracted histone samples are unique when compared to traditional samples, both from the perspective of identification and quantification, and thus require adapted sample preparation, including extraction, purification, and chemical derivatization of lysines.',

    'Often, manual expertise on the sequential elution of different peptidoforms is used to resolve [ambiguous spectrum annotations] in a targeted way.',

    'Unfortunately, the quantification of histone peptidoforms in DDA has to be performed at the MS1 level and therefore is highly susceptible to interference due to many coeluting isobaric peptides.',

    'Searching for many possible modifications using conventional database search tools has proven impractical.',

    'The most abundant PTMs on histones are methylations (me1/me2/me3) and acetylation (ac), which are usually analyzed without any requirement for enrichment. However, histones are modified in traces by almost any modification type discovered on any other protein to date. Because of this abundance of different PTMs, it is common to observe isobaric peptides, i.e., peptides carrying the same number and type of modifications but differently localized.',

    'Due to the high diversity, dynamic range, and combinatorial patterns of [histone] PTMs, their role in cellular biology is still not completely defined.',

    'There is no universally suitable approach for PTM characterization of all core histones and their variants, and the analytical setup must be tailored in accordance with the main experimental objectives.',

    'Although protocols including hydroxylamine treatment or boiling of propionylated peptides for the reversion of nonspecific O-acylation have been established, such approaches are not suitable for histone preparation as even natural PTMs are lost during hydrolysis of ester bonds.',

    'Proteomics data in general is noisy and complex, with peptides of interest often existing as low-intensity peaks in dense spectra. Histone proteomics is no different. Current software may not alert users to poor quality data and will report results for data too noisy to reproducibly quantify. In addition, popular processing methods can introduce errors during normalization and use simple statistics that can lead to many false-positive results.',

    'In high-quality samples, histone peptides are present as sharp peaks whose intensity is well above instrument noise. Low signal intensity leads to missing peaks, while high intensity leads to peak widening and inaccurate quantification. Peaks should also be measured with high mass accuracy and exhibit good peak separation by RPLC.',

    'Particular attention must be taken when considering relative amounts of PTMs within single experiments, as differently modified peptides might present dramatically different ionization efficiencies. Analysis of larger polypeptides or intact proteins, such as the middle-down or the top-down strategy, might overcome the issue of ionization variability.',

    'Histone PTMs affect chromatin structure that influences enzyme recruitment, gene regulation, DNA repair, and chromosome condensation.',

    'Regardless of the method used [DDA or DIA], histone analysis is still challenging becuase of the extremely high number of histone variants and complex combinatorial patterns of their modifications.',

    'Histone post-translational modifications (hPTMs) are epigentic marks that strongly effect numerous processes, including cell cycling and protein interaction ... together with DNA methylation and action of non-coding RNA, hPTMs are key factors in the regulation of processes that directly involve DNA, including gene expression, DNA repair, and replication.',

    'Although various chemical agents (and conditions) for labeling have been tested, none have met all the requirements for appropriate and straightforward quantification of isobaric peptides at the MS1 level',

    'The identification of peptides in coeluting peaks is hindered by intensity-based precursor selection in DDA mode. During peak elution, limited numbers of mass spectra can be obtained, and not all precursors in a coleuting peak will yield good fragmentation spectra that enable identification of the peptide forms.',

    'The ability to detect and quantify the differential distribution of histone PTMs makes MS an important tool not only in chromatin research but also in the discovery of novel pharmacological agents that inhibit histone-modifying enzymes such as histone deacetylases (HDACs).',

    'Very few data mining approaches exist for in-depth [histone] raw data analysis, and most published methods do not demonstrate any analytical figures of merit for reproducibility or their ability to assess the extent of modifications.',

    'targeting essential epigenetic mediators is being hampered because mapping of hPTMs using mass spectrometry (MS) is inconvenienced by complex data analysis and impractical raw data reuse.',

    'In essence, an LC-MS system measures the intensity and physicochemical properties of analytes like m/z, retention time (tR) and fragmentation pattern, turning a physical sample into a digital transcript. Yet, compared to mining sequencing data, extracting biological information from this raw data matrix is very challenging and is hypothesis-driven, especially for histones, given the combinatorial complexity of the histone code.',

    'The current way of [histone] data acquisition and sharing is not tailored to reuse, because (i) the datasets lack power because they are small, often limited by the nanoflow LC stability, (ii) data quality metrics are not always provided, (iii) the selection of hPTMs included in the search is research-driven and (iv) histones require extensive manual data curation, because of the many (near-)isobaric peptidoforms and ambiguous annotations they create.'

  )


  if(quotes_num == 'All'){

    print(excerpts)

  }else{

    quote = sample(excerpts,
           size = quotes_num,
           replace = FALSE)

    (quote)
  }


}


