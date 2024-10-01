#' Quotes from the histone literature
#' @description
#' Print quote(s) cherry-picked from the histone literature.
#'
#'
#' @param quotes_num Number of quotes to print out (\code{integer}). To read all the quotes, use "All".)
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
    } else {
      stop("Error: Argument must be either a number or 'All'.")
    }



  excerpts <- c(
    'The more PTMs the software is allowed to put on the peptide sequence, the great the error rate, The most modified peptides are the most likely to produce errors compared to the less modified peptide under the same global FDR',

    'The PSMs that bear the most PTMs are the class with the highest false discovery rate',

    'The more degrees of freedom given a PSM, the greater the proportion of false discoveries at a given score',

    'A boring explanation is more often correct than a brilliant explanation',

    'Biological databases contribute perspective to evaluate the reasonableness of PTMs',

    'PTM-rich samples require much computational resources and long search times',

    'Histone extracts consist of limited protein composition and only a small protein database suffices for feature annotation',

    'Propionylation did allow us to more clearly distinguish in-source decay from enzyme specificity',

    'PTM site(s) identified on peptides cannot be localized onto a protein when the peptide sequence itself is not unique to one protein',

    'The components of such computational pipelines [database search algorithms] dramatically affect the identification of PTM sites in the analysis',

    'The increase in the number of search parameters affects the scores of PDMs, which can itself complicate analysis',

    'Search for all PTMs is still an area of research and the most appropriate strategies that balance sensitivity and false discovery rates have not yet been fully explored although they are being tried',

    'Some PTMs are known to be altered during sample processing itself',

    'Peptides with PTMs are generally of low abundance whose detection is adversely affected by the vast excess of nonmodified peptides',

    'DIA method should conceptually have the advantage to comprehensively recover most of peptide ions for sequencing',

    'Use of multiple database search algorithms can help increase peptide identifications as well as PTM analysis as no single search engine is perfect',

    'A chance in abundance of the modified peptide does not necessarily imply a change in the extent of the PTM site (change can come from change of the protein abundance itself, addition of adjacent PTM, …)',

    'One should understand that a standard MS-based PTM analysis does not provide occupancy of the PTM site (i.e, stoichiometry), This is because modified and unmodified peptides behave unequally in the proteomic pipeline, making it difficult to compare them directly based on their relative MS responses',

    'MS-based histone analysis requires a distinct sample preparations, acquisition, and data analysis workflow when compared to traditional MS-based approaches',

    'The vast amount of highly similar peptides leads to several unconventional situations, including more coelution of isobaric peptidoforms during acquisition, Consequently, this results in a higher amount of chimeric MS/MS spetra and ambiguous spectrum annotations leading to problematic identifications',

    'From a quantitative point of view, the total histone protein abundance can be overall considered constant in a cell; at least their upper limit is proportional to the amount of DNA that is present in the cell, In fact, hPTMs induce changes in peptide abundance, thereby creating an intrinsic peptide-centric setting for histone analysis, Thus, while proteoforms are resolved using peptide data, peptidoforms need to be resolved by transition data',

    'Standard FDR control is inapplicable for searching combinatorial PTMs',

    'As far as false positives go, we believe that no accurate target-decoy strategy exists to assess FDR for combinatorial PTM searches in traditional DDA approaches',

    'Extracted histone samples are unique when compared to traditional samples, both from the perspective of identification and quantification, and thus require adapted sample preparation, including extraction, purification, and chemical derivatization of lysines',

    'Often, manual expertise on the sequential elution of different peptidoforms is used to resolve [ambiguous spectrum annotations] in a targeted way',

    'Unfortunately, the quantification of histone peptidoforms in DDA has to be perfomed at the MS1 level and therefore is highly susceptible to interference due to many coeluting isobaric peptides',

    'Searching for many possible modifications using conventional database search tools has proven impractical',

    'The most abundant PTMs on histones are methylations (me1/me2/me3) and acetylation (ac), which are usually analyzed without any requirement for enrichment, However, histones are modified in traces by almost any modification type discovered on any other protein to date, Because of this abundance of different PTMs, it is common to observe isobaric peptides, i.e., peptides carrying the same number and type of modifications but differently localized',

    'Due to the high diversity, dynamic range, and combinatorial patterns of [histone] PTMs, their role in cellular biology is still not completely defined',

    'There is no universally suitable approach for PTM characterization of all core histones abnd their variants, and the analytical setup must be tailored in accordance with the main experimental objectives',

    'Although protocols including hydroxylamine treatment or boiling of propionylated peptides for the reversion of nonspecific O-acylation have been established, such approaches are not suitable for histone preparation as even natural PTMs are lost during hydrolysis of ester bonds'
  )

  if(quotes_num == 'All'){

    print(excerpts)
  }else{

    sample(excerpts,
           size = quotes_num,
           replace = FALSE)
  }


}


