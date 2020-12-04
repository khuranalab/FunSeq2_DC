# FunSeq2 Pre-built data context
A modified version of FunSeq2 using new data context
## Background 
FunSeq2 is a tool to prioritize and annotate somatic variants from cancer whole genome sequencing. It uses a pre-built data context, which is generated from different genomic and cancer resources and updated on regular basis. We provide latest FunSeq2 scripts (with bug fixes) and updated data context. 
## Notes
### 2017_04_17
- FunSeq2_DC3 refers to FunSeq2 pre-built data context 3. It contains PCAWG OCT-2016 annotations 
- FunSeq2_DC2 refers to FunSeq2 pre-built data context 2. It contains PCAWG OCT-2015 annotations
- FunSeq2-1.0 refers to the original published work by Fu et al, 2012

Pre-built data context sources
- FunSeq2_DC3  (http://khuranalab.med.cornell.edu/data_DC3.html)
- FunSeq2_DC2  (http://khuranalab.med.cornell.edu/data.html)
- FunSeq2-1.0  (http://funseq2.gersteinlab.org/data/2.1.0)(http://khuranalab.med.cornell.edu/data.html) 

#### FunSeq2 modifications 

- Fix new spliceOverlap annotations from VAT
- Map current terms to sequence ontology terms
  - splice_variant means either splice_donor or splice_receptor
- Change enhancer annotations to include tissue-specificity
  - code is modified to incorporate specific format.
- Remove hierarchy and annotate all
- Do not remove germline
  - MAF=1

## References 
1.	Dhingra, P., Fu, Y., Gerstein, M. & Khurana, E. Using FunSeq2 for Coding and Non-Coding Variant Annotation and Prioritization. [Curr. Protoc. Bioinforma. 57, 15.11.1-15.11.17 (2017).](https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpbi.23)
2.	Khurana, E. et al. Integrative annotation of variants from 1092 humans: application to cancer genomics. [Science. 342(6154):1235587 (2013).](http://science.sciencemag.org/content/342/6154/1235587.full)
3.	Fu, Y. et al. FunSeq2: A framework for prioritizing noncoding regulatory variants in cancer. [Genome Biology 15, (2012).](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0480-5)
