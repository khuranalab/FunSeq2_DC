# FunSeq2 Pre-built data context
1. FunSeq2_DC3 refers to FunSeq2 pre-built data context 3. It contains PCAWG OCT-2016 annotations 
2. FunSeq2_DC2 refers to FunSeq2 pre-built data context 2. It contains PCAWG OCT-2015 annotations
3. FunSeq2-1.0 refers to the original published work by Fu et al, 2014

Pre-built data context sources
FunSeq2_DC3  (http://khuranalab.med.cornell.edu/data_DC3.html)
FunSeq2_DC2  (http://khuranalab.med.cornell.edu/data.html)
FunSeq2-1.0  (http://funseq2.gersteinlab.org/data/2.1.0) 

##### FunSeq2 modifications:
5. Fix new spliceOverlap annotations from VAT

4. Map current terms to sequence ontology terms
  * splice_variant means either splice_donor or splice_receptor
 
3. Change enhancer annotations to include tissue-specificity
  * code is modified to incorporate specific format.
 
2. Remove hierarchy and annotate all
 
1. Do not remove germline
  * MAF=1

