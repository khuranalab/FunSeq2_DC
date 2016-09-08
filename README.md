# FunSeq2 with data context version 3
a modified version of data context for funseq2 at PCAWG project

Pre-built data context sources
 1. used by [ FunSeq2_DC version 3 ] (http://khuranalab.med.cornell.edu/data_DC3.html)
 2. used by FunSeq2 publication 

##### FunSeq2 modifications:
5. Fix new spliceOverlap annotations from VAT

4. Map current terms to sequence ontology terms
  * splice_variant means either splice_donor or splice_receptor
 
3. Change enhancer annotations to include tissue-specificity
  * code is modified to incorporate specific format.
 
2. Remove hierarchy and annotate all
 
1. Do not remove germline
  * MAF=1

