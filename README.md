# FunSeq2 with data context version 2
a modified version of data context for funseq2 at PCAWG project

Pre-built data context sources
 1. used by [ FunSeq2_DC version 2 ] (http://khuranalab.med.cornell.edu/data.html)
 2. used by FunSeq2 publication 

##### FunSeq2 modifications:
1. Map current terms to sequence ontology terms
  * splice_variant means either splice_donor or splice_receptor
 
2. Change enhancer annotations to include tissue-specificity
  * code is modified to incorporate specific format.
 
3. Remove hierarchy and annotate all
 
4. Do not remove germline
  * MAF=1

