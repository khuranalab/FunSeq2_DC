# FunSeq_PCAWG
a modified version of funseq2 for PCAWG project

Pre-built data context sources
 1. used by FunSeq_PCAWG
 2. used by FunSeq2 publication 

##### FunSeq modifications:
1. Map current terms to sequence ontology terms
  * splice_variant means either splice_donor or splice_receptor
 
2. Change enhancer annotations to include tissue-specificity
  * code is modified to incorporate specific format.
 
3. Remove hierarchy and annotate all
 
4. Do not remove germline
  * MAF=1
