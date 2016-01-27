# FunSeq_PCAWG
a modified version of funseq2 for PCAWG project

Pre-built data context used from orginal Funseq2
 * http://funseq2.gersteinlab.org/data/2.1.2

##### FunSeq modifications:
1. Map current terms to sequence ontology terms
  * splice_variant means either splice_donor or splice_receptor
 
2. Change enhancer annotations to include tissue-specificity
  * code is modified to incorporate specific format.
 
3. Remove hierarchy and annotate all
 
4. Do not remove germline
  * MAF=1
