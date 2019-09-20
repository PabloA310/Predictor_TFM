#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: Obtain best k-mers

@author: Pablo Aliaga Gaspar
"""

#%%
# Initial stuff


import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer
import csv


def getKmers(sequence, k):
    """Get K-mers from a sequence, this function is used as auxiliar in
    --create_text function."""
    return [sequence[x:x + k].lower() for x in range(
            len(sequence) - k + 1)]

def create_text(seqFile, k):
    """Prepare data for --CountVectorizer() from table of sequences"""
    prot = pd.read_csv(seqFile, sep = "\t")
    prot['words'] = prot.apply(lambda x: getKmers(x['sequence'], k),
        axis = 1)
    prot = prot.drop('sequence', axis = 1)
    prot_texts = list(prot['words'])
    for item in range(len(prot_texts)):
        prot_texts[item] = ' '.join(prot_texts[item])
    y = prot.iloc[:,0].values
    return prot_texts, y

def get_top_n_words(sequences, n=None):
    """List the top n k-mers in one or various sequences."""
    vec = CountVectorizer().fit(sequences)
    bag_of_words = vec.transform(sequences)
    sum_words = bag_of_words.sum(axis=0) 
    words_freq = [(word, sum_words[0, idx]) for word, 
                  idx in vec.vocabulary_.items()]
    words_freq = sorted(words_freq, key = lambda x: x[1], reverse=True)
    return words_freq[:n]

#%%
# K-mers data
# Import and prepare data
prot_texts_alg, y_alg = create_text("algSeq_LTP.txt", k = 5) 
prot_texts_nlg, y_nlg = create_text("nlgSeq_LTP.txt", k = 5) 

# Obtain k-mers frequency
kmers_freq_alg = get_top_n_words(prot_texts_alg, 500)
kmers_freq_nlg = get_top_n_words(prot_texts_nlg, 500)

# Save lists
with open('Kmers_nlg_ltp', 'w') as myfile:
    wr = csv.writer(myfile, quoting = csv.QUOTE_ALL)
    wr.writerow(kmers_freq_nlg)
    
with open('Kmers_alg_ltp', 'w') as myfile2:
    wr = csv.writer(myfile2, quoting = csv.QUOTE_ALL)
    wr.writerow(kmers_freq_alg)

#%%
