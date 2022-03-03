#include "uniquekmers.hpp"
#include <jellyfish/mer_dna.hpp>

using namespace std;

//change to NodeKmers or something and remove the check for uniqueness in the constructor? (No advantage compared to checking simply all kmers for count==1 in the jf database?)

UniqueKmers::UniqueKmers(string read, unsigned int k) {
/*	jellyfish::mer_dna read_kmer("");	
	map<jellyfish::mer_dna, unsigned int> kmer_counts;
	for (int i = 0; i < read.size(); i++){
		char base = read[i];
		//compute the first kmer
		if (i < k-1) {
			read_kmer.shift_left(base);		
		}
		else {
			read_kmer.shift_left(base);
			kmer_counts[read_kmer]+=1;		
		}
	//	read_kmer = read.substr(i,k);
	//	jellyfish::mer_dna jf_kmer(read_kmer);
	//	kmer_counts[jf_kmer] += 1;
	//	uniqueKmerOccurrence[jf_kmer]	+= 1;
	}
	kmer_counts[read_kmer]+=1;*/
	int kmer_size = k;
	string allele = read;
	map<jellyfish::mer_dna, size_t> kmer_counts;
	size_t extra_shifts = kmer_size;
	jellyfish::mer_dna::k(kmer_size);
	jellyfish::mer_dna current_kmer("");
	for (size_t i = 0; i < allele.size(); ++i) {
		char current_base = allele[i];
		if (extra_shifts == 0) {
			kmer_counts[current_kmer] += 1;
		}
		if (  ( current_base != 'A') && (current_base != 'C') && (current_base != 'G') && (current_base != 'T') ) {
			extra_shifts = kmer_size + 1;
		}
		current_kmer.shift_left(current_base);
		if (extra_shifts > 0) extra_shifts -= 1;
	}
	kmer_counts[current_kmer] += 1;	
	
	
	//take those kmers which occur once in the string read
	for (auto& el: kmer_counts) {
		//cout << el.first << el.second << endl;
		if (el.second == 1) {
			uniqueKmerOccurrence[el.first]	= el.second;		
		}	
	}
	
}

unsigned int UniqueKmers::getKmerOccurrence(string kmer){
	jellyfish::mer_dna jf_kmer(kmer);
	return(uniqueKmerOccurrence[jf_kmer]);	
	
}

std::map<jellyfish::mer_dna, unsigned int> UniqueKmers::checkUniqueInGraph(JellyfishKmerCounter countReadKmers){
	std::map<jellyfish::mer_dna, unsigned int> unique_in_graph;
	for (auto& el: uniqueKmerOccurrence) {
		jellyfish::mer_dna kmer = el.first;
		if (countReadKmers.getKmerCount(kmer) == 1) {
			unique_in_graph[kmer] = 1;		
		}	
	}
	return(unique_in_graph);
}


/*map<jellyfish::mer_dna, size_t> counts;
	size_t extra_shifts = kmer_size;
	jellyfish::mer_dna::k(kmer_size);
	jellyfish::mer_dna current_kmer("");
	for (size_t i = 0; i < allele.size(); ++i) {
		char current_base = allele[i];
		if (extra_shifts == 0) {
			counts[current_kmer] += 1;
		}
		if (  ( current_base != 'A') && (current_base != 'C') && (current_base != 'G') && (current_base != 'T') ) {
			extra_shifts = kmer_size + 1;
		}
		current_kmer.shift_left(current_base);
		if (extra_shifts > 0) extra_shifts -= 1;
	}
	counts[current_kmer] += 1;

	// determine kmers unique to allele
	for (auto const& entry : counts) {
		if (entry.second == 1) occurences[entry.first].push_back(index);
	}*/
