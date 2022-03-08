#ifndef UNIQUEKMERS_HPP
#define UNIQUEKMERS_HPP

#include <iostream>
#include <map>

#include "jellyfishkmercounter.hpp"
#include <jellyfish/mer_dna.hpp>


class UniqueKmers{
	public:
		UniqueKmers(std::string read, unsigned int k);
		std::map<jellyfish::mer_dna, unsigned int> uniqueKmerOccurrence;
		
		unsigned int getKmerOccurrence(std::string kmer);
		
		// check which of the kmers in UniqueKmerOccurrence are unique in the graph
		std::map<jellyfish::mer_dna, unsigned int> checkUniqueInGraph(JellyfishKmerCounter countReadKmers);
	
};


#endif