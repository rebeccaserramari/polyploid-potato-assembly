#include "jellyfishkmercounter.hpp"
#include <bits/stdc++.h>

using namespace std;

// after the example from Jellyfish/examples/jf_count_dump/jf_count_dump.cc
// after the example from Jellyfish/examples/query_per_sequence/query_per_sequence.cc

typedef jellyfish::whole_sequence_parser<jellyfish::stream_manager<char**> > sequence_parser;

vector<char*> stringToChars(string filename);

using jellyfish::mer_dna;
using jellyfish::mer_dna_bloom_counter;
namespace err = jellyfish::err;

JellyfishKmerCounter::JellyfishKmerCounter(string readfile,unsigned int k,unsigned int threads) {

	jellyfish::mer_dna::k(k); // Set length of mers (k=25)
	const uint64_t hash_size    = 10000000; // Initial size of hash.
	const uint32_t num_reprobes = 126;
	const uint32_t num_threads  = threads; // Number of concurrent threads
	const uint32_t counter_len  = 7;  // Minimum length of counting field
	const bool     canonical    = true; // Use canonical representation

	// create the hash
	mer_hash_type* mer_hash = new mer_hash_type(hash_size, jellyfish::mer_dna::k()*2, counter_len, num_threads, num_reprobes);
	this->kmer_hash = mer_hash;

	vector<char*> chars_array;
	chars_array = stringToChars(readfile);


	// count the kmers
	mer_counter counter(num_threads, *kmer_hash, &chars_array[0],(&chars_array[0])+1, canonical);
	counter.exec_join(num_threads);

}

JellyfishKmerCounter::JellyfishKmerCounter(string readfile, string kmers, unsigned int k,unsigned int threads) {

	//count only those kmers which are present in kmers

	jellyfish::mer_dna::k(k); // Set length of mers (k=25)
	const uint64_t hash_size    = 10000000; // Initial size of hash.
	const uint32_t num_reprobes = 126;
	const uint32_t num_threads  = threads; // Number of concurrent threads
	const uint32_t counter_len  = 7;  // Minimum length of counting field
	const bool     canonical    = true; // Use canonical representation

	// create the hash
	mer_hash_type* mer_hash = new mer_hash_type(hash_size, jellyfish::mer_dna::k()*2, counter_len, num_threads, num_reprobes);
	this->kmer_hash = mer_hash;

	vector<char*> readfile_array;
	vector<char*> kmers_array;

	readfile_array = stringToChars(readfile);
	kmers_array = stringToChars(kmers);


	// initialize the mer_hash
	mer_initializer initializer(num_threads, *kmer_hash, &kmers_array[0],(&kmers_array[0])+1, canonical);
	initializer.exec_join(num_threads);
	
	//update the hash when reading the readfile
	mer_updater updater(num_threads, *kmer_hash, &readfile_array[0],(&readfile_array[0])+1, canonical);
	updater.exec_join(num_threads);

}

// copying the contents of the string to char array
std::vector<char*> stringToChars(std::string filename) {
	int n = filename.length();
	std::vector<char*> chars_array;
	std::string elem_str = filename;
	char* buffer = new char[elem_str.size() + 1];
	copy(elem_str.begin(), elem_str.end(), buffer);
	buffer[elem_str.size()] = 0;
	chars_array.push_back(buffer);
	chars_array.push_back(nullptr);
	return(chars_array);
}


void JellyfishKmerCounter::queryFromSequence(string readfile) {
	vector<char*> readfile_array;
	readfile_array = stringToChars(readfile);
	char** file_begin = &readfile_array[0];
	char** file_end = (&readfile_array[0])+1;
	sequence_mers mers(true);
	const sequence_mers mers_end(true);
	jellyfish::stream_manager<char**> streams(file_begin, file_end);
	//sequence_parser parser(4, 100, 1, streams);
	jellyfish::whole_sequence_parser<jellyfish::stream_manager<char**> > parser(4, 100, 1, streams); 
	while(true) {
		sequence_parser::job j(parser);
		vector<jellyfish::mer_dna> unique_mers;
		if(j.is_empty()) break;
		for(size_t i = 0; i < j->nb_filled; ++i) {
			mers = j->data[i].seq;
			if(mers != mers_end) {
				if (this->getKmerCount(*mers) == 1) {
					unique_mers.push_back(*mers);			
				}
				//std::cout << this->getKmerCount(*mers);
				++mers;
			}
			for( ; mers != mers_end; ++mers) {
				if (this->getKmerCount(*mers) == 1) {
					unique_mers.push_back(*mers);			
				}
			}
				//std::cout << " " << this->getKmerCount(*mers);
				//std::cout << "\n";
		}
		cout << "size of unique mers: " << unique_mers.size() << endl;
	}
	return;

}

void JellyfishKmerCounter::query_from_sequence_db(string readfile, string database, string kmerfile) {
	std::ifstream in(database, std::ios::in|std::ios::binary);
	const char * database_char = database.c_str();
	jellyfish::file_header header(in);
	if(!in.good())
		err::die(err::msg() << "Failed to parse header of file '" << database << "'");
	mer_dna::k(header.key_len() / 2);
	if(header.format() == "bloomcounter") {
		jellyfish::hash_pair<mer_dna> fns(header.matrix(1), header.matrix(2));
		mer_dna_bloom_counter filter(header.size(), header.nb_hashes(), in, fns);
		if(!in.good())
			err::die("Bloom filter file is truncated");
		in.close();
		query_from_db(readfile, filter, kmerfile);
	} else if(header.format() == binary_dumper::format) {
		jellyfish::mapped_file binary_map(database_char);
		binary_query bq(binary_map.base() + header.offset(), header.key_len(), header.counter_len(), header.matrix(),header.size() - 1, binary_map.length() - header.offset());
		query_from_db(readfile, bq, kmerfile);
	} else {
		err::die(err::msg() << "Unsupported format '" << header.format() << "'. Must be a bloom counter or binary list.");
	}

}

template<typename Database>
void JellyfishKmerCounter::query_from_db(string readfile, const Database& db, string kmerfile) {
	vector<char*> readfile_array;
	readfile_array = stringToChars(readfile);
	char** file_begin = &readfile_array[0];
	char** file_end = (&readfile_array[0])+1;
	sequence_mers mers(true);
	const sequence_mers mers_end(true);
	jellyfish::stream_manager<char**> streams(file_begin, file_end);
	//sequence_parser parser(4, 100, 1, streams);
	jellyfish::whole_sequence_parser<jellyfish::stream_manager<char**> > parser(24, 10, 20, streams); 
	ofstream uniquefile;
	uniquefile.open(kmerfile);
	while(true) {
		sequence_parser::job j(parser);
		vector<jellyfish::mer_dna> unique_mers;
		vector<jellyfish::mer_dna> specific_mers;
		if(j.is_empty()) break;
		int nodecounter = -1;
		for(size_t i = 0; i < j->nb_filled; ++i) {
			nodecounter += 1;
			mers = j->data[i].seq;
			int kmercount = 0;
			if(mers != mers_end) {
				if (this->getKmerCount(*mers) == 1) {
					unique_mers.push_back(*mers);		
					if (db.check(*mers) == 0) {
						specific_mers.push_back(*mers);
						uniquefile << ">" << j->data[i].header << "kmer";
						uniquefile << kmercount << "\n";
						uniquefile << *mers << "\n";
						kmercount++;				
					}	
				}
				//std::cout << this->getKmerCount(*mers);
				++mers;
			}
			for( ; mers != mers_end; ++mers) {
				if (this->getKmerCount(*mers) == 1) {
					unique_mers.push_back(*mers);	
					if (db.check(*mers) == 0) {
						specific_mers.push_back(*mers);	
						uniquefile << ">" << j->data[i].header << "kmer";
						uniquefile << kmercount << "\n";
						uniquefile << *mers << "\n";
						kmercount++;			
					}		
				}
			}
				//std::cout << " " << this->getKmerCount(*mers);
				//std::cout << "\n";
		}
		cout << "size of unique mers: " << unique_mers.size() << endl;
		cout << "size of specific mers: " << specific_mers.size() << endl;
	}
	uniquefile.close();
}


void JellyfishKmerCounter::queryFromSequence(string readfile, JellyfishKmerCounter shortReadKmers, string kmerfile) {
	vector<char*> readfile_array;
	readfile_array = stringToChars(readfile);
	char** file_begin = &readfile_array[0];
	char** file_end = (&readfile_array[0])+1;
	sequence_mers mers(true);
	const sequence_mers mers_end(true);
	jellyfish::stream_manager<char**> streams(file_begin, file_end);
	//sequence_parser parser(4, 100, 1, streams);
	jellyfish::whole_sequence_parser<jellyfish::stream_manager<char**> > parser(4, 100, 1, streams); 
	ofstream uniquefile;
	uniquefile.open(kmerfile);
	while(true) {
		sequence_parser::job j(parser);
		vector<jellyfish::mer_dna> unique_mers;
		vector<jellyfish::mer_dna> specific_mers;
		if(j.is_empty()) break;
		int nodecounter = -1;
		for(size_t i = 0; i < j->nb_filled; ++i) {
			nodecounter += 1;
			mers = j->data[i].seq;
			int kmercount = 0;
			if(mers != mers_end) {
				if (this->getKmerCount(*mers) == 1) {
					unique_mers.push_back(*mers);		
					if (shortReadKmers.getKmerCount(*mers) == 0) {
						specific_mers.push_back(*mers);
						uniquefile << ">" << j->data[i].header << "kmer";
						uniquefile << kmercount << "\n";
						uniquefile << *mers << "\n";
						kmercount++;				
					}	
				}
				//std::cout << this->getKmerCount(*mers);
				++mers;
			}
			for( ; mers != mers_end; ++mers) {
				if (this->getKmerCount(*mers) == 1) {
					unique_mers.push_back(*mers);	
					if (shortReadKmers.getKmerCount(*mers) == 0) {
						specific_mers.push_back(*mers);	
						uniquefile << ">" << j->data[i].header << "kmer";
						uniquefile << kmercount << "\n";
						uniquefile << *mers << "\n";
						kmercount++;			
					}		
				}
			}
				//std::cout << " " << this->getKmerCount(*mers);
				//std::cout << "\n";
		}
		cout << "size of unique mers: " << unique_mers.size() << endl;
		cout << "size of specific mers: " << specific_mers.size() << endl;
	}
	uniquefile.close();

}

uint64_t JellyfishKmerCounter::getKmerCount(jellyfish::mer_dna kmer){
	kmer.canonicalize();	
	
	const auto jf_ary = this->kmer_hash->ary();

	uint64_t val = 0;
	if (jf_ary->get_val_for_key(kmer, &val)) {
		return(val);
	}
	else {
		return(0);
	}
}



vector<jellyfish::mer_dna> JellyfishKmerCounter::dump(uint64_t count) {
	const auto jf_ary = this->kmer_hash->ary();
	vector<jellyfish::mer_dna> kmers;
 
	const auto end = jf_ary->end();
	for(auto it = jf_ary->begin(); it != end; ++it) {
		auto& key_val = *it;
		if (key_val.second == count) {
			kmers.push_back(key_val.first);		
		}
}

return kmers;
}

uint64_t JellyfishKmerCounter::getKmerCount(string kmer){
	jellyfish::mer_dna jkmer(kmer);	
	
	jkmer.canonicalize();	
	
	const auto jf_ary = this->kmer_hash->ary();

	uint64_t val = 0;
	if (jf_ary->get_val_for_key(jkmer, &val)) {
		return(val);
	}
	else {
		//handle the case that kmer is not found
		return(0);
	}
}