#ifndef JELLYFISHKMERCOUNTER_HPP
#define JELLYFISHKMERCOUNTER_HPP

// after the example from Jellyfish/examples/jf_count_dump/jf_count_dump.cc
// after the example from Jellyfish/examples/query_per_sequence/query_per_sequence.cc
#include <iostream>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/mer_overlap_sequence_parser.hpp>
#include <jellyfish/mer_iterator.hpp>

#include <jellyfish/err.hpp>
#include <jellyfish/file_header.hpp>
#include <jellyfish/stream_manager.hpp>
#include <jellyfish/whole_sequence_parser.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/jellyfish.hpp>
#include "sequence_mers.hpp"

typedef jellyfish::cooperative::hash_counter<jellyfish::mer_dna>                  mer_hash_type;
typedef jellyfish::mer_overlap_sequence_parser<jellyfish::stream_manager<char**>> sequence_parser_type;
typedef jellyfish::mer_iterator<sequence_parser_type, jellyfish::mer_dna>         mer_iterator_type;
typedef jellyfish::whole_sequence_parser<jellyfish::stream_manager<char**> > sequence_parser;


class mer_counter : public jellyfish::thread_exec {
  mer_hash_type&                    mer_hash_;
  jellyfish::stream_manager<char**> streams_;
  sequence_parser_type              parser_;
  const bool                        canonical_;

public:
  mer_counter(int nb_threads, mer_hash_type& mer_hash,
              char** file_begin, char** file_end,
              bool canonical)
    : mer_hash_(mer_hash)
    , streams_(file_begin, file_end)
    , parser_(jellyfish::mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 4096, streams_)
    , canonical_(canonical)
  { }
  
  virtual void start(int thid) {
    mer_iterator_type mers(parser_, canonical_);

    for( ; mers; ++mers)
      mer_hash_.add(*mers, 1);
    mer_hash_.done();
	 
  }
};

class mer_initializer : public jellyfish::thread_exec {
  mer_hash_type&                    mer_hash_;
  jellyfish::stream_manager<char**> streams_;
  sequence_parser_type              parser_;
  const bool                        canonical_;
 
public:
  mer_initializer(int nb_threads, mer_hash_type& mer_hash,
              char** file_begin, char** file_end,
              bool canonical)
    : mer_hash_(mer_hash)
    , streams_(file_begin, file_end)
    , parser_(jellyfish::mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 4096, streams_)
    , canonical_(canonical)
  { }
 
  virtual void start(int thid) {
    mer_iterator_type mers(parser_, canonical_);

		for( ; mers; ++mers)
      mer_hash_.set(*mers);
    mer_hash_.done();
  }
};

class mer_updater : public jellyfish::thread_exec {
  mer_hash_type&                    mer_hash_;
  jellyfish::stream_manager<char**> streams_;
  sequence_parser_type              parser_;
  const bool                        canonical_;
 
public:
  mer_updater(int nb_threads, mer_hash_type& mer_hash,
              char** file_begin, char** file_end,
              bool canonical)
    : mer_hash_(mer_hash)
    , streams_(file_begin, file_end)
    , parser_(jellyfish::mer_dna::k(), streams_.nb_streams(), 3 * nb_threads, 4096, streams_)
    , canonical_(canonical)
  { }
 
  virtual void start(int thid) {
    mer_iterator_type mers(parser_, canonical_);

		for( ; mers; ++mers)
      mer_hash_.update_add(*mers,1);
    mer_hash_.done();
  }
};

class JellyfishKmerCounter {
	public:
		JellyfishKmerCounter(std::string readfile,unsigned int k,unsigned int threads);
		JellyfishKmerCounter(std::string readfile,std::string kmers, unsigned int k,unsigned int threads);
		uint64_t getKmerCount(jellyfish::mer_dna kmer);
		uint64_t getKmerCount(std::string kmer);
		std::vector<jellyfish::mer_dna> dump(uint64_t count);
		void queryFromSequence(std::string readfile);
		void queryFromSequence(std::string readfile, JellyfishKmerCounter shortReadKmers, std::string kmerfile);
		template<typename Database>
		void query_from_db(std::string readfile, const Database& db, std::string kmerfile);
		void query_from_sequence_db(std::string readfile, std::string database, std::string kmerfile);
	private:
		mer_hash_type* kmer_hash;
		sequence_parser_type* sequence_parser;
		mer_iterator_type* mer_iterator;
}; 

#endif