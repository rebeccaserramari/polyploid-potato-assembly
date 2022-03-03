#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "argumentparser.hpp"
#include <jellyfish/mer_dna.hpp>
#include "jellyfishkmercounter.hpp"
#include "uniquekmers.hpp"

using namespace std;


int main(int argc, char* argv[]) {

	string readfile = "";
	string shortreadfile = "";
	string allelefile = "";
	string kmerfile = "";
	string countfile = "";

	cerr << "Polyploid genome assembly (I): Finding k-mers unique within Altus" << endl;
	cerr << "Author: Rebecca Serra Mari" << endl;

	ArgumentParser argparser;
	argparser.add_command("Polyassembly [options] -g <graph.gfa> -a <alignments.gaf>");
	argparser.add_optional_argument('r', "", "file with graph unitig sequences in fastq or fasta format");
	argparser.add_optional_argument('s', "", "file with short read data, in fastq or fasta format");
	argparser.add_optional_argument('k', "", "kmerfile to be written to during the kmer counting");
	argparser.add_optional_argument('c', "", "outfile for summed up unique kmer counts in short read sample");

		try {
		argparser.parse(argc, argv);
	} catch (const runtime_error& e) {
		argparser.usage();
		cerr << e.what() << endl;
		return 1;
	} catch (const exception& e) {
		return 0;
	}

	readfile = argparser.get_argument('r');
	shortreadfile = argparser.get_argument('s');
	kmerfile = argparser.get_argument('k');
	countfile = argparser.get_argument('c');
	
	unsigned int k = 71;
	unsigned int threads = 24;

	cout << "count kmers in Altus" << endl;
	JellyfishKmerCounter countReadKmers(readfile, k, threads);
	cout << "kmers in Altus counted" << endl;


	cout << "count kmers in Colomba" << endl;
	JellyfishKmerCounter shortReadKmers(shortreadfile, k, threads);
	cout << "kmers in Colomba counted" << endl;

	cout << "query unique Altus kmers not present in Colomba " << endl;
	countReadKmers.queryFromSequence(readfile, shortReadKmers, kmerfile);


	return(0);


}



