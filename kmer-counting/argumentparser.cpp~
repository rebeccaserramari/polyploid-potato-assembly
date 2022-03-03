#include "argumentparser.hpp"


using namespace std;

ArgumentParser::ArgumentParser ()
	:command(""),
	subcommand("")
{}

void ArgumentParser::add_command(string command) {
	this->command = command;
}

/*void ArgumentParser::allow_argument(char name, string command) {
//	this->option_to_commands[name].push_back(command);
	this->cmd_to_options[command].push_back(name);
}
*/

void ArgumentParser::add_subcommand(string command, vector<char> mandatory, vector<char> optional){
	this->subcommand = command;
	this->cmd_to_mandatory[command] = mandatory;
	this->cmd_to_optional[command] = optional;
	this->allowed_subcommands.push_back(command);

}

//void ArgumentParser::add_mandatory_argument(char name, string description, vector<string> subcommands) {
void ArgumentParser::add_mandatory_argument(char name, string description) {
	this->arg_description[name] = description;
	this->options += name;
	this->options += ':';
	this->mandatory.push_back(name);
	/*for (auto cmd: subcommands){
		this->cmd_to_options[cmd].push_back(name);	
	}*/
}

void ArgumentParser::add_optional_argument(char name, string default_val, string description) {
	this->arg_description[name] = description;
	this->options += name;
	this->options += ':';
	this->optional.push_back(name);
	this->optional_defaults[name] = default_val;
}


string ArgumentParser::get_subcommand(int argc, char* argv[]) {
	if (argc < 2) {
		ostringstream oss;
		oss << "Error: no subcommand given! Allowed subcommands: " << endl;
		for (auto it = this->allowed_subcommands.begin(); it != this->allowed_subcommands.end(); ++it) {		
			oss << "\t-" << *it << endl;
		}
		throw runtime_error(oss.str());	
	}
	else {
		return(string(argv[1]));
	}
}

void ArgumentParser::parse(int argc, char* argv[]) {

	this->subcommand = std::string(argv[1]);
	if (find(this->allowed_subcommands.begin(), this->allowed_subcommands.end(), this->subcommand) == this->allowed_subcommands.end()) {
		ostringstream oss;
		oss << "Error: the given subcommand " << this->subcommand << " is not allowed" << endl;
		throw runtime_error(oss.str());	
	}

	int c;
	vector<char> existing_options = cmd_to_mandatory[this->subcommand];
   existing_options.insert(existing_options.end(), cmd_to_optional[this->subcommand].begin(), cmd_to_optional[this->subcommand].end());
 	while ((c = getopt(argc, argv, this->options.c_str())) != -1) {
		// check if option exists
		if (find(existing_options.begin(), existing_options.end(), c) != existing_options.end()) {
				this->arguments[c] = string(optarg);			
		} 
		else {
			if (c == 'h') {
				ostringstream oss;
				oss << "help: " << endl;
				print_help();
			}
			else {
				ostringstream oss;
				oss << "Error: option" << (char) optopt << " does not exist. " << endl;
			}			
		}
	}
	// check if mandatory options are all given
	for (auto it = cmd_to_mandatory[this->subcommand].begin(); it != cmd_to_mandatory[this->subcommand].end(); ++it) {
		// check if it was seen
		if (this->arguments.find(*it) == this->arguments.end()) {
			ostringstream oss;
			oss << "Error: option -" << *it << " is mandatory." << endl;
			throw runtime_error(oss.str());
		}
	}
}

string ArgumentParser::get_arg_parameter(char name) {
	if (this->arguments.find(name) == this->arguments.end()) {
		return(optional_defaults[name]);	
	}
	else {
		return(this->arguments[name]);	
	}
}
void ArgumentParser::print_help() {
	cerr << this->command << endl;
	cerr << endl;
	cerr << "options:" << endl;
	cerr << "required arguments: " << endl;
	for (auto it = this->cmd_to_mandatory[this->subcommand].begin(); it != this->cmd_to_mandatory[this->subcommand].end(); it++) {
		auto ats = arg_description.find(*it);
		cerr << "\t-" << ats->first << " <VAL>\t" << ats->second << endl;
	}
	cerr << endl;
	cerr << "optional arguments: " << endl;
	for (auto it = this->cmd_to_optional[this->subcommand].begin(); it != this->cmd_to_optional[this->subcommand].end(); it++) {
		auto ats = arg_description.find(*it);
		cerr << "\t-" << ats->first << " <VAL>\t" << ats->second << endl;
	}
	cerr << endl;

}

