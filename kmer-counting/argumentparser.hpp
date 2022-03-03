#ifndef ARGUMENTPARSER_HPP
#define ARGUMENTPARSER_HPP

#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <unistd.h>

class ArgumentParser {
public:
	ArgumentParser();
	void add_command(std::string command);
	void add_subcommand(std::string command, std::vector<char> mandatory, std::vector<char> optional);
	void add_mandatory_argument(char name, std::string description);
	void add_optional_argument(char name, std::string default_val, std::string description);
	void parse(int argc, char* argv[]);
	std::string get_subcommand(int argc, char* argv[]);
	void allow_argument(char name, std::string command);
	std::string get_arg_parameter(char name);
	void print_help();

private:
	std::string command;
	std::string subcommand;
	std::string options;
	std::map<char, std::string> arg_description;
	std::map<char, std::string> arguments;
	//std::map<char, std::vector<std::string>> option_to_commands;
	std::map<std::string, std::vector<char>> cmd_to_mandatory;
	std::map<std::string, std::vector<char>> cmd_to_optional;
	std::vector<char> mandatory;
	std::vector<std::string> allowed_subcommands;
	std::vector<char> optional;
	std::map<char,std::string> optional_defaults;
};

#endif 
