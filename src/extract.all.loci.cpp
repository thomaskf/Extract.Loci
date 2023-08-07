#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>

#define FASTA_LEN 60

using namespace std;

void removeFrontSpace(string& str) {
     int i = 0;
     while (i < str.length() && str[i] < '!') {
     	i++;
     }
     if (i == 0)
     	return;
     if (i < str.length()) {
     	str = str.substr(i);
     } else {
     	str = "";
     }
     return;
}

string strUpper(string str) {
	string s = "";
	int i;
	for (i=0; i<str.length(); i++)
		s.append(1, toupper(str[i]));
	return s;
}

int main(int argc, char** argv) {
	if  (argc < 3) {
		cout << "Syntax: " << argv[0] << " [in nex data file] [out folder]" << endl;
		return 1;
	} 

	// for data
	vector<string> names;
	vector<string> seqs;
	
	// for loci
	vector<string> loci;
	vector<string> long_loci_name;
	vector<int> pos_fr;
	vector<int> pos_to;
	vector<char> isCodon;
	string locus;
	string loc;
	string pos_str;
	bool is_codon;
	string p_fr_str, p_to_str;
	string seq, f_name;
	int p_fr, p_to;
	
	bool in_data_block = false;
	bool in_set_block = false;
	bool in_loci = true;
	string aline;
	int i, j, k, l;
	int p1, p2;
	ifstream fin;
	fin.open(argv[1]);
	
	// check whether the file can be opened successfully
	if (!fin.is_open()) {
		cout << "Error opening the file: " << argv[1] << endl;
		exit(1);
	}
	
	//search for the DATA block
	while (getline(fin,aline)) {
		// remove the spaces in front
		removeFrontSpace(aline);
		
		// start the data block
		if (strUpper(aline) == "BEGIN DATA;") {
			in_data_block = true;
			continue;
		}
		
		// load the data block
		if (in_data_block) {
			// skip the line dimensions
			if (aline.length() >= 10 && strUpper(aline.substr(0,10)) == "DIMENSIONS")
				continue;
			// skip the line format
			if (aline.length() >= 6 && strUpper(aline.substr(0,6)) == "FORMAT")
				continue;
			// skip the line matrix
			if (aline.length() >= 6 && strUpper(aline.substr(0,6)) == "MATRIX")
				continue;
			// skip the line ;
			if (aline.length() >= 1 && aline.substr(0,1) == ";")
				continue;
			// reach end of the data block
			if (aline.length() >= 4 && strUpper(aline.substr(0,4)) == "END;") {
				in_data_block = false;
				continue;
			}
			// this should be the data block
			p1 = aline.find_first_of(" \t");
			p2 = aline.find_last_of(" \t");
			if (p1 != string::npos && p2 != string::npos && p1 > 0 && p2 < aline.length()-1) {
				names.push_back(aline.substr(0,p1));
				seqs.push_back(aline.substr(p2+1));
			}
		}
		
		// start the sets block
		if (strUpper(aline) == "BEGIN SETS;") {
			in_set_block = true;
			continue;
		}
		
		// load the sets block
		if (in_set_block) {
			if (aline.length() >= 6 && strUpper(aline.substr(0,6)) == "[LOCI]") {
				in_loci = true;
				continue;
			}
			if (aline.length() >= 9 && strUpper(aline.substr(0,9)) == "[GENOMES]") {
				in_loci = false;
				continue;
			}
			if (aline.length() >= 11 && strUpper(aline.substr(0,11)) == "[OUTGROUPS]") {
				in_loci = false;
				continue;
			}
			// reach end of the sets block
			if (aline.length() >= 4 && strUpper(aline.substr(0,4)) == "END;") {
				in_set_block = false;
				continue;
			}
			// for loci
			if (in_loci && aline.length() >= 7 && strUpper(aline.substr(0,7)) == "CHARSET") {
				p1 = aline.find_first_of(" =", 8);
				p2 = aline.find_last_of(" =");
				if (p1 != string::npos && p2 != string::npos && p1 > 8 && p2 < aline.length()-1) {
					locus = aline.substr(8, p1-8);
					pos_str = aline.substr(p2+1);
					is_codon = false;
					
					// for pos_str, remove the last ";"
					if (pos_str.length() > 1 && pos_str[pos_str.length()-1] == ';')
						pos_str = pos_str.substr(0, pos_str.length()-1);
						
					// 	check whether it is codon
					if (pos_str.length() > 2 && (pos_str.substr(pos_str.length()-2,2) == "\\3") ||  pos_str.substr(pos_str.length()-2,2) == "/3") {
						is_codon = true;
						pos_str = pos_str.substr(0, pos_str.length()-2);
					}
						
					// get the locus name if codon
					loc = locus;
					if (is_codon) {
						p1 = locus.find_first_of("_");
						if (p1 != string::npos && p1 > 0) {
							loc = locus.substr(0, p1);
						}
					}
					
					// get the start and end pos
					p1 = pos_str.find_first_of("-");
					if (p1 != string::npos && p1 > 0 && p1 < pos_str.length()-1) {
						p_fr_str = pos_str.substr(0,p1);
						p_to_str = pos_str.substr(p1+1,pos_str.length()-p1-1);
						p_fr = atoi(p_fr_str.c_str());
						p_to = atoi(p_to_str.c_str());
					} else {
						p_fr = p_to = atoi(pos_str.c_str());
					}
					
					// save to the vector
					loci.push_back(loc);
					long_loci_name.push_back(locus);
					pos_fr.push_back(p_fr);
					pos_to.push_back(p_to);
					isCodon.push_back((int) is_codon);
				}
			}
		}
	}

	// close the file
	fin.close();
	
	// output each loci
	i = 0;
	while (i<loci.size()) {
		loc = loci[i];
		p_fr = pos_fr[i]; // this is one-based (i.e. start from 1 instead of 0)
		p_to = pos_to[i];
		
		// output nex file
		ofstream fout;
		f_name = string(argv[2]) + "/" + loc + ".nex";
		fout.open(f_name.c_str());
		
		// check whether the file can be opened successfully
		if (!fout.is_open()) {
			cout << "Error creating/writing the file: " << f_name << endl;
			cout << "Please check whether the folder: " << argv[2] << " exists." << endl;
			exit(1);
		}
		
		fout << "#nexus" << endl;
		fout << "begin sets;" << endl;
		while (i < loci.size() && loci[i] == loc) {
			fout << "    charset " << long_loci_name[i] << " = " << pos_fr[i] - p_fr + 1 << "-" << pos_to[i] - p_fr + 1;
			if (isCodon[i]) {
				fout << "\\3";
			}
			fout << endl;
			i++;
		}
		fout << "end;" << endl;
		fout.close();
		
		// output fasta file
		ofstream fa_out;
		f_name = string(argv[2]) + "/" + loc + ".fa";
		fa_out.open(f_name.c_str());

		// check whether the file can be opened successfully
		if (!fa_out.is_open()) {
			cout << "Error creating/writing the file: " << f_name << endl;
			cout << "Please check whether the folder: " << argv[2] << " exists." << endl;
			exit(1);
		}

		for (j = 0; j < names.size(); j++) {
			fa_out << ">" << names[j] << endl;
			seq = seqs[j].substr(p_fr - 1, p_to - p_fr + 1); // because p_fr is zero-based
			k = 0;
			while (k < seq.length()) {
				l = FASTA_LEN;
				if (seq.length() - k < l)
					l = seq.length() - k;
				fa_out << seq.substr(k, l) << endl;
				k += l;
			}
		}
		fa_out.close();
	}
	
    return 0;
}

