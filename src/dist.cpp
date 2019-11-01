/*
 * main.cpp
 *
 *  Created on: 15 Mar 2010
 *      Author: sadam
 */

#include <iostream>
#include <tclap/CmdLine.h>
#include <string>
#include <sstream>
#include <fstream>
#include "Fasta.h"
#include "ModelFactoryWag.h"
#include "ModelFactoryPlusF.h"
#include "ModelFactoryDarwin.h"
#include "DistanceFactoryAlign.h"
#include "DistanceFactoryAngle.h"
#include "DateCompiled.h"
#include "Alphabet.h"

#include "main.h"

// ugly global variable
cmdlineopts_t cmdlineopts;

int main(int argc, char** argv)
{
	//command line parsing
	try {

		TCLAP::CmdLine cmd(
				"ProGraphMSA, fast multiple sequence alignment", ' ', LogoDate);

		TCLAP::ValueArg<std::string> outputArg("o", "output", "Output file name",
				false, "/dev/null", "filename");
		cmd.add(outputArg);

		TCLAP::SwitchArg darwinArg("w","darwin","use darwin's model of evolution (instead of WAG)");
		cmd.add(darwinArg);

		TCLAP::SwitchArg aafreqsArg("F", "estimate_aafreqs", "estimate equilibrium amino acid frequencies from input data");
		cmd.add(aafreqsArg);

		TCLAP::ValueArg<double> pcountArg("C", "aafreqs_pseudocount", "pseudocount for estimating equilibrium amino acid frequencies", false, 1125.0, "count");
		cmd.add(pcountArg);

		TCLAP::UnlabeledValueArg<std::string> seqArg("sequences", "input sequences", true, "", "fasta file");
		cmd.add(seqArg);

		cmd.parse(argc, argv);

		cmdlineopts.output_file = outputArg.getValue();
		cmdlineopts.sequence_file = seqArg.getValue();
		cmdlineopts.pseudo_count = pcountArg.getValue();
		cmdlineopts.darwin_flag = darwinArg.getValue();
		cmdlineopts.aafreqs_flag = aafreqsArg.getValue();

		cmdlineopts.mldist_flag = false;
		cmdlineopts.mldist_gap_flag = false;

		std::map<std::string,std::string> seqs = FastaLib(cmdlineopts.sequence_file.c_str()).readAll();

		std::ofstream custom_out;
		std::ostream *out;

		if(outputArg.isSet()) {
			custom_out.open(cmdlineopts.output_file.c_str());
			out = &custom_out;
		} else {
			out = &std::cout;
		}

		ModelFactory<AA> *model_factory;
		if(!cmdlineopts.darwin_flag) {
                   model_factory = new ModelFactoryWag();
		} else {
                   model_factory = new ModelFactoryDarwin();
		}

                std::map<std::string,sequence_t<AA> > seqs2;
                for(std::map<std::string,std::string>::const_iterator it = seqs.begin(); it != seqs.end(); ++it) {
                   seqs2[it->first] = sequenceFromString<AA>(it->second);
                }

		if(cmdlineopts.aafreqs_flag) {
			std::vector<sequence_t<AA> > seq_vector;
			for (std::map<std::string, sequence_t<AA> >::iterator it = seqs2.begin(); it != seqs2.end(); ++it) {
				seq_vector.push_back((*it).second);
			}
			ModelFactoryPlusF<AA> *plusF = new ModelFactoryPlusF<AA>(model_factory);
			plusF->estimateFreqs(seq_vector);
			model_factory = plusF;
		}

                std::vector<std::string> seqs_order(seqs.size());

                index_t i=0;
                for(std::map<std::string,sequence_t<AA> >::const_iterator it = seqs2.begin(); it != seqs2.end(); ++it, ++i) {
                   seqs_order[i] = it->first;
                }

                DistanceFactory<AA> *dist_factory_align = new DistanceFactoryAlign<AA>(model_factory);
                DistanceFactory<AA> *dist_factory_angle = new DistanceFactoryAngle<AA>();

                DistanceMatrix dist_p = dist_factory_align->computePwDistances(seqs2,seqs_order);
                cmdlineopts.mldist_flag = true;
                DistanceMatrix dist_angle = dist_factory_align->computePwDistances(seqs2,seqs_order);
                DistanceMatrix dist_ml = dist_factory_angle->computePwDistances(seqs2,seqs_order);

                for(index_t i=0; i<seqs.size(); ++i) {
                   for(index_t j=i+1; j<seqs.size(); ++j) {
                      *out << dist_p.distances(i,j) << "\t" << dist_ml.distances(i,j) << "\t" << dist_angle.distances(i,j) << "\t";
                      *out << dist_p.variances(i,j) << "\t" << dist_ml.variances(i,j) << "\t" << dist_angle.variances(i,j) << std::endl;
                   }
                }

                delete dist_factory_align;
                delete dist_factory_angle;
                delete model_factory;
        }
        catch (TCLAP::ArgException &e)
        {
           std::cerr << "Command line error: " << e.error() << " for arg " << e.argId() << std::endl;
        }

        return 0;
}
