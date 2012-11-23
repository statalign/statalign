
StatAlign v2.0
http://statalign.github.com/statalign/
 
 
### INTRODUCTION ###

StatAlign is an extendable software package for Bayesian analysis of Protein, DNA and RNA
sequences. Multiple alignments, phylogenetic trees and evolutionary parameters are
co-estimated in a Markov Chain Monte Carlo framework, allowing for reliable measurement
of the accuracy of the results. 

This approach eliminates common artifacts that traditional methods suffer from, at the cost
of increased computational time. These artifacts include the dependency of the constructed
phylogeny on a single (probably suboptimal) alignment and bias towards the guide tree upon
which the alignment relies. 

The models behind the analysis permit the comparison of evolutionarily distant sequences:
the TKF92 insertion-deletion model can be coupled to an arbitrary substitution model.
A broad range of models for nucleotide and amino acid data is included in the package and
the plug-in management system ensures that new models can be easily added.


### USAGE ###

StatAlign is written in Java and thus requires no installation, although a Java Virtual
Machine must be available on the host computer. If you do not have one, visit

  www.java.com/getjava/ 

The easiest way to run StatAlign is then by double clicking the statalign.jar file that
is part of the program package. This launches the graphical interface where sequences
to analyse can be easily loaded from the menus -- see the Help page for more information.

Alternatively, StatAlign can also be run from the command line, which is recommended for
automated, script-driven scenarios. To get help on console options, type:

  java -jar statalign.jar -help

An example command line setup:

  java -Xmx512m -jar statalign.jar -ot=Fasta -mcmc=10k,100k,1k seqs.fasta

The -Xmx512m option is a standard JVM argument that sets the memory limit of the program
at 512 MiBs. This is our recommended minimum, and increase it as necessary with large
inputs. The -ot options selects the output alignment format, -mcmc sets MCMC parameters
such as the number of burn-in steps and the number of samples, see the user documentation
in the Help menu and on-line for tips on how to set these values. The input file must
contain the sequences to align in Fasta format.


### LICENSE ###

StatAlign is distributed under the GNU General Public License Version 3. A copy of the
license is found in LICENSE.txt.
