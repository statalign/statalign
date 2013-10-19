StatAlign v3.0a

http://statalign.github.io/
 
 
### INTRODUCTION ###

StatAlign is an extendable software package for Bayesian analysis of
Protein, DNA and RNA sequences. Multiple alignments, phylogenetic trees
and evolutionary parameters are co-estimated in a Markov Chain Monte
Carlo framework, allowing for reliable measurement of the accuracy of
the results. 

This approach eliminates common artifacts that traditional methods suffer
from, at the cost of increased computational time. These artifacts
include the dependency of the constructed phylogeny on a single (probably
suboptimal) alignment and bias towards the guide tree upon which the
alignment relies.

The models behind the analysis permit the comparison of evolutionarily
distant sequences: the TKF92 insertion-deletion model can be coupled to
an arbitrary substitution model. A broad range of models for nucleotide
and amino acid data is included in the package and the plug-in management
system ensures that new models can be easily added.


### USAGE ###

StatAlign is written in Java and thus requires no installation, but you
must have Java 6 or newer on your system. If you do not have such a Java
framework, please download the most recent one from:
  
  http://www.java.com/getjava/ 

Once you have Java, extract the StatAlign zip archive, and double-click
the StatAlign.jar file. This will launch the graphical interface where
sequences to analyse can be easily loaded from the menus -- see the Help
page for more information.

Alternatively, StatAlign can also be run from a console, which is
recommended for automated, script-driven scenarios. To print the list of
command line options, use the following command:

  java -jar StatAlign.jar -help

An example command line setup:

  java -Xmx512m -jar StatAlign.jar -ot=Fasta -mcmc=10k,100k,1k seqs.fasta

The -Xmx512m option is a standard JVM argument that sets the memory limit
of the program at 512 MiBs. This is our recommended minimum, and increase
it as necessary for large inputs. The -ot options selects the output
alignment format, -mcmc sets MCMC parameters such as the number of
burn-in steps and the number of samples, see the user documentation in
the Help menu and on-line for tips on how to set these values. The input
file must contain the sequences to align in Fasta format.


### LICENSE ###

StatAlign is distributed under the GNU General Public License Version 3.
A copy of the license is found in LICENSE.txt.
