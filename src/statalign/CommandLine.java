package statalign;

import java.io.File;
import java.lang.reflect.Modifier;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ml.options.OptionData;
import ml.options.OptionSet;
import ml.options.Options;
import ml.options.Options.Multiplicity;
import ml.options.Options.Separator;
import statalign.base.AutomateParameters;
import statalign.base.MainManager;
import statalign.base.Tree;
import statalign.base.Utils;
import statalign.io.DataType;
import statalign.io.RawSequences;
import statalign.model.ext.ModelExtension;
import statalign.model.subst.RecognitionError;
import statalign.model.subst.SubstitutionModel;
import statalign.model.subst.plugins.Dayhoff;
import statalign.model.subst.plugins.Kimura3;
import statalign.postprocess.Postprocess;
import statalign.postprocess.PostprocessManager;

public class CommandLine {

	// Variables

	/** Does this instance belong to a parallel run? */
	private boolean isParallel;

	/** Specifies whether the instance should output info to the stdout. */
	private boolean verbose;

	private List<String> substModNames = new ArrayList<String>();
	private Map<String, Integer> postprocAbbr = new HashMap<String, Integer>();
	public ArrayList<String> pluginParameters;

	// Functions

	/**
	 * Creates a <tt>CommandLine</tt> object with <tt>verbose</tt> set to false.
	 * 
	 * @param isParallel
	 *            does this instance belong to a parallel run?
	 */
	public CommandLine(boolean isParallel) {
		this.isParallel = isParallel;
		verbose = false;
	}

	/**
	 * Fills run-time parameters using a list of command-line arguments. If
	 * error occurs displays error message or usage information.
	 * 
	 * @param args
	 *            list of command-line arguments
	 * @return 0 on success, 1 if usage info and 2 if error msg has been
	 *         displayed
	 */
	public int fillParams(String[] args, MainManager manager) {
		initArrays(manager);

		/*
		 * ArrayList<String> parsedArgs = new ArrayList<String>(); for(int i = 0
		 * ; i < args.length ; i++) { if(!args[i].startsWith("plugin:")) {
		 * System.out.println(args[i]); parsedArgs.add(args[i]); } }
		 */

		Options opt = new Options(args, Multiplicity.ZERO_OR_ONE, 0,
				Integer.MAX_VALUE);
		opt.addSet("run")
				.addOption("debug", Separator.EQUALS)
				.addOption("verbose", Separator.EQUALS)
				.addOption("reportBurnin", Separator.EQUALS)
				.addOption("subst", Separator.EQUALS)
				.addOption("mcmc", Separator.EQUALS)
				.addOption("seed", Separator.EQUALS)
				.addOption("usealign", Separator.EQUALS)
				.addOption("usetree", Separator.EQUALS)
				.addOption("ot", Separator.EQUALS)
				.addOption("log", Separator.EQUALS)
				.addOption("plugin", Separator.COLON, Multiplicity.ZERO_OR_MORE)
				.addOption("help", Separator.COLON, Multiplicity.ZERO_OR_MORE)
				.addOption("list", Separator.COLON, Multiplicity.ZERO_OR_MORE)
				.addOption("automate", Separator.EQUALS, Multiplicity.ZERO_OR_ONE);

		OptionSet set;
		if ((set = opt.getMatchingSet(false, false)) == null) {
			return usage(null,manager);
		}
		if (set.isSet("debug")) {
			Utils.DEBUG = true;
			System.out.println("Debug mode enabled.");
		}
		if (set.isSet("verbose")) {
			Utils.VERBOSE = true;
			System.out.println("Verbose output enabled.");
		}
		if (set.isSet("reportBurnin")) {
			manager.inputData.doReportDuringBurnin = true;
			System.out.println("Enabled printing of logging information during burnin.");
		}
		if (set.isSet("help")) {
			OptionData help = set.getOption("help");
			ArrayList<String> argsVector = new ArrayList<String>();
			for (int i = 0; i < help.getResultCount(); i++) {
				argsVector.add(help.getResultValue(i)); 
			}
			return usage(argsVector,manager);
		}
		else if (set.isSet("list")) {
			OptionData list = set.getOption("list");
			for (int i = 0; i < list.getResultCount(); i++) {
				if (list.getResultValue(i).equals("subst")) {
					System.out.println("Available substitition models:\n");
					for (String model : substModNames) {
						try {
							String simpleName = Class.forName(model).getSimpleName();
							System.out.print("\t" + simpleName);
							if (simpleName.equals("Kimura3")) {
								System.out.print("\t\t(Default for DNA)");							
							}
							else if (simpleName.equals("Dayhoff")) {
								System.out.print("\t\t(Default for proteins)");
							}
						}
						catch (ClassNotFoundException e) {
							System.out.println("Substitution model "+model+" has no matching class.");
						}
						System.out.println();						
					}
				}
				else if (list.getResultValue(i).equals("plugin") || list.getResultValue(i).equals("plugins")) {
					System.out.println("Available plugins:\n");
					List<ModelExtension> pluginList = Utils.findPlugins(ModelExtension.class);
					for (ModelExtension plugin : pluginList) {
						System.out.println(plugin.getPluginID()+" ("+(plugin.isActive() ? "active" : "inactive")+")\n");
					}
				}
			}
			return 1;
		}

//		FastaReader f = new FastaReader();
		try {
			for (String inputFile : set.getData()) {
//				try {
//					RawSequences seq = f.read(inputFile);
//					int errorNum = f.getErrors();
//					if (errorNum > 0)
//						warning(errorNum + " errors occured reading "
//								+ inputFile + ", please check file format");
//					if (verbose) {
//						System.out.println("INFO: read " + seq.size()
//								+ " sequences from " + inputFile);
//					}
//					manager.inputData.seqs.add(seq);
//				} catch (IOException e) {
//					return error("error reading input file " + inputFile);
//				}
				File file = new File(inputFile);
				if(!file.exists())
					return error("input file does not exist: "+inputFile);
				DataType data = manager.dataMan.read(file);
				if(data == null) {
					return error("input file does not appear to be in a known format: "+inputFile);
				} else if(data.getSeqs() != null) {
					manager.inputData.seqs.add(data.getSeqs());
				} else if(data instanceof Tree) {
					if(manager.inputData.tree != null)
						return error("more than one initial tree has been given!");
					manager.inputData.tree = (Tree) data;
				} 
				if (!(data instanceof RawSequences)) {
					manager.inputData.auxData.add(data);
				}
			}
			
			manager.fullPath = set.getData().get(0);
			manager.inputData.setBaseFile(new File(set.getData().get(0)));

			if (set.isSet("subst")) {
				String modelName = set.getOption("subst").getResultValue(0);
				for (String model : substModNames) {
					Class<?> cl = Class.forName(model);
					if (cl.getSimpleName().equalsIgnoreCase(modelName)) {
						manager.inputData.model = (SubstitutionModel) cl
								.newInstance();
						try {
							manager.inputData.model
									.acceptable(manager.inputData.seqs);
						} catch (RecognitionError e) {
							return error("The substitution model "
									+ modelName
									+ " cannot be used with the given input sequences!");
						}
						break;
					}
				}
				if (manager.inputData.model == null) {
					return error("Unknown substitution model: " + modelName
							+ "\n");
				}
			} else {
				manager.inputData.model = new Kimura3();
				try {
					manager.inputData.model.acceptable(manager.inputData.seqs);
					if (verbose) {
						System.out.println("Automatically selected "
								+ manager.inputData.model.getClass()
										.getSimpleName()
								+ " as substitution model.");
					}
				} catch (RecognitionError e) {
					try {
						manager.inputData.model = new Dayhoff();
						manager.inputData.model
								.acceptable(manager.inputData.seqs);
						if (verbose) {
							System.out.println("Automatically selected "
									+ manager.inputData.model.getClass()
											.getSimpleName()
									+ " as substitution model.");
						}
					} catch (RecognitionError ee) {
						return error("Default substitution model "
								+ manager.inputData.model.getClass()
										.getSimpleName()
								+ " does not accept the given input sequences:\n"
								+ ee.message);
					}
				}
			}

			if (set.isSet("mcmc")) {
				String mcmcPars = set.getOption("mcmc").getResultValue(0);
				String[] pars = mcmcPars.split(",");

				if (pars.length != 3 && pars.length != 4 && pars.length != 5) {
					return error("MCMC parameters not recognized: " + mcmcPars);
				}

				manager.inputData.pars.burnIn = parseValue(pars[0]);
				manager.inputData.pars.cycles = parseValue(pars[1]);
				manager.inputData.pars.sampRate = parseValue(pars[2]);
				if (pars.length == 4) {
					if (!isParallel) {
						manager.inputData.pars.randomisationPeriod = parseValue(pars[3]);
					}
					else {
						manager.inputData.pars.swapRate = parseValue(pars[3]);
						if (manager.inputData.pars.swapRate < 0) {
							return error("MCMC parameter not recognized: "
									+ mcmcPars);
						}
					}
				}
				if (pars.length == 5) {
					if (!isParallel) {
						return error("Unrecognized MCMC parameters for non-parallel version.");
					}
					else {
						manager.inputData.pars.randomisationPeriod = parseValue(pars[3]);
						manager.inputData.pars.swapRate = parseValue(pars[4]);
						if (manager.inputData.pars.swapRate < 0) {
							return error("MCMC parameter not recognized: "
									+ mcmcPars);
						}
					}
				}

				if (manager.inputData.pars.burnIn < 0
						|| manager.inputData.pars.cycles < 0
						|| manager.inputData.pars.sampRate < 0) {
					return error("MCMC parameters not recognized: " + mcmcPars);
				}
			}

			if (set.isSet("seed")) {
				String seedPar = set.getOption("seed").getResultValue(0);
				try {
					manager.inputData.pars.seed = Integer.parseInt(seedPar);
				} catch (NumberFormatException e) {
					return error("error parsing seed parameter: " + seedPar);
				}
			}
			
			if(set.isSet("usealign")) {				
				return error("Specification of initial alignment is not yet supported");
/* 				
				String useAlign = set.getOption("usealign").getResultValue(0);
				if(useAlign.length() != 1 || (manager.inputData.useAlign =
						"NIF".indexOf(Character.toUpperCase(useAlign.charAt(0)))) == -1)
					return error("bad format for parameter -usealign: "+useAlign);
//				if(manager.inputData.useAlign == 2)
//					return error("fixing the alignment is currently not supported");
				if(manager.inputData.useAlign > 0 && !manager.inputData.seqs.isAligned())
					return error("input sequences must be aligned when -usealign="+useAlign+" is used");
				if(manager.inputData.useAlign == 1)
					return error("initial alignment is not yet supported");
*/
			}
			
			if(set.isSet("usetree")) {
				String useTree = set.getOption("usetree").getResultValue(0);
				if(useTree.length() != 1 || (manager.inputData.useTree =
						"NITE".indexOf(Character.toUpperCase(useTree.charAt(0)))) == -1)
					return error("bad format for parameter -usetree: "+useTree);
				if(manager.inputData.useTree > 0 && manager.inputData.tree == null)
					return error("a tree must be given when -usetree="+useTree+" is used");
				manager.modelExtMan.addToFilenameExtension("useTree"+manager.inputData.useTree);
			} else {
				// default depends on whether a tree file was provided
				manager.inputData.useTree = manager.inputData.tree == null ? 0 : 1;
			}
			
			if (set.isSet("ot")) {
				String outType = set.getOption("ot").getResultValue(0);
				int i;
				for (i = 0; i < MainManager.alignmentTypes.length; i++) {
					if (outType.equalsIgnoreCase(MainManager.alignmentTypes[i])) {
						manager.inputData.currentAlignmentType = i;
						break;
					}
				}
				if (i == MainManager.alignmentTypes.length) {
					return error("Unknown output type: " + outType + "\n");
				}
			}

			// retrieve all parameters starting with plugin:
			OptionData plugins = set.getOption("plugin");
			ArrayList<String> argsVector = new ArrayList<String>();
			for (int i = 0; i < plugins.getResultCount(); i++) {
				argsVector.add(plugins.getResultValue(i));
			}
			
			//Postprocess.pluginParameters = new PluginParameters(argsVector);
			pluginParameters = argsVector;
			//parsePluginParameters(argsVector,manager);

			// TODO allow rnaMode to be switched off even for RNA sequences (plugin param)
			// TODO move rnaMode to a "RNA container" plugin
			if(manager.inputData.seqs.isRNA()) {
//				PostprocessManager.rnaMode = true;
//				System.out.println("RNA mode activated.");
			}

			AutomateParameters.setAutomateBurnIn(false);
			AutomateParameters.setAutomateStepRate(false);
			AutomateParameters.setAutomateNumberOfSamples(false);

			OptionData automation = set.getOption("automate");
			if (automation != null) {
				if (automation.getResultCount() > 0) {
					String[] split = automation.getResultValue(0).split(",");
					ArrayList<String> values = new ArrayList<String>();
					for (int i = 0; i < split.length; i++) {
						values.add(split[i].trim().toLowerCase());
					}
					//System.out.println(values);
					if (values.contains("burn")) {
						AutomateParameters.setAutomateBurnIn(true);
					}
					if (values.contains("rate")) {
						AutomateParameters.setAutomateStepRate(true);
					}
					if (values.contains("cycl")) {
						AutomateParameters.setAutomateNumberOfSamples(true);
					}
				}
				/*
				else if (automation.getResultCount() == 0) {
					System.out.println("automating all parameters");
					AutomateParameters.setAutomateBurnIn(true);
					AutomateParameters.setAutomateStepRate(true);
					AutomateParameters.setAutomateNumberOfSamples(true);
				}*/
				
				if (set.isSet("log")) {
					String log = set.getOption("log").getResultValue(0);
					String[] keys = log.split(",");
					List<Postprocess> pps = manager.postProcMan.getPlugins();

					for (Postprocess pp : pps)
						pp.sampling = false;

					for (String key : keys) {
						try {
							pps.get(postprocAbbr.get(key.toUpperCase())).sampling = true;
						} catch (NullPointerException e) {
							return error("Log file entry code list not recognised: "
									+ log);
						}
					}
				}
			}

		} catch (Exception e) {
			e.printStackTrace();
			return error("Unknown error: " + e.getLocalizedMessage());
		}

		return 0;
	}

	private String getUsageString(ArrayList<String> args, MainManager man) {
		StringBuilder sb = new StringBuilder();
		
		if (args == null || args.size() == 0) {

			sb.append("Usage:\n\n");

			if (isParallel) {
				// TODO: this.
				sb.append("    <depends> statalign.jar [options] seqfile1 [seqfile2 ...] [treefile]\n\n\n");
			} else {
				sb.append("    java -Xmx512m -jar statalign.jar [options] seqfiles [treefile]\n\n");
			}
	
			sb.append("Description:\n\n")
			  .append("    StatAlign can be used for Bayesian analysis of protein, DNA and RNA\n")
			  .append("    sequences. Multiple alignments, phylogenetic trees and evolutionary\n")
			  .append("    parameters are co-estimated in a Markov Chain Monte Carlo framework.\n\n")
			
			  .append("Parameters:\n\n")
			  .append("    seqfiles\n")
			  .append("        Files with the input sequences to align, they must be in Fasta format.\n")
			  .append("    treefile\n")
			  .append("        An optional tree file in Newick format. When given, the tree is used\n")
			  .append("          as the initial tree of the MCMC analysis.\n")
			  .append("        Default: the initial tree is built using the NJ algorithm\n\n")
	
			  .append("Options:\n\n")
	
			  .append("    -subst=MODEL\n")
			  .append("        Select from the present substitution models (for full list, see \"-list:subst\").\n")
			  .append("        Default: " + Kimura3.class.getSimpleName()
					+ " (for DNA/RNA data), " + Dayhoff.class.getSimpleName()
					+ " (for protein data)\n\n");
			
			if (isParallel) {
				sb.append("    -mcmc=burn,cycl,samprate[,randomisationPeriod],swaprate\n");
				sb.append("        Sets MCMC parameters: burn-in, cycles after burn-in, sampling rate, swap rate,\n");
				sb.append("        and (optionally) initial randomisation period.\n");
				sb.append("          Abbreviations k and m mean 1e3 and 1e6 factors.\n");
				sb.append("        Default: 20k,50k,100,0,100\n\n");
			} else {
				sb.append("    -mcmc=burn,cycl,rate[,randomisationPeriod]\n");
				sb.append("        Sets MCMC parameters: burn-in, cycles after burn-in, sampling rate,\n");
				sb.append("        and (optionally) initial randomisation period.\n");
				sb.append("          Abbreviations k and m mean 1e3 and 1e6 factors.\n");
				sb.append("        Default: 20k,50k,100,0\n\n");
			}
			
			sb.append("    -automate=burn,cycl,rate\n")
			  .append("        Automate MCMC parameters: burn-in, cycles after burn-in, sampling rate.\n")
			  .append("        Select which parameters to automate by listing one or more of: burn, cycl, rate\n\n")
			
//			  .append("    -usealign=N|I|F\n")
//			  .append("        Specifies how the alignment of the input sequences is used\n")
//			  .append("          N: it is ignored and sequences are re-aligned from scratch\n")
//			  .append("          I: it is used as the initial alignment and sampled during MCMC\n")
//			  .append("          F: it is fixed throughout the MCMC analysis\n")
//			  .append("        Default: N\n\n")
			
			  .append("    -usetree=N|I|T|E\n")
			  .append("        Defines how the specified tree should be used\n")
			  .append("          N: the tree is constructed from scratch using NJ\n")
			  .append("          I: it is used as the initial tree and sampled during MCMC\n")
			  .append("          T: the topology is fixed throughout the MCMC analysis, edges are sampled\n")
			  .append("          E: both topology and edge lengths are fixed during the MCMC\n")
			  .append("        Default: N (without tree file) or I (when tree file is given)\n\n")
			
			  .append("    -plugin:PLUGIN_NAME\n")
			  .append("        Specify additional plugins to be run and the corresponding parameters for each.\n")
			  .append("        Each plugin should be specified seperately, e.g. -plugin:ppfold -plugin:structal\n")
			  .append("        For additional usage information for a specific plugin, use -help:PLUGIN_NAME.\n\n")
			
			  .append("    -help:PLUGIN_NAME\n")
			  .append("        Prints usage information for the specified plugin.\n\n")
	
			  .append("    -verbose=true\n")
			  .append("        Enables printing of additional MCMC information, such as locations\n")
			  .append("        of topology switches, and acceptance rates at regular intervals.\n")
	
			  .append("    -list:subst\n")
			  .append("        Prints the list of available substition models.\n\n")
	
			  .append("    -list:plugin\n")
			  .append("        Prints the list of available model extension plugins.\n\n")
			
			  .append("    -seed=value\n")
			  .append("        Sets the random seed (same value will reproduce same results for\n")
			  .append("          identical input and settings)\n")
			  .append("        Default: 1\n\n")
	
			  .append("    -reportBurnin=true\n")
			  .append("        Enables the printing of logging information during the burnin.\n")
			  .append("        Default: false\n\n")
	
			  
			  .append("    -ot=OUTTYPE\n")
			  .append("        Sets output alignment type.\n")
			  .append("          (One of: "
					+ Utils.joinStrings(MainManager.alignmentTypes, ", ") + ")\n")
			  .append("        Default: " + MainManager.alignmentTypes[0] + "\n\n")
	
			  .append("    -log=["
					+ Utils.joinStrings(postprocAbbr.keySet().toArray(), "][,")
					+ "]\n")
			  .append("        Lets you customise what is written into the log file (one entry\n")
			  .append("        for each sample).\n")
			  .append(buildPpListStr(man, "          "))
			  .append("        Default: " + buildDefPpList(man) + "\n\n");
		}
		else {
			List<ModelExtension> pluginList = Utils.findPlugins(ModelExtension.class);
			// TODO Maybe also want to display info for postprocessing plugins?
			for (int i=0; i<args.size(); i++) {
				System.out.println("help:"+args.get(i)+"\n");
				boolean pluginFound = false;
				for (ModelExtension plugin : pluginList) {
					if (args.get(i).equals(plugin.getPluginID())) {
						sb.append(plugin.getUsageInfo());
						pluginFound = true;
						continue;
					}
				}
				if (!pluginFound) {
					throw new IllegalArgumentException("No information available for '"+args.get(i)+"'\n");
				}
			}
		}
		return sb.toString();
	}

	private static int parseValue(String string) {
		if (string.isEmpty())
			return -1;
		int factor = 1;
		switch (Character.toUpperCase(string.charAt(string.length() - 1))) {
		case 'K':
			factor = 1000;
			break;
		case 'M':
			factor = 1000000;
			break;
		}
		if (factor > 1)
			string = string.substring(0, string.length() - 1);
		int result = -1;
		try {
			result = factor * Integer.parseInt(string);
		} catch (NumberFormatException e) {
		}
		return result;
	}

	private int usage(ArrayList<String> args, MainManager man) {
		try {
			System.out.println(getUsageString(args,man));
		}
		catch (IllegalArgumentException e) {
			System.out.println(e.getMessage());
		}
		return 1;
	}

	private static int error(String msg) {
		System.out.println("statalign: " + msg);
		return 2;
	}

	private static void warning(String msg) {
		System.out.println("warning: " + msg);
	}

	private void initArrays(MainManager man) {
		findSubstMods();
		fillPostprocAbbr(man.postProcMan);
	}

	private void fillPostprocAbbr(PostprocessManager man) {
		List<Postprocess> plugins = man.getPlugins();

		final String[] keys = new String[plugins.size()];
		Integer[] sorted = new Integer[plugins.size()];

		for (int i = 0; i < plugins.size(); i++) {
			keys[i] = plugins.get(i).getTabName().toUpperCase().replace(' ', '_');
			sorted[i] = i;
		}

		Arrays.sort(sorted, new Comparator<Integer>() {
			@Override
			public int compare(Integer o1, Integer o2) {
				return keys[o1].compareTo(keys[o2]);
			}
		});

		int prevOverlap = 0;
		for (int i = 0; i < plugins.size(); i++) {
			int nextOverlap = 0;
			if (i < plugins.size() - 1)
				nextOverlap = checkOverlap(keys[i], keys[i + 1]);
			String key = keys[i].substring(0,
					Math.max(prevOverlap, nextOverlap) + 1);
			postprocAbbr.put(key, i);
			prevOverlap = nextOverlap;
		}

		// LinkedList<Integer> list = new LinkedList<Integer>();
		// for (int i = 0; i < man.plugins.length; i++)
		// list.add(i);
		//
		// for (int len = 1; list.size() > 0; len++) {
		// for (ListIterator<Integer> it = list.listIterator(); it.hasNext();) {
		// int pp = it.next();
		// String str = man.plugins[pp].getTabName().substring(0, len)
		// .toUpperCase().replace(' ', '_');
		// Integer val;
		// if ((val = postprocAbbr.get(str)) == null) { // empty slot
		// postprocAbbr.put(str, pp);
		// it.remove(); // pp is done
		// } else if (val >= 0) { // first collision
		// postprocAbbr.put(str, -1);
		// // it.add(val); // pp is not done
		// }
		// }
		// }
		//
		// for (String key : postprocAbbr.keySet()) { // remove mappings for
		// // collisions
		// if (postprocAbbr.get(key) == -1)
		// postprocAbbr.remove(key);
		// }
	}

	private int checkOverlap(String str1, String str2) {
		int i = 0;
		while (i < str1.length() && i < str2.length()
				&& str1.charAt(i) == str2.charAt(i))
			i++;
		return i;
	}

	private void findSubstMods() {
		for (String model : Utils.classesInPackage(SubstitutionModel.class
				.getPackage().getName() + ".plugins")) {
			try {
				Class<?> cl = Class.forName(model);
				if (!Modifier.isAbstract(cl.getModifiers())
						&& SubstitutionModel.class.isAssignableFrom(cl))
					substModNames.add(model);
			} catch (Exception e) { // handle class access exceptions etc.
				// e.printStackTrace(System.err);
			}
		}
	}

	private String buildPpListStr(MainManager man, String linePrefix) {
		StringBuilder build = new StringBuilder();
		List<Postprocess> plugins = man.postProcMan.getPlugins();
		for (String key : postprocAbbr.keySet()) {
			build.append(linePrefix);
			build.append(key);
			build.append(": ");
			build.append(plugins.get(postprocAbbr.get(key))
					.getTip());
			build.append("\n");
		}
		return build.toString();
	}

	private String buildDefPpList(MainManager man) {
		List<Postprocess> plugins = man.postProcMan.getPlugins();
		StringBuilder build = new StringBuilder();
		for (String key : postprocAbbr.keySet())
			if (plugins.get(postprocAbbr.get(key)).sampling) {
				build.append(key);
				build.append(",");
			}
		build.deleteCharAt(build.length() - 1);
		return build.toString();
	}

	public boolean isVerbose() {
		return verbose;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

}
