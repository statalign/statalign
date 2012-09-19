package statalign;

import java.io.IOException;
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
import statalign.base.Utils;
import statalign.io.RawSequences;
import statalign.io.input.plugins.FastaReader;
import statalign.model.subst.RecognitionError;
import statalign.model.subst.SubstitutionModel;
import statalign.model.subst.plugins.Dayhoff;
import statalign.model.subst.plugins.Kimura3;
import statalign.postprocess.PluginParameters;
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
	
	// Functions

	/** 
	 * Creates a <tt>CommandLine</tt> object with <tt>verbose</tt> set to false.
	 * @param isParallel does this instance belong to a parallel run?
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
		ArrayList<String> parsedArgs = new ArrayList<String>();
		for(int i = 0 ; i < args.length ; i++)
		{
			if(!args[i].startsWith("plugin:"))
			{
				System.out.println(args[i]);
				parsedArgs.add(args[i]);
			}
		}*/
		
		Options opt = new Options(args, Multiplicity.ZERO_OR_ONE, 1,
				Integer.MAX_VALUE);
		opt.addSet("run").addOption("subst", Separator.EQUALS)
				.addOption("mcmc", Separator.EQUALS)
				.addOption("seed", Separator.EQUALS)
				.addOption("ot", Separator.EQUALS)
				.addOption("log", Separator.EQUALS).addOption("plugin", Separator.COLON, Multiplicity.ZERO_OR_MORE);

		OptionSet set;
		if ((set = opt.getMatchingSet(false, false)) == null) {
			return usage(manager);
		}

		FastaReader f = new FastaReader();
		try {
			for (String inputFile : set.getData()) {
				try {
					RawSequences seq = f.read(inputFile);
					int errorNum = f.getErrors();
					if (errorNum > 0)
						warning(errorNum + " errors occured reading "
								+ inputFile + ", please check file format");
					if (verbose) {
						System.out.println("INFO: read " + seq.size()
								+ " sequences from " + inputFile);
					}
					manager.inputData.seqs.add(seq);
				} catch (IOException e) {
					return error("error reading input file " + inputFile);
				}
			}

			// TODO: Change? If dealing with merged datasets this does not make
			// any sense.
			manager.fullPath = set.getData().get(0);

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
							return error("Substitution model "
									+ modelName
									+ " does not accept the given input sequences!");
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

				if (pars.length != 3 && pars.length != 4) {
					return error("MCMC parameters not recognized: " + mcmcPars);
				}

				manager.inputData.pars.burnIn = parseValue(pars[0]);
				manager.inputData.pars.cycles = parseValue(pars[1]);
				manager.inputData.pars.sampRate = parseValue(pars[2]);

				if (pars.length == 4) {
					if (!isParallel) {
						return error("Unrecognized MCMC parameters for non-parallel version.");
					}

					manager.inputData.pars.swapRate = parseValue(pars[3]);
					if (manager.inputData.pars.swapRate < 0) {
						return error("MCMC parameter not recognized: "
								+ mcmcPars);
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

			if (set.isSet("log")) {
				String log = set.getOption("log").getResultValue(0);
				String[] keys = log.split(",");
				Postprocess[] pps = manager.postProcMan.plugins;

				for (Postprocess pp : pps)
					pp.sampling = false;

				for (String key : keys) {
					try {
						pps[postprocAbbr.get(key.toUpperCase())].sampling = true;
					} catch (NullPointerException e) {
						return error("Log file entry code list not recognised: "
								+ log);
					}
				}
			}

			// retrieve all parameters starting with plugin:
			OptionData plugins = set.getOption("plugin");
			ArrayList<String> argsVector = new ArrayList<String>();
			for(int i = 0 ; i < plugins.getResultCount() ; i++)
			{
				argsVector.add(plugins.getResultValue(i));
			}
			Postprocess.pluginParameters = new PluginParameters(argsVector);
			
			AutomateParameters.setAutomateBurnIn(false);
			AutomateParameters.setAutomateStepRate(false);
			AutomateParameters.setAutomateNumberOfSamples(false);

		} catch (Exception e) {
			e.printStackTrace();
			return error("Unknown error: " + e.getLocalizedMessage());
		}

		return 0;
	}

	private String getUsageString(MainManager man) {
		StringBuilder sb = new StringBuilder();

		sb.append("Usage:\n\n");
		if (isParallel) {
			// TODO: this.
			sb.append("    <depends> statalign.jar [options] seqfile1 [seqfile2 ...]\n\n\n");
		} else {
			sb.append("    java -Xmx512m -jar statalign.jar [options] seqfile1 [seqfile2 ...]\n\n\n");
		}

		sb.append("Description:\n\n");
		sb.append("    StatAlign can be used for Bayesian analysis of protein, DNA and RNA\n");
		sb.append("    sequences. Multiple alignments, phylogenetic trees and evolutionary\n");
		sb.append("    parameters are co-estimated in a Markov Chain Monte Carlo framework.\n");
		sb.append("    The input sequence files must be in Fasta format.\n\n\n");

		sb.append("Options:\n\n");

		sb.append("    -subst=MODEL\n");
		sb.append("        Lets you select from the present substitution models (see list below)\n");
		sb.append("        Default: " + Kimura3.class.getSimpleName()
				+ " (for DNA/RNA data), " + Dayhoff.class.getSimpleName()
				+ " (for protein data)\n\n");

		if (isParallel) {
			sb.append("    -mcmc=burn,cycl,samprate,swaprate\n");
			sb.append("        Sets MCMC parameters: burn-in, cycles after burn-in, sampling rate, swap rate.\n");
			sb.append("          Abbreviations k and m mean 1e3 and 1e6 factors.\n");
			sb.append("        Default: 10k,100k,1k,100\n\n");
		} else {
			sb.append("    -mcmc=burn,cycl,rate\n");
			sb.append("        Sets MCMC parameters: burn-in, cycles after burn-in, sampling rate.\n");
			sb.append("          Abbreviations k and m mean 1e3 and 1e6 factors.\n");
			sb.append("        Default: 10k,100k,1k\n\n");
		}

		sb.append("    -seed=value\n");
		sb.append("        Sets the random seed (same value will reproduce same results for\n");
		sb.append("          identical input and settings)\n");
		sb.append("        Default: 1\n\n");

		sb.append("    -ot=OUTTYPE\n");
		sb.append("        Sets output alignment type.\n");
		sb.append("          (One of: " + Utils.joinStrings(MainManager.alignmentTypes, ", ") + ")\n");
		sb.append("        Default: " + MainManager.alignmentTypes[0] + "\n\n");

		sb.append("    -log=["
				+ Utils.joinStrings(postprocAbbr.keySet().toArray(), "][,")
				+ "]\n");
		sb.append("        Lets you customise what is written into the log file (one entry\n");
		sb.append("        for each sample).\n");
		sb.append(buildPpListStr(man, "          "));
		sb.append("        Default: " + buildDefPpList(man) + "\n");

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

	private int usage(MainManager man) {
		System.out.println(getUsageString(man));
		System.out.println("\nList of available substitution models:");
		for (String model : substModNames) {
			try {
				System.out.println("\t" + Class.forName(model).getSimpleName());
			} catch (Exception e) {
			}
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
		Postprocess[] plugins = man.plugins;
		
		final String[] keys = new String[plugins.length];
		Integer[] sorted = new Integer[plugins.length];
		
		for(int i = 0; i < plugins.length; i++) {
			keys[i] = plugins[i].getTabName().toUpperCase().replace(' ', '_');
			sorted[i] = i;
		}
			
		Arrays.sort(sorted, new Comparator<Integer>() {
			@Override
			public int compare(Integer o1, Integer o2) {
				return keys[o1].compareTo(keys[o2]);
			}
		});
		
		int prevOverlap = 0;
		for(int i = 0; i < plugins.length; i++) {
			int nextOverlap = 0;
			if(i < plugins.length-1)
				nextOverlap = checkOverlap(keys[i], keys[i+1]);
			String key = keys[i].substring(0, Math.max(prevOverlap, nextOverlap)+1);
			postprocAbbr.put(key, i);
			prevOverlap = nextOverlap;
		}
		
//		LinkedList<Integer> list = new LinkedList<Integer>();
//		for (int i = 0; i < man.plugins.length; i++)
//			list.add(i);
//
//		for (int len = 1; list.size() > 0; len++) {
//			for (ListIterator<Integer> it = list.listIterator(); it.hasNext();) {
//				int pp = it.next();
//				String str = man.plugins[pp].getTabName().substring(0, len)
//						.toUpperCase().replace(' ', '_');
//				Integer val;
//				if ((val = postprocAbbr.get(str)) == null) { // empty slot
//					postprocAbbr.put(str, pp);
//					it.remove(); // pp is done
//				} else if (val >= 0) { // first collision
//					postprocAbbr.put(str, -1);
//					// it.add(val); // pp is not done
//				}
//			}
//		}
//
//		for (String key : postprocAbbr.keySet()) { // remove mappings for
//													// collisions
//			if (postprocAbbr.get(key) == -1)
//				postprocAbbr.remove(key);
//		}
	}

	private int checkOverlap(String str1, String str2) {
		int i = 0;
		while(i < str1.length() && i < str2.length() && str1.charAt(i) == str2.charAt(i))
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
		for (String key : postprocAbbr.keySet()) {
			build.append(linePrefix);
			build.append(key);
			build.append(": ");
			build.append(man.postProcMan.plugins[postprocAbbr.get(key)]
					.getTip());
			build.append("\n");
		}
		return build.toString();
	}

	private String buildDefPpList(MainManager man) {
		StringBuilder build = new StringBuilder();
		for (String key : postprocAbbr.keySet())
			if (man.postProcMan.plugins[postprocAbbr.get(key)].sampling) {
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
