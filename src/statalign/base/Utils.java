package statalign.base;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;

import statalign.model.subst.plugins.Dayhoff;
import statalign.ui.ErrorMessage;

/**
 * 
 * This class contains multi-purpose static functions.
 * 
 * @author miklos, novak, herman
 *
 */
public class Utils{

	private Utils(){
		//... so that no objects of this class can be created
	}

	/**
	 * Debugging mode (various consistency checks done if on)
	 */
	public static final boolean DEBUG = true;
	
	/**
	 * The random number generator used throughout the program.
	 * A new generator is constructed at each MCMC run using the seed in the
	 * corresponding MCMCPars object.
	 */
	public static Random generator = new Random(1);
	
	/**
	 * During the burnin, the SPAN variables for all continuous parameters
	 * are adjusted in order to ensure that the average acceptance rate is between 
	 * 0.1 and 0.4 where possible. This is done by repeatedly multiplying the SPAN
	 * by SPAN_MULTIPLIER until the acceptance falls within the desired range,
	 * specified by [MIN_ACCEPTANCE,MAX_ACCEPTANCE].
	 */
	public static final double SPAN_MULTIPLIER = 0.7;
	public static final double MIN_ACCEPTANCE = 0.15;
	public static final double MAX_ACCEPTANCE = 0.4;
	
	public static final double WINDOW_CHANGE_FACTOR = 0.8;
	public static final double MIN_WINDOW_MULTIPLIER = 0.1;
	public static final double MAX_WINDOW_MULTIPLIER = 2;
	public static double WINDOW_MULTIPLIER = 1.0;

	public static final double MIN_SAMPLES_FOR_ACC_ESTIMATE = 10;

	/**
	 * When a new edge length is proposed in MCMC, it is drawn uniformly from the
	 * neighborhood of the current value. The neighborhood size has this span.
	 */
	public static double EDGE_SPAN = 0.1;
	/**
	 * When a new R of the TKF92 model is proposed in MCMC, it is drawn uniformly from the
	 * a neighborhood of the current value. The neighborhood size has this span.
	 */
	public static double R_SPAN = 0.05;
	/**
	 * When a new lambda of the TKF92 model is proposed in MCMC, it is drawn uniformly from the
	 * a neighborhood of the current value. The neighborhood size has this span.
	 */
	public static double LAMBDA_SPAN = 0.01;
	/**
	 * When a new mu of the TKF92 model is proposed in MCMC, it is drawn uniformly from the
	 * a neighborhood of the current value. The neighborhood size has this span.
	 */
	public static double MU_SPAN = 0.01;
	/**
	 * log(0) is set to Double.NEGATIVE_INFINITY. This is used in logarithmic adding.
	 * The logarithm of an empty sum is set to this value.
	 * 
	 */
	public static final double log0 = Double.NEGATIVE_INFINITY;
	
	private static double[] tempDoubleArray;

	/**
	 * This function selects a random integer such that the expected squared value will be O(length)
	 * The probability of the selection of that particular index is returned in selectLike.
	 * 
	 * (MuDouble is used to allow for another return value, in C++ a double pointer/reference
	 * could be used instead)
	 * 
	 * @param length The length of the array we need.
	 * @param selectLike A mutable double object to return the selection probability
	 * @return A random integer as described above
	 */
	public static int linearizerWeight(int length, MuDouble selectLike, double expectedLength){
		if(tempDoubleArray == null || tempDoubleArray.length < length){
			tempDoubleArray = new double[length];
		}
		double p = 1.0 - 1.0/expectedLength;
		tempDoubleArray[0] = 1.0/expectedLength;
		double sum = tempDoubleArray[0];
		for(int i = 1; i < length; i++){
			tempDoubleArray[i] = tempDoubleArray[i-1] * p;
			sum += tempDoubleArray[i];
		}
		double w = generator.nextDouble() * sum;
		int k = 0;
		double x = tempDoubleArray[0];
		while(x < w && k < length - 1){
			k++;
			x += tempDoubleArray[k];
		}

		//	System.out.println((tempDoubleArray[k]/sum));
		selectLike.value = (tempDoubleArray[k]/sum); 

		return k;
	}

	/**
	 * 
	 * This function returns the probability of choosing a particular index with linearizerWeight.
	 * The value returned is equal to mu.value when linearizerWeight(length, mu) returns 'index'.
	 * 
	 * @param length Distribution parameter as in linearizerWeight
	 * @param index Selected index
	 * @return Probability of the selection
	 */
	public static double linearizerWeightProb(int length, int index, double expectedLength){
		if(tempDoubleArray == null || tempDoubleArray.length < length){
			tempDoubleArray = new double[length];
		}
		double p = 1.0 - 1.0/expectedLength;
		tempDoubleArray[0] = 1.0/expectedLength;
		double sum = tempDoubleArray[0];
		for(int i = 1; i < length; i++){
			tempDoubleArray[i] = tempDoubleArray[i-1] * p;
			sum += tempDoubleArray[i];
		}

		return tempDoubleArray[index]/sum;
	}

	/**
	 * This function returns a random index, weighted by the weights in the array `weights'
	 */
	public static int weightedChoose(int[] weights){
		int sum = 0;

		for(int i = 0; i < weights.length; i++){
			sum += weights[i];
		}

		int w = generator.nextInt(sum);
		int k = 0;
		sum = 0;
		while((sum += weights[k]) <= w){
			k++;
		}

		return k;
	}

	/**
	 * Similar to weightedChoose(weights), but the log-probability of the selection
	 * will be subtracted from the mutable double object selectLogLike
	 * (reason: proposal is in the denominator of acceptance ratio)
	 * 
	 * (MuDouble is used to allow for another return value, in C++ a double pointer/reference
	 * could be used instead)
	 * 
	 */
	public static int weightedChoose(double[] weights, MuDouble selectLogLike){
		double sum = 0.0;

		for(int i = 0; i < weights.length; i++){
			sum += weights[i];
		}
		if(selectLogLike != null)
			selectLogLike.value += Math.log(sum);

		double w = generator.nextDouble() * sum;
		int k = 0;
		sum = 0.0;
		while(k < weights.length-1 && (sum += weights[k]) <= w){
			k++;
		}
		if(selectLogLike != null)
			selectLogLike.value -= Math.log(weights[k]);

		assert (weights[k] > 1e-5) : "weightedChoose error";

		return k;
	}

	/**
	 * Behaves exactly like weightedChoose(new double[]{1-prob,prob}, selectLogLike), but faster
	 */
	public static int chooseOne(double prob, MuDouble selectLogLike) {
		if(generator.nextDouble() < prob) {
			selectLogLike.value -= Math.log(prob);
			return 1;
		}
		selectLogLike.value -= Math.log(1.0-prob);
		return 0;
	}

	/**
	 * Equivalent to weightedChoose(weights, selectLogLike) where
	 * logWeights[i] = Math.log(weights[i]), but avoids overflows that might result from
	 * exponentiation. 

	 * (MuDouble is used to allow for another return value, in C++ a double pointer/reference
	 * could be used instead)
	 */
	public static int logWeightedChoose(double[] logWeights, MuDouble selectLogLike){
		double logSum = log0;

		for(int i = 0; i < logWeights.length; i++){
			logSum = logAdd(logSum, logWeights[i]);
		}
		if(selectLogLike != null)
			selectLogLike.value += logSum;

		double w = Math.log(generator.nextDouble()) + logSum;
		int k = 0;
		logSum = log0;
		while(k < logWeights.length-1 && (logSum = logAdd(logSum, logWeights[k])) <= w){
			k++;
		}
		if(selectLogLike != null)
			selectLogLike.value -= logWeights[k];

		assert (logWeights[k] > log0) : "logWeightedChoose error";

		return k;
	}

	/**
	 * Logarithmically add two numbers
	 * @param a log(x)
	 * @param b log(y)
	 * @return log(x+y)
	 */
	public static double logAdd(double a, double b) {
		if(a == b)
			return Math.log(2)+a;
		if(a < b)
			return b+Math.log(Math.exp(a-b)+1);
		return a+Math.log(Math.exp(b-a)+1);
	}

    static void calcFelsen(double[] res, double[] fel1, double[][] prob1, double[] fel2, double[][] prob2) {
		double s;
		int i, j, len = res.length;

		for(i = 0; i < len; i++) {
			if(fel1 != null) {
				s = 0.0;
				for(j = 0; j < len; j++)
					s += prob1[i][j]*fel1[j];
				res[i] = s;
			} 
			else{
				res[i] = 1.0;
			}
			if(fel2 != null) {
				s = 0.0;
				for(j = 0; j < len; j++){
					s += prob2[i][j]*fel2[j];
				}
				res[i] *= s;
			}
		}
	}

   static boolean calcFelsenWithCheck(double[] res, double[] fel1, double[][] prob1, double[] fel2, double[][] prob2) {
		double s;
		int i, j, len = res.length;

		double[] oldres = new double[res.length];
		for(i = 0; i < res.length; i++){
			oldres[i] = res[i];
		}

		for(i = 0; i < len; i++) {
			if(fel1 != null) {
				s = 0.0;
				for(j = 0; j < len; j++)
					s += prob1[i][j]*fel1[j];
				res[i] = s;
			} 
			else{
				res[i] = 1.0;
			}
			if(fel2 != null) {
				s = 0.0;
				for(j = 0; j < len; j++){
					s += prob2[i][j]*fel2[j];
				}
				res[i] *= s;
			}
		}
		boolean match = true;
		for(i = 0; i < res.length && match; i++){
			match = Math.abs(oldres[i]/res[i] - 1.0) < 0.0001;
		}
		return match;
	}

	/**
	 *  Calculates emission probability from Felsenstein likelihoods
	 */
   static double calcEmProb(double fel[], double aaEquDist[]) {
		double p = 0;

		for(int i = 0; i < fel.length; i++)
			p += fel[i]*aaEquDist[i];
		return p;
	}
	
	/**
	 * Makes Enumeration iterable.
	 * 
	 * @param <T> Enumeration element type
	 * @param en  the Enumeration
	 * @return an Iterable that can iterate through the elements of the Enumeration
	 */
	public static <T> Iterable<T> iterate(final Enumeration<T> en) {
		final Iterator<T> iterator = new Iterator<T>() {
			@Override
			public boolean hasNext() {
				return en.hasMoreElements();
			}
			@Override
			public T next() {
				return en.nextElement();  
			}
			@Override
			public void remove() {
				throw new UnsupportedOperationException();
			}
		};
		return new Iterable<T>() {
			@Override
			public Iterator<T> iterator() {
				return iterator;
			}
		};
	}
	
	/**
	 * Joins strings using a separator string. Accepts any <code>Object</code>s
	 * converting them to strings using their <code>toString</code> method.
	 * 
	 * @param strs
	 *            strings to join
	 * @param separator
	 *            the separator string
	 * @return a string made up of the strings separated by the separator
	 */
	public static String joinStrings(Object[] strs, String separator) {
		StringBuilder result = new StringBuilder();
		int i, l = strs.length;
		if (l > 0) {
			result.append(strs[0]);
		}
		for (i = 1; i < l; i++) {
			result.append(separator);
			result.append(strs[i]);
		}
		return result.toString();
	}

	/**
	 * Joins strings using a prefix and a separator string. Accepts any
	 * <code>Object</code>s converting them to strings using their
	 * <code>toString</code> method.
	 * 
	 * @param strs
	 *            strings to join
	 * @param prefix
	 *            prefix for each string
	 * @param separator
	 *            the separator string
	 * @return a string made up of the strings with the given prefix and
	 *         separated by the separator
	 */
	public static String joinStrings(Object[] strs, String prefix, String separator) {
		StringBuilder result = new StringBuilder();
		int i, l = strs.length;
		if (l > 0) {
			result.append(prefix);
			result.append(strs[0]);
		}
		for (i = 1; i < l; i++) {
			result.append(separator);
			result.append(prefix);
			result.append(strs[i]);
		}
		return result.toString();
	}

	/**
	 * Finds all classes in a given package and all of its subpackages by walking
	 * through class path. Handles both directories and jar files.
	 * 
	 * @param packageName  the package in which the classes are searched for
	 * @return array of found class names (with full package prefixes)
	 */
	public static String[] classesInPackage(String packageName) {
		ArrayList<String> classNames = new ArrayList<String>();
		String packageDir = packageName.replace('.', '/');

		CircularArray<File> dirs = new CircularArray<File>();
		CircularArray<String> pkgPrefs = new CircularArray<String>();
		ClassLoader.getSystemClassLoader();
//		ErrorMessage.showPane(null, System.getProperty("java.class.path"), false);
        // TODO: FIX the classpath stuff.
//		String[] paths = new String[]{ System.getProperty("user.dir") + "/bin" };
		String[] paths = System.getProperty("java.class.path").split(File.pathSeparator);
		for(String path : paths) {
			try {
				File file = new File(path);

				if(!file.isDirectory()) {
					JarFile jarFile = new JarFile(file);

					Enumeration<JarEntry> fileNames = jarFile.entries();
					while(fileNames.hasMoreElements()) {
						String entryName = fileNames.nextElement().getName();
						if(entryName.startsWith(packageDir+"/") && entryName.endsWith(".class")) {
							classNames.add(entryName.substring(0, entryName.length()-6).replace('/', '.'));
						}
					}
				} else {
					dirs.push(new File(path+"/"+packageDir));
					pkgPrefs.push(packageName);
					while((file = dirs.shift()) != null) {
						File[] children = file.listFiles();
						String pkgPrefix = pkgPrefs.shift()+".";

						if(children != null) {
							for(File child : children) {
								String chName = child.getName();
								if(child.isDirectory()) {
									dirs.push(child);
									pkgPrefs.push(pkgPrefix+chName);
								} else if(chName.endsWith(".class")) {
									classNames.add(pkgPrefix+chName.substring(0, chName.length()-6));
								}
							}
						} // if children
					} // for dirs
				} // if !cpFile.isDir
			} catch(IOException e) {
				ErrorMessage.showPane(null, e, true);
			}
		} // for classpath
		
		return classNames.toArray(new String[classNames.size()]);
	}
	
	/**
	 * Locates all plugins that are descendants of the specified plugin superclass. The plugins are
	 * expected to be in the package <code><i>root</i>.plugins</code> where <code><i>root</i></code> refers
	 * to the package of the superclass.
	 * @param superClass the ancestral plugin class
	 * @return list of plugins found
	 */
	@SuppressWarnings("unchecked")
	public static<T> List<T> findPlugins(Class<T> superClass) {
		String[] pluginNames = Utils.classesInPackage(superClass.getPackage().getName()+".plugins");
		List<T> pluginList = new ArrayList<T>();
		for(int i = 0; i < pluginNames.length; i++) {
			try {
				Class<?> cl = Class.forName(pluginNames[i]);
				if(!superClass.isAssignableFrom(cl))
					continue;
				pluginList.add((T)cl.newInstance());
			} catch (Exception e) {
			}
		}
		return pluginList;
	}

	/**
	 * Transforms an alignment into the prescribed format
	 * 
	 * @param s String array containing the alignment in StatAlign format
	 * @param type The name of the format, might be "StatAlign", "Clustal", "Fasta", "Phylip", "Nexus"
	 * @param input The input data. Needed for the Nexus format that needs
	 *             a name of the alignment (set to input.title) the type of the alignment (either
	 *             nucleotide or protein, read from input.model) and the list of characters
	 *             in the substitution model (also read from input.model).
	 * @return String array containing the alignment in the prescribed format
	 */
   public static String[] alignmentTransformation(String[] s, String type, InputData input){
		if(type.equals("StatAlign")){
			return s;
		}
		String[] returnAlignment = null;
		//here we need different alignment types
		String[] names = new String[s.length];
		String[] alignment = new String[s.length];
		int leavesNumber = 0;
		for(int i = 0; i < names.length; i++){
			names[i] = s[i].substring(0,s[i].indexOf('\t'));
			if(names[i].charAt(0) != ' '){
				leavesNumber++;
			}
			alignment[i] = s[i].substring(s[i].indexOf('\t')+1,s[i].length());
		}
		String[] onlyLeavesAlignment = new String[leavesNumber];
		String[] newNames =  new String[leavesNumber];
		leavesNumber = 0;
		for(int i = 0; i < names.length; i++){
			if(names[i].charAt(0) != ' '){
				onlyLeavesAlignment[leavesNumber] = s[i].substring(s[i].indexOf('\t')+1);
				newNames[leavesNumber] = names[i];
				newNames[leavesNumber] = newNames[leavesNumber].replaceAll(" ", "");
//				newNames[leavesNumber] = newNames[leavesNumber].replaceAll("_", " ");
				newNames[leavesNumber] = newNames[leavesNumber].replaceAll("\\{", "(");
				newNames[leavesNumber] = newNames[leavesNumber].replaceAll("\\}", ")");
				leavesNumber++;
			}
		}
		onlyLeavesAlignment = deleteAllGaps(onlyLeavesAlignment);
		if(type.equals("Clustal")){
			for(int i = 0; i < newNames.length; i++){
				newNames[i] = newNames[i].replaceAll("\\(", "_");
				newNames[i] = newNames[i].replaceAll("\\)", "_");
				if(newNames[i].indexOf(' ') == -1){
					newNames[i] = newNames[i].substring(0,Math.min(newNames[i].length(),30));
					while(newNames[i].length() < 30){
						newNames[i] += " ";
					}
				}
				else{
					newNames[i] = newNames[i].substring(0,Math.min(30, newNames[i].indexOf(' ')));
					while(newNames[i].length() < 30){
						newNames[i] += " ";
					}				
				}

			}
			//	System.out.println("Making Fasta alignment");
			returnAlignment = new String[(onlyLeavesAlignment.length + 3) *(onlyLeavesAlignment[0].length()/60 + (onlyLeavesAlignment[0].length() % 60 == 0 ? 0 : 1)) + 1];
		//	System.out.println("only leaves number: "+onlyLeavesAlignment.length+
			//				   "alignment length: "+onlyLeavesAlignment[0].length()+
				//			   	"new array length: "+returnAlignment.length);
			int j = 1;
			returnAlignment[0] = "CLUSTAL X type multiple sequence alignment generated by StatAlign";
			for(int i = 0; i < onlyLeavesAlignment[0].length()/60; i++){
				returnAlignment[j] = "";
				returnAlignment[j+1] = "";
				j+=2;
				//System.out.println(returnAlignment[j]);
				for(int k = 0; k < newNames.length; k++){
					returnAlignment[j] = newNames[k]+"\t"+onlyLeavesAlignment[k].substring(i*60, (i+1) * 60);
					//System.out.println(returnAlignment[j]);
					j++;
				}
				returnAlignment[j] = "                              \t";
				for(int k = i*60; k < (i+1) * 60; k++){
					returnAlignment[j] += clustalCharacter(onlyLeavesAlignment,k);
				}
				j++;
			}
			//last line, if needed
			if(onlyLeavesAlignment[0].length() % 60 != 0){
				returnAlignment[j] = "";
				returnAlignment[j+1] = "";
				j+=2;
				//System.out.println(returnAlignment[j]);
				for(int k = 0; k < newNames.length; k++){
					returnAlignment[j] = newNames[k]+"\t"+onlyLeavesAlignment[k].substring((onlyLeavesAlignment[k].length()/60)*60);
					//System.out.println(returnAlignment[j]);
					j++;
				}
				returnAlignment[j] = "                              \t";
				for(int k = (onlyLeavesAlignment[0].length()/60)*60; k < onlyLeavesAlignment[0].length(); k++){
					returnAlignment[j] += clustalCharacter(onlyLeavesAlignment,k);
				}
			}
		}
		else if(type.equals("Fasta")){
		//	System.out.println("Making Fasta alignment");
			returnAlignment = new String[onlyLeavesAlignment.length *(onlyLeavesAlignment[0].length()/60 + (onlyLeavesAlignment[0].length() % 60 == 0 ? 1 : 2))];
		//	System.out.println("only leaves number: "+onlyLeavesAlignment.length+
			//				   "alignment length: "+onlyLeavesAlignment[0].length()+
				//			   	"new array length: "+returnAlignment.length);
			int j = 0;
			for(int i = 0; i < newNames.length; i++){
				returnAlignment[j] = ">"+newNames[i]; 
				//System.out.println(returnAlignment[j]);
				j++;
				for(int k = 0; k < onlyLeavesAlignment[i].length()/60; k++){
					returnAlignment[j] = onlyLeavesAlignment[i].substring(k*60, (k+1) * 60);
					//System.out.println(returnAlignment[j]);
					j++;
				}
				if(onlyLeavesAlignment[0].length() % 60 != 0){
					returnAlignment[j] = onlyLeavesAlignment[i].substring((onlyLeavesAlignment[i].length()/60)*60);
					//System.out.println(returnAlignment[j]);
					j++;
				}
			}
		}
		else if(type.equals("Phylip")){
			returnAlignment = new String[(onlyLeavesAlignment.length + 1)*(onlyLeavesAlignment[0].length()/60 + (onlyLeavesAlignment[0].length() % 60 == 0 ? 0 : 1))];
			returnAlignment[0] = " "+onlyLeavesAlignment.length+" "+onlyLeavesAlignment[0].length();
			for(int i = 0; i < newNames.length; i++){
				newNames[i] = newNames[i].replaceAll("\\(", "_");
				newNames[i] = newNames[i].replaceAll("\\)", "_");
				newNames[i] = newNames[i].substring(0,Math.min(newNames[i].length(),10));
				while(newNames[i].length() < 13){
					newNames[i] += " ";
				}
			}
			int j = 1;
			for(int k = 0; k < newNames.length; k++){
				returnAlignment[j] = newNames[k];
				for(int l = 0; l < 6; l++){
					returnAlignment[j] += (onlyLeavesAlignment[k].length() >= (l+1)*10 ? 
										onlyLeavesAlignment[k].substring(l*10,(l+1)*10)+" " :
											"");
				}
				j++;
			}
			for(int i = 1; i < onlyLeavesAlignment[0].length()/60; i++){
				returnAlignment[j] = "";
				j++;
				for(int k = 0; k < newNames.length; k++){
					returnAlignment[j] = "             ";
					for(int l = 0; l < 6; l++){
						returnAlignment[j] += onlyLeavesAlignment[k].substring(i*60+l*10,i*60+(l+1)*10)+" ";
					}
					j++;
				}
				
			}
			//last line, if needed
			if(onlyLeavesAlignment[0].length() % 60 != 0 && onlyLeavesAlignment[0].length() > 60){
				returnAlignment[j] = "";
				j++;
				for(int k = 0; k < newNames.length; k++){
					returnAlignment[j] = "             ";
					int l;
					for(l = 0; l < (onlyLeavesAlignment[0].length() - (onlyLeavesAlignment[0].length()/60)*60)/10; l++){
						returnAlignment[j] += onlyLeavesAlignment[k].substring((onlyLeavesAlignment[0].length()/60)*60+l*10,
								(onlyLeavesAlignment[0].length()/60)*60+(l+1)*10)+" ";
					}
					//and the last piece...
					returnAlignment[j] += onlyLeavesAlignment[k].substring((onlyLeavesAlignment[0].length()/60)*60+l*10);
					j++;
				}
			}
		}
		else if(type.equals("Nexus")){
			returnAlignment = new String[(onlyLeavesAlignment.length + 1)*((onlyLeavesAlignment[0].length()/50) + (onlyLeavesAlignment[0].length() % 50 == 0 ? 0 : 1)) + 10];
			returnAlignment[0] = "#NEXUS";
			returnAlignment[1] = "[TITLE: "+input.title+"]";
			returnAlignment[2] = "";
			returnAlignment[3] = "begin data;";
			returnAlignment[4] = "dimensions ntax="+onlyLeavesAlignment.length+" nchar="+onlyLeavesAlignment[0].length();
			returnAlignment[5] = "format interleave datatype="+input.model.getType()+"   gap=- symbols=\"";
			for(int i = 0; i < input.model.alphabet.length; i++){
				returnAlignment[5] += input.model.alphabet[i];
			}
			returnAlignment[5]+="\";";
			returnAlignment[6] = "";
			returnAlignment[7] = "matrix";
			for(int i = 0; i < newNames.length; i++){
				newNames[i] = newNames[i].replaceAll("\\(", "_");
				newNames[i] = newNames[i].replaceAll("\\)", "_");
				newNames[i] = newNames[i].substring(0,Math.min(newNames[i].length(),24));
				while(newNames[i].length() < 25){
					newNames[i] += " ";
				}
			}
			int j = 8;
			for(int i = 0; i < onlyLeavesAlignment[0].length()/50; i++){
				//System.out.println(returnAlignment[j]);
				for(int k = 0; k < newNames.length; k++){
					returnAlignment[j] = newNames[k]+onlyLeavesAlignment[k].substring(i*50, (i+1) * 50);
					//System.out.println(returnAlignment[j]);
					j++;
				}
				returnAlignment[j] = "";
				j++;
			}
			if(onlyLeavesAlignment[0].length() % 50 != 0){
				for(int k = 0; k < newNames.length; k++){
					returnAlignment[j] = newNames[k]+onlyLeavesAlignment[k].substring((onlyLeavesAlignment[k].length()/50)*50);
					//System.out.println(returnAlignment[j]);
					j++;
				}
			}
			returnAlignment[j] = ";";
			returnAlignment[j+1] = "end;";
		
		}
		return returnAlignment;
	}
	
	private static char clustalCharacter(String[] alignment, int pos){
		char c = '*';
		for(int i = 0; i < alignment.length - 1; i++){
			if(alignment[i].charAt(pos) != alignment[i+1].charAt(pos)){
				c = ' ';
				break;
			}
		}
		return c;
	}
	
	private static String[] deleteAllGaps(String[] alignment){
		String[] gapFreeAlignment = new String[alignment.length];
		for(int i = 0; i < gapFreeAlignment.length; i++){
			gapFreeAlignment[i] = "";
		}
		for(int i = 0; i < alignment[0].length(); i++){
			if(!allSpace(alignment,i)){
				for(int j = 0; j < alignment.length; j++){
					gapFreeAlignment[j] += alignment[j].charAt(i);
				}
			}
		}
		return gapFreeAlignment;
	}
	
   static boolean allSpace(String[] s, int p){
		boolean b = true;
		for(int i = 0; i < s.length && b; i++){
			b = s[i].charAt(p) == ' ';
		}
		return b;
	}
   
   public static char[] copyOf(char[] array) {
		int len = array.length;
		char[] copy = new char[len];
		System.arraycopy(array, 0, copy, 0, len);
		return copy;
   }
   
	public static int[] copyOf(int[] array) {
		int len = array.length;
		int[] copy = new int[len];
		System.arraycopy(array, 0, copy, 0, len);
		return copy;
	}
	
	public static double[] copyOf(double[] array) {
		int len = array.length;
		double[] copy = new double[len];
		System.arraycopy(array, 0, copy, 0, len);
		return copy;
	}


	/**
	 * Merely for testing purposes....
	 * 
	 * @param args Arguments are not used.
	 */
	public static void main(String[] args){
		for(String s : Utils.classesInPackage(Dayhoff.class.getPackage().getName()))
			System.out.println(s);
//		double[] w = new double[] {1,0,3,4,5};
//		for(int i = 0; i < 10000; i++){
//			MuDouble x = new MuDouble();
//			int index = linearizerWeight(1000,x);
//			System.out.println("I chose index "+index+" probability: "+x.value+" backprobability: "+
//					linearizerWeightProb(1000,index));
//		}
//		for(int i = 0; i < 1000; i++){
//			Double x = new Double(0.0);
//
//		}
	}
}

/**
 * Mutable Double class.
 * 
 * A MuDouble object can be used in place of a double pointer/reference that is missing
 * from the Java language.
 */
class MuDouble {
	public double value;
	
	public MuDouble() {
	}
	
	public MuDouble(double value) {
		this.value = value;
	}
	
	@Override
	public String toString() {
		return Double.toString(value);
	}
}