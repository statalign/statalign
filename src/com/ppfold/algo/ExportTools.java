package com.ppfold.algo;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Stack;

/**
 * Static methods for exporting output to various formats.
 * 
 * @author Z.Sukosd
 */

public class ExportTools {

	public static int writeTree(String dir_fil, String filename, String extension, Tree tree) throws IOException{

		new File(dir_fil).mkdir();
		PrintWriter out;
		String newname = filename.replaceAll("/", "-");
		newname = newname.replaceAll("\\*", "-");
	
		try {
			out = new PrintWriter(new FileWriter(dir_fil + "/" + newname
					+ extension));
			writeChildrenNodeToFile(tree.getRoot(),null,out);
			out.write(";\n");
			out.close();
			System.out.println("Written file " + dir_fil + "/" + newname
					+ extension);
			return 0;
		}
		catch(IOException e){
			System.err.println("Error writing file: " + dir_fil + "/" + newname
					+ extension);
			throw new IOException(e);
		}
	}
	
	private static void writeChildrenNodeToFile(Node thisnode, Node parent, PrintWriter out){
		//Writes tree in Newick format 
		if(thisnode.getChildren().size()>0){
			out.write("(");
		}
		int cnt = 0;
		for(Node child:thisnode.getChildren()){
			writeChildrenNodeToFile(child, thisnode, out);
			if(cnt<thisnode.getChildren().size()-1){
				//if there are more children
				out.write(",");
			}
			cnt++;
		}
		if(thisnode.getChildren().size()>0){
			out.write(")");
		}
		if(thisnode.getName()!=null){
			//a node with name
			out.write(thisnode.getName());
		}
		if(thisnode.getDistanceFromParent()!=0){
			out.format(":%.4f",thisnode.getDistanceFromParent());
		}
	}
	
	public static int writeTabbedMatrix(String dir_fil, String filename,
			String extension, float[][] matrix) throws IOException {

		new File(dir_fil).mkdir();

		PrintWriter out;
		String newname = filename.replaceAll("/", "-");
		newname = newname.replaceAll("\\*", "-");

		try {
			out = new PrintWriter(new FileWriter(dir_fil + "/" + newname
					+ extension));
			int i;
			for (i = 0; i < matrix.length; i++) {
				float[] row = matrix[i];
				boolean dlm = false;
				int j;
				for (j = 0; j < row.length; j++) {
					if (dlm) {
						out.print('\t');
					} else {
						// out.print(i + "\t"); //prints line number
						dlm = true;
					}
					out.format("%7.3e", row[j]);
				}
				out.println();
			}
			out.close();
			System.out.println("Written file " + dir_fil + "/" + newname
					+ extension);
			return 0;
		} catch (IOException e) {
			System.err.println("Error writing file: " + dir_fil + "/" + newname
					+ extension);
			throw new IOException(e);
		}
	}
	
	
	public static int writeTabbedMatrixToStream(byte[][] matrix, BufferedOutputStream out) throws IOException {

		try {
			int i;
			for (i = 0; i < matrix.length; i++) {
				byte[] row = matrix[i];
				boolean dlm = false;
				int j;
				for (j = 0; j < row.length; j++) {
					if (dlm) {
						out.write("\t".getBytes());
					} else {
						// out.print(i + "\t"); //prints line number
						dlm = true;
					}
					out.write((String.format("%7.2e", (float)row[j]/100)).getBytes());
				}
				out.write("\n".getBytes());
			}
			out.close();
			System.out.println("Written basepair plot");
			return 0;
		} catch (IOException e) {
			System.err.println("Error writing basepair plot");
			throw new IOException(e);
		}
	}

	public static int writeStructureReliability(String dir_fil,
			String filename, String extension, char[] structure,
			float[] reliability) throws IOException {
		new File(dir_fil).mkdir();

		PrintWriter out;
		String newname = filename.replaceAll("/", "-");
		newname = newname.replaceAll("\\*", "-");

		try {
			out = new PrintWriter(new FileWriter(dir_fil + "/" + newname
					+ extension));
			for (int c = 0; c < structure.length; c++) {
				out.format((c + 1) + "\t " + structure[c] + "\t %.4f \n",
						reliability[c]);
			}
			out.close();
			System.out.println("Written file " + dir_fil + "/" + newname
					+ extension);
			return 0;
		} catch (IOException e) {
			System.err.println("Error writing file: " + dir_fil + "/" + newname
					+ extension);
			throw new IOException(e);
		}

	}
	
	public static int writeStructureReliability_toStream(char[] structure,
			float[] reliability, BufferedOutputStream out) throws IOException {
		try {
			for (int c = 0; c < structure.length; c++) {
				out.write(String.format((c + 1) + "\t " + structure[c] + "\t %.4f \n",
						reliability[c]).getBytes());
			}
			out.close();
			System.out.println("Written .st file");
			return 0;
		} catch (IOException e) {
			System.err.println("Error writing .st file");
			throw new IOException(e);
		}

	}
	
	public static int writeFastaStructure(String dir_fil, String filename,
			String extension, char[] structure, List<char[]> columns,
			List<String> names) throws IOException {
		new File(dir_fil).mkdir();

		PrintWriter out;
		String newname = filename.replaceAll("/", "-");
		newname = newname.replaceAll("\\*", "-");

		try {
			out = new PrintWriter(new FileWriter(dir_fil + "/" + newname
					+ extension));
			out.print("pairingmask " + String.valueOf(structure));
			for (int seq = 0; seq < columns.get(0).length; seq++) {
				out.print("\n");
				String thisline = "";
				for (int c = 0; c < columns.size(); c++) {
					thisline = thisline
							.concat((structure[c] == '(' || structure[c] == ')') ? String
									.valueOf(columns.get(c)[seq]).toUpperCase()
									: String.valueOf(columns.get(c)[seq])
											.toLowerCase());
				}
				out.print(names.get(seq) + " " + thisline);
			}
			out.close();
			System.out.println("Written file " + dir_fil + "/" + newname
					+ extension);
			return 0;
		} catch (IOException e) {
			System.err.println("Error writing file: " + dir_fil + "/" + newname
					+ extension);
			throw new IOException(e);
		}

	}
	
	public static int writeFastaStructure_toStream(char[] structure, List<char[]> columns,
			List<String> names, BufferedOutputStream out) throws IOException {

		try {
			out.write(("pairingmask " + String.valueOf(structure)).getBytes());
			for (int seq = 0; seq < columns.get(0).length; seq++) {
				out.write("\n".getBytes());
				String thisline = "";
				for (int c = 0; c < columns.size(); c++) {
					thisline = thisline
							.concat((structure[c] == '(' || structure[c] == ')') ? String
									.valueOf(columns.get(c)[seq]).toUpperCase()
									: String.valueOf(columns.get(c)[seq])
											.toLowerCase());
				}
				out.write((names.get(seq) + " " + thisline).getBytes());
			}
			out.close();
			System.out.println("Written seq format");
			return 0;
		} catch (IOException e) {
			System.err.println("Error writing seq format");
			throw new IOException(e);
		}
	}
	

	public static int writeFullFastaStructure(String dir_fil, String filename,
			String extension, char[] structure, List<char[]> columns,
			List<String> names) throws IOException {
		// will write following format:
		// consensus pairing mask
		// sequence 1 sequence
		// sequence 1 structure
		// sequence 2 sequence
		// sequence 2 structure
		// etc.
		new File(dir_fil).mkdir();

		PrintWriter out;
		String newname = filename.replaceAll("/", "-");
		newname = newname.replaceAll("\\*", "-");

		try {
			out = new PrintWriter(new FileWriter(dir_fil + "/" + newname
					+ extension));
			out.print("consensus pairingmask " + String.valueOf(structure)
					+ "\n");
			for (int seq = 0; seq < columns.get(0).length; seq++) {
				out.print("\n");
				String thisline = "";
				for (int c = 0; c < columns.size(); c++) {
					thisline = thisline
				//			.concat((structure[c] == '(' || structure[c] == ')') ? String
				//					.valueOf(columns.get(c)[seq])
				//					: String.valueOf(columns.get(c)[seq]));
							.concat(String.valueOf(columns.get(c)[seq]));
				}
				//System.out.println(names.get(seq));
				//this is where one can do stem extension if required
				out.print(names.get(seq) + " " + reducesequence(thisline)
						+ "\n");
				out.print("pairingmask " + reducestructure(thisline, structure)
						+ "\n");
			}
			out.close();
			System.out.println("Written file " + dir_fil + "/" + newname
					+ extension);
			return 0;
		} catch (IOException e) {
			System.err.println("Error writing file: " + dir_fil + "/" + newname
					+ extension);
			throw new IOException(e);
		}

	}

	public static int writeFullFastaStructure_toStream(char[] structure, List<char[]> columns,
			List<String> names, BufferedOutputStream out) throws IOException {
		// will write following format:
		// consensus pairing mask
		// sequence 1 sequence
		// sequence 1 structure
		// sequence 2 sequence
		// sequence 2 structure
		// etc.

		try {
			out.write(("consensus pairingmask " + String.valueOf(structure)
					+ "\n").getBytes());
			for (int seq = 0; seq < columns.get(0).length; seq++) {
				out.write(("\n").getBytes());
				String thisline = "";
				for (int c = 0; c < columns.size(); c++) {
					thisline = thisline
							.concat((structure[c] == '(' || structure[c] == ')') ? String
									.valueOf(columns.get(c)[seq])
									: String.valueOf(columns.get(c)[seq]));
				}
				System.out.println(names.get(seq));
				//this is where one can do stem extension if required
				out.write((names.get(seq) + " " + reducesequence(thisline)
						+ "\n").getBytes());
				out.write(("pairingmask " + reducestructure(thisline, structure)
						+ "\n").getBytes());
			}
			out.close();
			System.out.println("Written .seq file");
			return 0;
		} catch (IOException e) {
			System.err.println("Error writing .seq file");
			throw new IOException(e);
		}

	}

	public static float[] reducereliabilities(String gappedsequencein, float[] gappedprobabilities) {
		int cnt = 0;
		char [] gappedsequence = gappedsequencein.toCharArray();
		for (char c : gappedsequence) {
			if (!MatrixTools.isGap(c)) {
				cnt++;
			}
		}
		int pos = 0;
		float[] reducedprobabilities = new float[cnt];
		for (int i = 1; i<gappedsequence.length; i++) {
			char c = gappedsequence[i];
			if (!MatrixTools.isGap(c)) {
				reducedprobabilities[pos] = gappedprobabilities[i];
				pos++;
			}
		}
		return reducedprobabilities;
	}
	
	
	
	public static String reducesequence(String gappedsequencein) {
		char[] gappedsequence = gappedsequencein.toCharArray();
		int cnt = 0;
		for (char c : gappedsequence) {
			if (!MatrixTools.isGap(c)) {
				cnt++;
			}
		}
		int pos = 0;
		char[] reducedsequence = new char[cnt];
		for (char c : gappedsequence) {
			if (!MatrixTools.isGap(c)) {
				reducedsequence[pos] = c;
				pos++;
			}
		}
		return new String(reducedsequence);
	}

	public static String reducestructure(String gappedsequencein,
			char[] gappedstructurein) {
		char[] gappedsequence = gappedsequencein.toCharArray();
		char[] gappedstructure = gappedstructurein.clone(); //make a copy so we don't affect the original!
		int cnt = 0;
		for (char c : gappedsequence) {
			if(!MatrixTools.isGap(c)) {
				cnt++;
			}
		}
		
		//System.out.println(gappedsequencein.length()
		//		+ " is length of input sequence ");
		//System.out.println(cnt + " letters in total ");

		// create pairing array for gapped structure
		Stack<Integer> s = new Stack<Integer>();
		int[] pairing = new int[gappedstructure.length];
		for (int i = 0; i < pairing.length; i++) {
			pairing[i] = -1;
		}
		for (int i = 0; i < gappedstructure.length; i++) {
			if (gappedstructure[i] == '(') {
				s.push(i);
			}
			if (gappedstructure[i] == ')') {
				if (s.isEmpty()) {
					System.err.println("Unmatched right paranthesis at pos "
							+ i);
				}
				int leftPos = s.pop();
				pairing[leftPos] = i;
				pairing[i] = leftPos;
			}
		}
		if (!s.isEmpty()) {
			System.err.println("Unmatched left parantheses at pos " + s.pop());
		}
		for (int i = 0; i < gappedsequencein.length(); i++) {
			// turn all positions with gaps to unpaired
			if (pairing[i] != -1 && MatrixTools.isGap(gappedsequence[i])) {
				gappedstructure[i] = '.';
				gappedstructure[pairing[i]] = '.';
			}
		}

		char[] reducedstructure = new char[cnt];
		int pos = 0;

		for (int i = 0; i < gappedsequencein.length(); i++) {
			if (!MatrixTools.isGap(gappedsequence[i])) {
				reducedstructure[pos] = gappedstructure[i];
				pos++;
			}
		}
		//System.out.println(new String(reducedstructure));
		return new String(reducedstructure);
	}

	public static int writeCTFormat(String dir_fil, String filename,
			String extension, List<char[]> columns, char[] structure,
			float[] reliability) throws IOException {

		//System.out.println(dir_fil);

		// Create consensus RNA structure pairing array.
		Stack<Integer> s = new Stack<Integer>();
		int[] pairing = new int[structure.length];
		for (int i = 0; i < pairing.length; i++) {
			pairing[i] = -1;
		}
		for (int i = 0; i < structure.length; i++) {
			if (structure[i] == '(') {
				s.push(i);
			}
			if (structure[i] == ')') {
				if (s.isEmpty()) {
					System.err
							.println("Unmatched right parantheses at position "
									+ i + " in " + new String(structure));
				}
				int leftPos = s.pop();
				pairing[leftPos] = i;
				pairing[i] = leftPos;
			}
		}
		if (!s.isEmpty()) {
			System.err.println("Unmatched left parantheses at pos " + s.pop()
					+ " in " + new String(structure));
		}

		new File(dir_fil).mkdir();

		PrintWriter out;
		String newname = filename.replaceAll("/", "-");
		newname = newname.replaceAll("\\*", "-");

		try {
			out = new PrintWriter(new FileWriter(dir_fil + "/" + newname
					+ extension));
			out.print("\t" + columns.size() + "  ENERGY = "
					+ calcAverage(reliability) + " \n");

			for (int c = 0; c < structure.length; c++) {
				// out.print("\t" + (c+1) + " " +
				// getMaximumConsensus(columns.get(c)) + "     " +
				// c + "     " + (c+2) + "\t" + (pairing[c]+1) + "     " + (c+1)
				// + "\n" );

				out.format("\t %d " + getMaximumConsensus(columns.get(c)).toLowerCase()
						+ "%6d" + "%6d\t%6d%6d\n", (c + 1), c, (c + 2),
						(pairing[c] + 1), (c + 1));
			}
			out.close();
			System.out.println("Written file " + dir_fil + "/" + newname
					+ extension);
			return 0;
		} catch (IOException e) {
			System.err.println("Error writing file: " + dir_fil + "/" + newname
					+ extension);
			throw new IOException(e);
		}
	}

	
	
	public static int writeCTFormat_toStream(List<char[]> columns, int[] pairing,
			float[] reliability, BufferedOutputStream out) throws IOException {
		try {
			out.write(("\t" + columns.size() + "  ENERGY = "
					+ calcAverage(reliability) + " \n").getBytes());

			for (int c = 0; c < pairing.length; c++) {
				// out.print("\t" + (c+1) + " " +
				// getMaximumConsensus(columns.get(c)) + "     " +
				// c + "     " + (c+2) + "\t" + (pairing[c]+1) + "     " + (c+1)
				// + "\n" );

				out.write( String.format("\t %d " + getMaximumConsensus(columns.get(c)).toLowerCase()
						+ "%6d" + "%6d\t%6d%6d\n", (c + 1), c, (c + 2),
						(pairing[c] + 1), (c + 1)).getBytes() );
			}
			out.close();
			System.out.println("Written .ct file ");
			return 0;
		} catch (IOException e) {
			System.err.println("Error writing .ct file!");
			throw new IOException(e);
		}
	}

	
	
	private static double calcAverage(float[] reliability) {
		float sum = 0;
		for (int i = 0; i < reliability.length; i++) {
			sum += reliability[i];
		}
		float average = sum / reliability.length;
		return average;
	}

	private static String getAbsoluteConsensus(char[] column) {
		boolean isa = false, isg = false, isc = false, isu = false;

		for (int i = 0; i < column.length; i++) {
			column[i] = Character.toLowerCase(column[i]);
		}

		for (int i = 0; i < column.length; i++) {
			if (!isa && column[i] == 'a') {
				isa = true;
			} else if (!isc && column[i] == 'c') {
				isc = true;
			} else if (!isg && column[i] == 'g') {
				isg = true;
			} else if (!isu && (column[i] == 't' || column[i] == 'u')) {
				isu = true;
			} else {
			}
		}
		return getSymbol(isa, isg, isc, isu);
	}

	private static String getRelativeConsensus(char[] column) {
		float threshold = 0.5f; // when should variations be ignored;
		// minimum=0.5
		boolean isa = false, isg = false, isc = false, isu = false;
		int ca = 0, cg = 0, cc = 0, cu = 0;

		for (int i = 0; i < column.length; i++) {
			column[i] = Character.toLowerCase(column[i]);
		}

		for (int i = 0; i < column.length; i++) {
			if (column[i] == 'a') {
				ca++;
			} else if (column[i] == 'c') {
				cc++;
			} else if (column[i] == 'g') {
				cg++;
			} else if (column[i] == 't' || column[i] == 'u') {
				cu++;
			} else {
			}
		}

		if ((float) ca / (float) column.length > threshold) {
			isa = true;
		} else if ((float) cc / (float) column.length > threshold) {
			isc = true;
		} else if ((float) cg / (float) column.length > threshold) {
			isg = true;
		} else if ((float) cg / (float) column.length > threshold) {
			isu = true;
		} else {
			return getAbsoluteConsensus(column);
		}

		return getSymbol(isa, isg, isc, isu);
	}

	private static String getMaximumConsensus(char[] column) {
		int ca = 0, cg = 0, cc = 0, cu = 0;

		for (int i = 0; i < column.length; i++) {
			column[i] = Character.toLowerCase(column[i]);
		}

		for (int i = 0; i < column.length; i++) {
			if (column[i] == 'a') {
				ca++;
			} else if (column[i] == 'c') {
				cc++;
			} else if (column[i] == 'g') {
				cg++;
			} else if (column[i] == 't' || column[i] == 'u') {
				cu++;
			} else {
			}
		}

		int maxarg = getMaximumArgument(ca, cg, cc, cu);

		if (maxarg == 1) {
			return getSymbol(true, false, false, false);
		} else if (maxarg == 2) {
			return getSymbol(false, true, false, false);
		} else if (maxarg == 3) {
			return getSymbol(false, false, true, false);
		} else if (maxarg == 4) {
			return getSymbol(false, false, false, true);
		} else {
			return getRelativeConsensus(column);
		}
	}

	private static String getSymbol(boolean isa, boolean isg, boolean isc,
			boolean isu) {

		if (isa && !isg && !isc && !isu) {
			return "A";
		} else if (!isa && !isg && !isc && isu) {
			return "U";
		} else if (!isa && isg && !isc && !isu) {
			return "G";
		} else if (!isa && !isg && isc && !isu) {
			return "C";
		}

		else if (isa && isg && !isc && !isu) {
			return "R";
		} else if (!isa && !isg && isc && isu) {
			return "Y";
		} else if (!isa && isg && isc && !isu) {
			return "S";
		} else if (isa && !isg && !isc && isu) {
			return "W";
		} else if (!isa && isg && !isc && isu) {
			return "K";
		} else if (isa && !isg && isc && !isu) {
			return "M";
		}

		else if (!isa && isg && isc && isu) {
			return "B";
		} else if (isa && isg && !isc && isu) {
			return "D";
		} else if (isa && !isg && isc && isu) {
			return "H";
		} else if (isa && isg && isc && !isu) {
			return "V";
		}

		else {
			return "N";
		}

	}

	private static int getMaximumArgument(int arg1, int arg2, int arg3, int arg4) {
		// returns the order nr. of the argument that is maximum.
		if ((arg1 > arg2) && (arg1 > arg3) && (arg1 > arg4)) {
			return 1;
		} else if ((arg2 > arg1) && (arg2 > arg3) && (arg2 > arg4)) {
			return 2;
		} else if ((arg3 > arg1) && (arg3 > arg2) && (arg3 > arg4)) {
			return 3;
		} else if ((arg4 > arg1) && (arg4 > arg2) && (arg4 > arg3)) {
			return 4;
		} else {
			return 0;
		}
	}

	public static int writeRpToStream(String sequence, int[] pairing,
			float[] reliability, BufferedOutputStream out) throws IOException {
		try {
			for (int c = 0; c < pairing.length; c++) {
				// out.print("\t" + (c+1) + " " +
				// getMaximumConsensus(columns.get(c)) + "     " +
				// c + "     " + (c+2) + "\t" + (pairing[c]+1) + "     " + (c+1)
				// + "\n" );

				out.write( String.format("\t %d " 
						+ "\t%6d\t%5.4f\n", (c + 1), 
						(pairing[c] + 1), reliability[c]).getBytes());
			}
			out.close();
			System.out.println("Written .rp file ");
			return 0;
		} catch (IOException e) {
			System.err.println("Error writing .rp file!");
			throw new IOException(e);
		}
	}

	public static int writeTabbedMatrixD(String dir_fil, String filename,
			String extension, double[][] matrix) throws IOException {

		new File(dir_fil).mkdir();

		PrintWriter out;
		String newname = filename.replaceAll("/", "-");
		newname = newname.replaceAll("\\*", "-");

		try {
			out = new PrintWriter(new FileWriter(dir_fil + "/" + newname
					+ extension));
			int i;
			for (i = 0; i < matrix.length; i++) {
				double[] row = matrix[i];
				boolean dlm = false;
				int j;
				for (j = 0; j < row.length; j++) {
					if (dlm) {
						out.print('\t');
					} else {
						// out.print(i + "\t"); //prints line number
						dlm = true;
					}
					out.format("%7.3e", row[j]);
				}
				out.println();
			}
			out.close();
			System.out.println("Written file " + dir_fil + "/" + newname
					+ extension);
			return 0;
		} catch (IOException e) {
			System.err.println("Error writing file: " + dir_fil + "/" + newname
					+ extension);
			throw new IOException(e);
		}
		
	}
	
}
