package com.ppfold.algo;

/**
 * Methods for matrix and vector operations. Alphabet is defined here
 * 
 * @author Z.Sukosd
 * @see CYKJob
 * @see Outside
 * @see ExpectationMatrixCalc
 */

public class MatrixTools {

	final static String alphabet = "augcryswkmbdhvn";
	
	public static int[] convertColumn(char[] column_char) {
		int[] result = new int[column_char.length];

		for (int i = 0; i < column_char.length; i++) {
			result[i] = convertChar(column_char[i]);
		}
		return result;
	}

	
	public static double[][] createSingleVectors() {
		final double[][] matrix = new double[alphabet.length()][4];
		for (int i = 0; i < alphabet.length(); i++) {
			final char c = alphabet.charAt(i);
			matrix[i] = createSVector(c);
		}
		return matrix;
	}

	public static double[][][] createDoubleVectors() {
		final double[][][] matrix = new double[alphabet.length()][alphabet
				.length()][16];
		double[] tmplongvector = new double[16];
		final double[][] tmpmatrix = new double[4][4];
		for (int i = 0; i < alphabet.length(); i++) {
			for (int j = 0; j < alphabet.length(); j++) {
				final char a = alphabet.charAt(i);
				final char b = alphabet.charAt(j);
				tmplongvector = createDVector(a, b, tmpmatrix, tmplongvector);
				copyFromTo(tmplongvector, matrix[i][j]);
			}
		}
		return matrix;
	}

	public static double[] createVector(final int n, final double value) {
		final double[] result = new double[n];
		for (int i = 0; i < n; i++) {
			result[i] = value;
		}
		return result;
	}
	
	public static void resetVector(double [] input, double value){
		for(int i = 0; i<input.length; i++){
			input[i] = value;
		} 
	}

	public static double[][] multiply(final double[][] input1,
			final double[][] input2) {
		final double[][] result = new double[input1.length][input2[0].length];
		for (int i = 0; i < input1.length; i++) {
			for (int j = 0; j < input2[0].length; j++) {
				for (int k = 0; k < input1[0].length; k++) {
					result[i][j] += input1[i][k] * input2[k][j];
				}
			}
		}

		return result;
	}

	/**
	 * Calculates matrix product of matrix and vector
	 * Sets result into the first argument. 
	 * Third argument should be initialized to 0. 
	 */		
	public static void multiplyVectorMatrix(final double[] vector, final double[][] matrix,
			 final double[] tmpresult) {
		//avoids the creation of a new object.
		for (int i = 0; i < matrix[0].length; i++) {
			for (int k = 0; k < vector.length; k++) {
				tmpresult[i] += matrix[i][k] * vector[k];
			}
		}
		copyFromTo(tmpresult, vector);
	}
	
	/**
	 * Calculates matrix product of matrix and vector
	 * Sets result into the second argument. 
	 * Third argument should be initialized to 0. 
	 */	
	public static void multiplyMatrixVector(final double[][] input1,
			final double[] input2, final double[] result) {
		//avoids the creation of a new object.
		for (int i = 0; i < input1.length; i++) {
			for (int k = 0; k < input1[0].length; k++) {
				result[i] += input1[i][k] * input2[k];
			}
		}
		copyFromTo(result, input2);
	}

	public static double[] multiplyMatrixVector(final double[][] input1,
			final double[] input2) {
		//creates a new object
		double result[] = new double[input1.length];
		for (int i = 0; i < input1.length; i++) {
			for (int k = 0; k < input1[0].length; k++) {
				result[i] += input1[i][k] * input2[k];
			}
		}
		return result;
	}
	
	
	public static void copyFromTo(final double[] from, final double[] to) {
		for (int i = 0; i < from.length; i++) {
			to[i] = from[i];
		}
	}

	public static double[][] multiplyVectorVector(final double[] input1,
			final double[] input2, final double[][] result) {
		for (int i = 0; i < input1.length; i++) {
			for (int j = 0; j < input2.length; j++) {
				result[i][j] = 0;
			}
		}

		for (int i = 0; i < input1.length; i++) {
			for (int k = 0; k < input2.length; k++) {
				result[i][k] += input1[i] * input2[k];
			}
		}
		return result;
	}

	public static void multiplySeries(final double[] result,
			final double[] input2) {
		for (int i = 0; i < result.length; i++) {
			result[i] = result[i] * input2[i];
		}
	}


	public static double scalarProduct(final double[] input1,
			final double[] input2) {
		double result = 0;

		for (int i = 0; i < input1.length; i++) {
			result += input1[i] * input2[i];
		}
		return result;
	}

	public static double[][] expTD(final double[][] matrix, final double t) {
		// takes diagonal input, returns diagonal output
		final double[][] result = new double[matrix.length][matrix.length];
		for (int i = 0; i < matrix.length; i++) {
			result[i][i] = Math.exp(matrix[i][i] * t);
		}
		return result;
	}

	public static double[][] expRT(final double[][] D, final double t,
			final double[][] V, final double[][] V1) {
		double[][] result = new double[V.length][V1[0].length];
		result = expTD(D, t);
		result = multiply(result, V1);
		result = multiply(V, result);
		return result;
	}

	public static void print(final double[][] input) {
		for (int r = 0; r < input.length; r++) {
			for (int c = 0; c < input[r].length; c++) {
				System.out.print("\t" + input[r][c]);
			}
			System.out.println();
		}
	}

	public static void print(final int[][] input) {
		for (int r = 0; r < input.length; r++) {
			for (int c = 0; c < input[r].length; c++) {
				System.out.print("\t" + input[r][c]);
			}
			System.out.println();
		}
	}
	
	public static void prints(final double[] input) {
		System.out.print("[");
		for (int r = 0; r < input.length; r++) {
			System.out.print(" " + input[r]);
		}
		System.out.println("]");
	}


	public static void print(final int[] input) {
		for (int r = 0; r < input.length; r++) {
			System.out.print(" " + input[r]);
		}
		System.out.println();
	}
	
	public static double[] serializeMatrix(final double[][] input,
			final double[] result) {
		int cnt = 0;
		for (int i = 0; i < input.length; i++) {
			for (int j = 0; j < input[0].length; j++) {
				result[cnt] = input[i][j];
				cnt++;
			}
		}
		return result;
	}
	
	public static char[] createNts(char nt) {
		nt = Character.toLowerCase(nt);
		char[] nts;
		if (nt == 'a') {
			nts = new char[1];
			nts[0] = 'a';
		} else if (nt == 'u' || nt == 't') {
			nts = new char[1];
			nts[0] = 'u';
		} else if (nt == 'g') {
			nts = new char[1];
			nts[0] = 'g';
		} else if (nt == 'c') {
			nts = new char[1];
			nts[0] = 'c';
		} else if (nt == 'r') {
			nts = new char[2];
			nts[0] = 'a';
			nts[1] = 'g';
		} else if (nt == 'y') {
			nts = new char[2];
			nts[0] = 'c';
			nts[1] = 't';
		} else if (nt == 's') {
			nts = new char[2];
			nts[0] = 'g';
			nts[1] = 'c';
		} else if (nt == 'w') {
			nts = new char[2];
			nts[0] = 'a';
			nts[1] = 't';
		} else if (nt == 'k') {
			nts = new char[2];
			nts[0] = 'g';
			nts[1] = 't';
		} else if (nt == 'm') {
			nts = new char[2];
			nts[0] = 'a';
			nts[1] = 'c';
		} else if (nt == 'b') {
			nts = new char[3];
			nts[0] = 'c';
			nts[1] = 'g';
			nts[2] = 't';
		} else if (nt == 'd') {
			nts = new char[3];
			nts[0] = 'a';
			nts[1] = 'g';
			nts[2] = 't';
		} else if (nt == 'h') {
			nts = new char[3];
			nts[0] = 'a';
			nts[1] = 'c';
			nts[2] = 't';
		} else if (nt == 'v') {
			nts = new char[3];
			nts[0] = 'a';
			nts[1] = 'c';
			nts[2] = 'g';
		} else if (nt == 'n' || MatrixTools.isGap(nt)) {
			nts = new char[4];
			nts[0] = 'a';
			nts[1] = 'c';
			nts[2] = 'g';
			nts[3] = 'u';
		} else {
			System.out.println("Illegal character: " + nt);
			return null;
		}
		return nts;
	}

	public static double[] createNtVector(final char nt) {
		final double[] vector = createVector(4, 0);
		final char[] nts = createNts(nt);

		for (final char n : nts) {

			if (n == 'a') {
				vector[0] = 1;
			}
			if (n == 'u') {
				vector[1] = 1;
			}
			if (n == 'g') {
				vector[2] = 1;
			}
			if (n == 'c') {
				vector[3] = 1;
			}
		}
		return vector;
	}

	
	public static double[] createSVector(final char nt) {
		final double[] vector = createVector(4, 0.01); //uncertainty in nucleotide
		//final double[] vector = createVector(4, 0.00);

		final char[] nts = createNts(nt);

		for (final char n : nts) {

			if (n == 'a') {
				vector[0] = 1;
			}
			if (n == 'u') {
				vector[1] = 1;
			}
			if (n == 'g') {
				vector[2] = 1;
			}
			if (n == 'c') {
				vector[3] = 1;
			}
		}
		return vector;
	}

	public static double[] createDVector(final char nt1, final char nt2,
			final double[][] tmpmatrix, final double[] tmplongvector) {

		final double[] vector1 = createSVector(nt1);
		final double[] vector2 = createSVector(nt2);
		return serializeMatrix(
				multiplyVectorVector(vector1, vector2, tmpmatrix),
				tmplongvector);
	}

	static int convertChar(char c) {
		// REMEMBER: Alphabets in the whole file must match!!
		// Takes in a character, returns a number that represents its position in the alphabet
		c = Character.toLowerCase(c);
		
		//convert to u
		if(c=='t'){c='u';}
		
		//how to handle gaps? 
		if(isGap(c)){c='n';}
		
		for(int i = 0; i<alphabet.length(); i++){
			if(c==alphabet.charAt(i)){
				return i;
			}
		}
		return -1;
	}


	public static boolean isGap(char c) {
		if(c == '-' || c == '.'){
			return true;
		}
		else{
			return false;
		}
	}

}
