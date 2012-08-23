package statalign.distance;

import java.util.ArrayList;
import java.util.Comparator;

/**
 * A class that makes it possible to compare two pairs where the first element of
 * each pair is an  ArrayList of Strings (alignment) and the second element is a String
 * 
 * The String (the second element of the pair) determines which is bigger
 * 
 * @author Ingolfur
 * 
 *
 */

public class ComplexPair implements Comparator < Pair < ArrayList <String>,String> >  {
	@Override
	public int compare(Pair<ArrayList<String>,String>  o1, Pair<ArrayList<String>,String>  o2) {
		return o1.getRight().compareTo(o2.getRight());
	}
}
