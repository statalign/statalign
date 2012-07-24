package statalign.distance;

import java.util.ArrayList;
import java.util.Comparator;

public class ComplexPair implements Comparator < Pair < ArrayList <String>,String> >  {
	@Override
	public int compare(Pair<ArrayList<String>,String>  o1, Pair<ArrayList<String>,String>  o2) {
		return o1.getRight().compareTo(o2.getRight());
	}
}
