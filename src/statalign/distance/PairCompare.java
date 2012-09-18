package statalign.distance;

import java.util.Comparator;

/**
 * Makes it possible to compare two pairs. 
 * The element to the left is the dominant one.
 *  
 * 
 * @author ingolfur
 *
 */

public class PairCompare implements Comparator < Pair <String,Integer> >  {

	@Override
	public int compare(Pair<String,Integer>  o1, Pair<String,Integer>  o2) {
		return o1.getLeft().compareTo(o2.getLeft());
	}
}
